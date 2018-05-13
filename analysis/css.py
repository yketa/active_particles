"""
Module css calculates shear strain and displacement vorticity as well as
their correlations.

Files are saved according to active_particles.naming.Css (shear strain) and
active_particles.naming.Ccc (displacement vorticity) standards. Data in these
files can be plotted with active_particles.analysis.css_plot.

Environment parameters
----------------------
DATA_DIRECTORY : string
	Data directory.
	DEFAULT: current working directory
PARAMETERS_FILE : string
	Simulation parameters file.
	DEFAULT: DATA_DIRECTORY/active_particles.naming.parameters_file
WRAPPED_FILE : string
	Wrapped trajectory file. (.gsd)
	DEFAULT: DATA_DIRECTORY/active_particles.naming.wrapped_trajectory_file
UNWRAPPED_FILE : string
	Unwrapped trajectory file. (.dat)
	NOTE: .dat files defined with active_particles.dat
	DEFAULT: DATA_DIRECTORY/active_particles.naming.unwrapped_trajectory_file
INITIAL_FRAME : int
	Frame to consider as initial.
	NOTE: INITIAL_FRAME < 0 will be interpreted as the initial frame being
	the middle frame of the simulation.
	DEFAULT: -1
TIME : int
	Lag time for displacement.
	NOTE: TIME < 0 will be interpreted as a lag time corresponding to the total
	number of simulation frames - INITIAL_FRAME + TIME.
	DEFAULT: -1
INTERVAL_MAXIMUM : int
	Maximum number of intervals of length dt considered in correlations
	calculations.
	DEFAULT: 1
R_CUT : float
	Cut-off radius for coarse graining function.
	NOTE: Value of R_CUT is then multiplied by average particle diameter from
	simulation parameters file.
	DEFAULT: 2
SIGMA : float
	Length scale of the spatial extent of the coarse graining function.
	DEFAULT: R_CUT
N_CASES : int
	Number of boxes in each direction to compute the shear strain and
	displacement vorticity grid.
	DEFAULT: smallest integer value greater than or equal to the square root of
		the number of particles from the simulation parameters file.

Environment parameters to set when plotting data (SHOW=True)
------------------------------------------------
R_MAX : float
	Half size of the box showed for 2D correlation.
	NOTE: R_MAX < 0 will be interpreted as the box shown being the actual
	simulation box.
	DEFAULT: 20
DISPLAY_GRID : int
	Index of map in list of variable maps to display.
	DEFAULT : 0

Output
------
> Prints execution time.
> Saves the list of computed shear strain maps and the averaged shear strain
correlation according to active_particles.naming.Css standards in
DATA_DIRECTORY.
> Saves the list of computed displacement vorticity maps and the averaged
displacement vorticity correlations according to active_particles.naming.Ccc
standards in DATA_DIRECTORY.
> Plots data map and correlation for both shear strain and displacement
vorticity if SHOW=True is set in the environment.
	> These figures are saved in DATA_DIRECTORY if SAVE=True is set in the
	environment.
"""

import active_particles.naming as naming

from active_particles.init import get_env
from active_particles.exponents import float_to_letters
from active_particles.dat import Dat

from active_particles.analysis.neighbours import NeighboursGrid
from active_particles.analysis.correlations import corField2D_scalar_average
from active_particles.analysis.coarse_graining import CoarseGraining

from os import getcwd as current_dir

from math import ceil

import gsd
import gsd.hoomd
import gsd.pygsd

import numpy as np

import pickle

from operator import itemgetter

from collections import OrderedDict

from datetime import datetime

import matplotlib as mpl
if get_env('SHOW', default=False, vartype=bool):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
cmap = plt.cm.jet

def define_variables():
	"""
	Defines variables which are needed both to compute shear strain and
	displacement vorticity correlations or to plot already computed shear
	strain and displacement vorticity correlations from
	active_particles.analysis.Css_plot.

	NOTE: Variables are defined globally inside the Css module.
	"""

	global data_dir ; data_dir = get_env('DATA_DIRECTORY',
		default=current_dir())	# data directory

	global dt ; dt = get_env('TIME', default=-1, vartype=int)	# lag time for displacement

	global init_frame ; init_frame = get_env('INITIAL_FRAME', default=-1,
		vartype=int)	# frame to consider as initial
	global int_max ; int_max = get_env('INTERVAL_MAXIMUM', default=1,
		vartype=int)	# maximum number of intervals of length dt considered in correlations calculations

	global r_max ; r_max = get_env('R_MAX', default=20, vartype=float)	# half size of the box showed for 2D correlation

	parameters_file = get_env('PARAMETERS_FILE',
		default=data_dir + '/' + naming.parameters_file)			# simulation parameters file
	with open(parameters_file, 'rb') as param_file:
		global parameters ; parameters = pickle.load(param_file)	# parameters hash table

	global r_cut ; r_cut = parameters['a']*get_env('R_CUT', default=2,
		vartype=float)														# cut-off radius for coarse graining function
	global sigma ; sigma = get_env('SIGMA', default=r_cut, vartype=float)	# length scale of the spatial extent of the coarse graining function

	r_max = parameters['box_size']/2 if r_max < 0 else r_max	# half size of the box showed for 2D correlation

	global prep_frames ; prep_frames = ceil(parameters['prep_steps']
		/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

	global Ncases ; Ncases = get_env('N_CASES',
		default=ceil(np.sqrt(parameters['N'])), vartype=int)	# number of boxes in each direction to compute the shear strain and displacement vorticity grid

	global Nentries ; Nentries =\
		parameters['N_steps']//parameters['period_dump'] 			# number of time snapshots in velocity and position files
	init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame
	global Nframes ; Nframes = Nentries - init_frame 				# number of frames available for the calculation

	dt = Nframes + dt if dt <= 0 else dt	# length of the interval of time for which displacements are calculated

	attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'init_frame': init_frame, 'dt': dt,
		'int_max': int_max, 'Ncases': Ncases, 'r_cut': r_cut, 'sigma': sigma}	# attributes displayed in filenames
	naming_Css = naming.Css()													# Css naming object
	global Css_filename ; Css_filename, = naming_Css.filename(**attributes)
	naming_Ccc = naming.Ccc()													# Ccc naming object
	global Ccc_filename ; Ccc_filename, = naming_Ccc.filename(**attributes)

	global display_grid ; display_grid = get_env('DISPLAY_GRID', default=0,
		vartype=int)	# index of map in list of variable maps to display

def strain_vorticity(point, time, dt, positions, u_traj, sigma, r_cut,
	box_size, neighbours_grid):
	"""
	From coarse-grained displacement field u, calculates
	> linearised shear strain, \\epsilon_{xy} =
		1/2 (\\frac{\\partial u_y}{\\partial x}
		     + \\frac{\\partial u_x}{\\partial y})
	> and displacement vorticity, \\omega =
		\\frac{\\partial u_y}{\\partial x} - \\frac{\\partial u_x}{\\partial y}

	Resorts to neighbours grids to speed up calculations.

	Parameters
	----------
	point : array like
		Coordinates of point at which to calculate shear strain and
		displacement vorticity.
	time : int
		Frame at which shear strain will be calculated.
	dt : int
		Length of the interval of time for which the displacements are
		calculated.
	positions : (N, 2) shaped array like
		Array of wrapped particle positions.
	u_traj : dat.Dat
		Unwrapped trajectory object.
	sigma : float
		Length scale of the spatial extent of the coarse graining function.
	r_cut : float
		Cut-off radius for coarse graining function.
		Also cut-off radius with which neighbours grid has been computed.
	box_size : float
		Length of the system's square box.
	neighbours_grid : active_particles.analysis.neighbours.NeighboursGrid
		Neighbours grid.

	Returns
	-------
	strain : float
		Linearised shear strain at point.
	vorticity : float
		Displacement vorticity at point.
	"""

	wrcut = neighbours_grid.get_neighbours(point)	# particles indexes within r_cut of point
	Nwrcut = len(wrcut)								# number of neighbouring particles
	if Nwrcut == 0:									# if there is no particles within r_cut of point
		return 0, 0

	pos_wrcut = (np.array(itemgetter(*wrcut)(positions), ndmin=2)
		- point + box_size/2)%box_size - box_size/2 # position at time time (with boundary conditions) with point as the centre of the frame
	pos0_wrcut = u_traj.position(time, *wrcut)		# positions at time time (without periodic boundary conditions)
	pos1_wrcut = u_traj.position(time + dt, *wrcut)	# positions at time time + dt (without periodic boundary conditions)
	dis_wrcut = pos1_wrcut - pos0_wrcut				# displacements of the particles between time and time + dt

	coarse_graining = CoarseGraining(np.sqrt(np.sum(pos_wrcut**2, axis=-1)),
		sigma, r_cut)														# coarse graining object
	rho = coarse_graining.average(np.full(len(pos_wrcut), fill_value=1))	# averaged density
	if rho == 0: return 0, 0												# if there is no particles within r_cut

	Ax, Ay, Aux, Auy, Auxy, Auyx = tuple(map(coarse_graining.average,
		[pos_wrcut[:, 0], pos_wrcut[:, 1], dis_wrcut[:, 0], dis_wrcut[:, 1],
		dis_wrcut[:, 0]*pos_wrcut[:, 1], dis_wrcut[:, 1]*pos_wrcut[:, 0]]))	# coarse grained x, y, u_x, u_y, u_x * y, u_y * x
	strain = 0.5*((Ay*Aux + Ax*Auy)/((rho*sigma)**2)
		- (Auxy + Auyx)/(rho*(sigma**2)))									# linearised shear strain
	vorticity = (Ax*Auy - Ay*Aux)/((rho*sigma)**2)\
		- (Auyx - Auxy)/(rho*(sigma**2))									# displacement vorticity

	return strain, vorticity

def strain_vorticity_grid(box_size, Ncases, grid_points, time, dt,
	prep_frames, w_traj, u_traj, sigma, r_cut):
	"""
	Calculates grids of (linearised) shear strain and displacement vorticity
	from coarse-grained displacement field.

	Resorts to neighbours grids to speed up calculations.

	Parameters
	----------
	box_size : float
		Length of the system's square box.
	Ncases : int
		Number of boxes in each direction to compute the shear strain and
		displacement vorticity grid.
	grid_points : array like of coordinates
		Grid points at which shear strain will be evaluated.
	time : int
		Frame at which shear strain will be calculated.
	dt : int
		Length of the interval of time for which the displacements are
		calculated.
	prep_frames : int
		Number of preparation frames.
	w_traj : gsd.hoomd.HOOMDTrajectory
		Wrapped trajectory variable.
	u_traj : active_particles.dat.Dat
		Unwrapped trajectory object.
	sigma : float
		Length scale of the spatial extent of the coarse graining function.
	r_cut : float
		Cut-off radius for coarse graining function.

	Returns
	-------
	Sgrid : 2D array like
		Shear strain grid.
	Cgrid : 2D array like
		Displacement vorticity grid.

	Output
	------
	Prints neighbours grid computation time.
	"""

	# NEIGHBOURS GRID

	startTime0 = datetime.now()	# start time for neighbours grid computation

	positions = w_traj[int(prep_frames + time +
		dt*get_env('ENDPOINT', default=False, vartype=bool)
		)].particles.position[:, :2]								# array of wrapped particle positions
	neighbours_grid = NeighboursGrid(positions, box_size, r_cut)	# neighbours grid

	print("Neighbours grid computation time (time = %e): %s" %
		(time, datetime.now() - startTime0))	# neighbours grid computation time

	# SHEAR STRAIN AND DISPLACEMENT VORTICITY GRIDS CALCULATION

	Sgrid, Cgrid = tuple(np.transpose(list(map(
		lambda point: strain_vorticity(point, time, dt, positions,
		u_traj, sigma, r_cut, box_size, neighbours_grid),
		grid_points))))	# shear strain and displacement vorticity lists

	correct_grid = lambda Grid: np.transpose(
		np.reshape(Grid, (Ncases, Ncases)))[::-1]	# get grids with the same orientation as positions
	return correct_grid(Sgrid), correct_grid(Cgrid)	# shear strain and displacement vorticity grids

def plot(Grid, Corr, var):
	"""
	Plots variable grid and correlations.

	Parameters
	----------
	Grid : list of 2D array-like
		List of variable grids.
	Corr : 2D array-like
		Variable correlation grid.
	var : string
		Name of variable.
	"""

	fig, ax = plt.subplots(1, 2)	# 1 x 2 figure

	fig.set_size_inches(16, 16)		# figure size
	fig.subplots_adjust(wspace=0.4)	# width space
	fig.subplots_adjust(hspace=0.3)	# height space

	fig.suptitle(r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
		% (parameters['N'], parameters['density'], parameters['vzero'],
		parameters['dr']) + '\n' +
		r'$S_{init}=%.2e, \Delta t=%.2e$' % (init_frame,
		dt*parameters['period_dump']*parameters['time_step']) +
		r'$, S_{max}=%.2e, N_{cases}=%.2e$' % (int_max, Ncases) +
		r'$, r_{cut}=%.2e, \sigma=%.2e$' % (r_cut, sigma))

	# VARIABLE GRID

	S2D = Grid[display_grid]	# displaying only the first computed variable grid

	Smin = -2*np.std(S2D)
	Smax = 2*np.std(S2D)

	SvNorm = colors.Normalize(vmin=Smin, vmax=Smax)
	SscalarMap = cmx.ScalarMappable(norm=SvNorm, cmap=cmap)

	ax[0].imshow(S2D, cmap=cmap, norm=SvNorm,
		extent=[-parameters['box_size']/2, parameters['box_size']/2,
			-parameters['box_size']/2, parameters['box_size']/2])

	ax[0].set_xlabel(r'$x$')
	ax[0].set_ylabel(r'$y$')
	ax[0].set_title('2D ' + r'$%s$' % var)

	divider0 = make_axes_locatable(ax[0])
	cax0 = divider0.append_axes("right", size="5%", pad=0.05)
	cb0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap, norm=SvNorm,
		orientation='vertical')
	cb0.set_label(r'$%s$' % var, labelpad=20, rotation=270)

	# 2D CORRELATIONS

	C = 'C_{%s%s}' % (var, var)

	Cmin = np.min(Corr)
	Cmax = np.max(Corr)

	CvNorm = colors.Normalize(vmin=Cmin, vmax=Cmax)
	CscalarMap = cmx.ScalarMappable(norm=CvNorm, cmap=cmap)

	r_max_cases = int(r_max*(Ncases/parameters['box_size']))
	C2D_display = np.roll(np.roll(Corr, int(Ncases/2), axis=0),
		int(Ncases/2), axis=1)[int(Ncases/2) - r_max_cases:
		int(Ncases/2) + r_max_cases + 1, int(Ncases/2) - r_max_cases:
		int(Ncases/2) + r_max_cases + 1]	# part of variable correlations to display

	ax[1].imshow(C2D_display, cmap=cmap, norm=CvNorm,
		extent=[-r_max, r_max, -r_max, r_max])

	ax[1].set_xlabel(r'$x$')
	ax[1].set_ylabel(r'$y$')
	ax[1].set_title('2D ' + r'$%s$' % C)

	divider1 = make_axes_locatable(ax[1])
	cax1 = divider1.append_axes("right", size="5%", pad=0.05)
	cb1 = mpl.colorbar.ColorbarBase(cax1, cmap=cmap, norm=CvNorm,
		orientation='vertical')
	cb1.set_label(r'$%s$' % C, labelpad=20, rotation=270)

if __name__ == '__main__':	# executing as script

	startTime = datetime.now()

	# VARIABLE DEFINITIONS

	define_variables()	# define some of variables needed to compute shear strain and displacement vorticity correlations

	wrap_file_name = get_env('WRAPPED_FILE',
		default=data_dir + '/' + naming.wrapped_trajectory_file)	# wrapped trajectory file (.gsd)
	unwrap_file_name = get_env('UNWRAPPED_FILE',
		default=data_dir + '/' + naming.unwrapped_trajectory_file)	# unwrapped trajectory file (.dat)

	grid_points = np.array([(x, y) for x in\
		np.linspace(- parameters['box_size']*(1 - 1./Ncases)/2,
		parameters['box_size']*(1 - 1./Ncases)/2, Ncases, endpoint=True)\
		for y in\
		np.linspace(- parameters['box_size']*(1 - 1./Ncases)/2,
		parameters['box_size']*(1 - 1./Ncases)/2, Ncases, endpoint=True)
		])	# grid points at which shear strain will be evaluated

	times = np.array(list(OrderedDict.fromkeys(map(
		lambda x: int(x),
		np.linspace(init_frame, Nentries - dt - 1, int_max)
		))))	# frames at which shear strain will be calculated

	# SHEAR STRAIN, DISPLACEMENT VORTICITY AND THEIR CORRELATIONS

	with gsd.pygsd.GSDFile(open(wrap_file_name, 'rb')) as wrap_file,\
		open(unwrap_file_name, 'rb') as unwrap_file:	# opens wrapped and unwrapped trajectory files

		w_traj = gsd.hoomd.HOOMDTrajectory(wrap_file);	# wrapped trajectory object
		u_traj = Dat(unwrap_file, parameters['N'])		# unwrapped trajectory object
		Sgrid, Cgrid = tuple(np.transpose(list(map(lambda time:
			strain_vorticity_grid(parameters['box_size'], Ncases, grid_points,
			time, dt, prep_frames, w_traj, u_traj, sigma, r_cut)
			, times)), (1, 0, 2, 3)))					# lists of shear strain and displacement vorticity correlations

	Css2D, Ccc2D = tuple(map(corField2D_scalar_average, [Sgrid, Cgrid]))	# shear strain and displacement vorticity fields correlations

	# SAVING

	with open(data_dir + '/' + Css_filename, 'wb') as Css_dump_file,\
		open(data_dir + '/' + Ccc_filename, 'wb') as Ccc_dump_file:
		pickle.dump([Sgrid, Css2D], Css_dump_file)
		pickle.dump([Cgrid, Ccc2D], Ccc_dump_file)

	# EXECUTION TIME

	print("Execution time: %s" % (datetime.now() - startTime))

	# PLOT

	if get_env('SHOW', default=False, vartype=bool):

		plot(Sgrid, Css2D, '\epsilon_{xy}')	# plotting shear strain map and correlation
		plot(Cgrid, Ccc2D, '\omega')		# plotting displacement vorticity map and correlation

		plt.show()
