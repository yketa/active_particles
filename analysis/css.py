"""
Module css calculates or plots shear strain and displacement vorticity as well
as their correlations.

Files are saved according to active_particles.naming.Css (shear strain) and
active_particles.naming.Ccc (displacement vorticity) standards.

Environment modes
-----------------
COMPUTE : bool
	Compute shear strain and displacement vorticity.
	DEFAULT: False
PLOT : bool
	Plot saved shear strain and displacement vorticity as well as their
	correlations.
	DEFAULT: False
SHOW [COMPUTE or PLOT mode] : bool
	Show graphs.
	DEFAULT: False
SAVE [COMPUTE or PLOT mode] : bool
	Save graphs.
	DEFAULT: False

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
	DEFAULT: active_particles.analysis.css._r_cut
SIGMA : float
	Length scale of the spatial extent of the coarse graining function.
	DEFAULT: R_CUT
N_CASES : int
	Number of boxes in each direction to compute the shear strain and
	displacement vorticity grid.
	DEFAULT: smallest integer value greater than or equal to the square root of
	         the number of particles from the simulation parameters file.
BOX_SIZE : float
	Size of the square box to consider.
	DEFAULT: simulation box size
X_ZERO : float
	1st coordinate of the centre of the square box to consider.
	DEFAULT: 0
Y_ZERO : float
	2nd coordinate of the centre of the square box to consider.
	DEFAULT: 0
R_MAX [PLOT or SHOW mode] : float
	Half size of the box showed for 2D correlation.
	NOTE: R_MAX < 0 will be interpreted as the box shown being the actual
	      simulation box.
	DEFAULT: active_particles.analysis.css._r_max
DISPLAY_GRID [COMPUTE mode] : int
	Index of map in list of variable maps to display.
	DEFAULT : 0

Output
------
[COMPUTE MODE]
> Prints neigbours grid computation time and execution time.
> Saves a computed map of shear strain and the averaged shear strain
correlation according to active_particles.naming.Css standards in
DATA_DIRECTORY.
> Saves a computed map of displacement vorticity and the averaged displacement
vorticity correlations according to active_particles.naming.Ccc standards
in DATA_DIRECTORY.
[SHOW or PLOT mode]
> Plots data map and correlation for both shear strain and displacement
vorticity.
[SAVE mode]
> Saves data map and correlation figures in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output
from active_particles.dat import Dat, Gsd
from active_particles.maths import relative_positions

from active_particles.analysis.neighbours import NeighboursGrid
from active_particles.analysis.correlations import corField2D_scalar_average
from active_particles.analysis.correlations import CorGrid
from active_particles.analysis.coarse_graining import GaussianCG,\
	CoarseGraining

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

from math import ceil

import numpy as np

import pickle

from operator import itemgetter

from collections import OrderedDict

from datetime import datetime

import matplotlib as mpl
if not(get_env('SHOW', default=False, vartype=bool)):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable

from active_particles.plot.mpl_tools import GridCircle

def strain_vorticity(point, time, dt, positions, displacements,
	sigma, r_cut, box_size, neighbours_grid):
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
	displacements : (N, 2) shaped array like
		Array of particle displacements.
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

	pos_wrcut = relative_positions(
		np.array(itemgetter(*wrcut)(positions), ndmin=2), point, box_size)	# position at time time (with boundary conditions) with point as the centre of the frame
	dis_wrcut = np.array(itemgetter(*wrcut)(displacements), ndmin=2)		# displacements of the particles between time and time + dt

	coarse_graining = CoarseGraining(GaussianCG(sigma, r_cut).factors,
		pos_wrcut)															# coarse graining object
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

def strain_vorticity_grid(box_size, Ncases, grid_points, time, dt, w_traj,
	u_traj, sigma, r_cut):
	"""
	Calculates grids of (linearised) shear strain and displacement vorticity
	from coarse-grained displacement field.

	Resorts to neighbours grids to speed up calculations.

	Parameters
	----------
	box_size : float
		Length of the considered system's square box.
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
	w_traj : active_particles.dat.Gsd
		Wrapped trajectory object.
	u_traj : active_particles.dat.Dat
		Unwrapped trajectory object.
	sigma : float
		Length scale of the spatial extent of the coarse graining function.
	r_cut : float
		Cut-off radius for coarse graining function.

	Returns
	-------
	sgrid : 2D array like
		Shear strain grid.
	cgrid : 2D array like
		Displacement vorticity grid.

	Output
	------
	Prints neighbours grid computation time.
	"""

	# NEIGHBOURS GRID

	startTime0 = datetime.now()	# start time for neighbours grid computation

	positions = w_traj.position(time +
		dt*get_env('ENDPOINT', default=False, vartype=bool))		# array of wrapped particle positions
	neighbours_grid = NeighboursGrid(positions, box_size, r_cut)	# neighbours grid

	displacements = u_traj.displacement(time, time + dt)	# array of particle displacements

	print("Neighbours grid computation time (time = %e): %s" %
		(time, datetime.now() - startTime0))	# neighbours grid computation time

	# SHEAR STRAIN AND DISPLACEMENT VORTICITY GRIDS CALCULATION

	sgrid, cgrid = tuple(np.transpose(list(map(
		lambda point: strain_vorticity(point, time, dt, positions,
		displacements, sigma, r_cut, box_size, neighbours_grid),
		grid_points))))	# shear strain and displacement vorticity lists

	correct_grid = lambda grid: np.transpose(
		np.reshape(grid, (Ncases, Ncases)))[::-1]	# get grids with the same orientation as positions
	return correct_grid(sgrid), correct_grid(cgrid)	# shear strain and displacement vorticity grids

def plot(grid, corr, box_size, var, naming_standard):
	"""
	Plots variable grid and correlations.

	Parameters
	----------
	grid : 2D array-like
		Variable grid.
	corr : 2D array-like
		Variable correlation grid.
	box_size : float
		Length of the box in one dimension.
	var : string
		Name of variable.
	naming_standard : active_particles.naming standard
		Standard naming object.

	Returns
	-------
	fig : matplotlib figure
		Figure.
	ax : matplotlib axis
		Axis.
	fl : active_particles.plot.mpl_tools.GridCircle
		Grid circle object.
	"""

	# GRID AND CORRELATION FIGURE

	cmap = plt.cm.jet

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

	Smin = -2*np.std(grid)
	Smax = 2*np.std(grid)

	SvNorm = colors.Normalize(vmin=Smin, vmax=Smax)
	SscalarMap = cmx.ScalarMappable(norm=SvNorm, cmap=cmap)

	ax[0].imshow(grid, cmap=cmap, norm=SvNorm,
		extent=[-box_size/2, box_size/2, -box_size/2, box_size/2])

	ax[0].set_xlabel(r'$x$')
	ax[0].set_ylabel(r'$y$')
	ax[0].set_title('2D ' + r'$%s$' % var)

	divider0 = make_axes_locatable(ax[0])
	cax0 = divider0.append_axes('right', size='5%', pad=0.05)
	cb0 = mpl.colorbar.ColorbarBase(cax0, cmap=cmap, norm=SvNorm,
		orientation='vertical')
	cb0.set_label(r'$%s$' % var, labelpad=20, rotation=270)

	# 2D CORRELATIONS

	C = 'C_{%s%s}' % (var, var)

	Cmin = np.min(corr)
	Cmax = np.max(corr)

	CvNorm = colors.Normalize(vmin=Cmin, vmax=Cmax)

	cgrid = CorGrid(corr, box_size, display_size=2*r_max)

	ax[1].imshow(cgrid.display_grid.grid, cmap=cmap, norm=CvNorm,
		extent=[-r_max, r_max, -r_max, r_max])

	ax[1].set_xlabel(r'$x$')
	ax[1].set_ylabel(r'$y$')
	ax[1].set_title('2D ' + r'$%s$' % C)

	divider1 = make_axes_locatable(ax[1])
	cax1 = divider1.append_axes("right", size="5%", pad=0.05)
	cb1 = mpl.colorbar.ColorbarBase(cax1, cmap=cmap, norm=CvNorm,
		orientation='vertical')
	cb1.set_label(r'$%s$' % C, labelpad=20, rotation=270)

	# SAVING

	if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
		image_name, = naming_standard.image().filename(**attributes)
		fig.savefig(joinpath(data_dir, image_name))

	# GRID CIRCLE FIGURE

	gridcircle = GridCircle(cgrid.display_grid.grid,
		extent=[-r_max, r_max, -r_max, r_max])
	fig_gc, (ax_grid, ax_plot), cb_gc = gridcircle.get_fig_ax_cmap()

	fig_gc.set_size_inches(16, 16)		# figure size
	fig_gc.subplots_adjust(wspace=0.4)	# width space
	fig_gc.subplots_adjust(hspace=0.3)	# height space

	fig_gc.suptitle(r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
		% (parameters['N'], parameters['density'], parameters['vzero'],
		parameters['dr']) + '\n' +
		r'$S_{init}=%.2e, \Delta t=%.2e$' % (init_frame,
		dt*parameters['period_dump']*parameters['time_step']) +
		r'$, S_{max}=%.2e, N_{cases}=%.2e$' % (int_max, Ncases) +
		r'$, r_{cut}=%.2e, \sigma=%.2e$' % (r_cut, sigma))

	ax_grid.set_xlabel(r'$x$')
	ax_grid.set_ylabel(r'$y$')
	ax_grid.set_title('2D ' + r'$%s$' % C)
	cb_gc.set_label(r'$%s$' % C, labelpad=20, rotation=270)

	ax_plot.set_xlabel(r'$\theta$')
	ax_plot.set_ylabel(r'$%s(r, \theta)$' % C)

	# RETURN

	return fig, ax, gridcircle

# DEFAULT VARIABLES

_r_cut = 2	# default cut-off radius for coarse graining function
_r_max = 20	# default half size of the box showed for 2D correlation

if __name__ == '__main__':	# executing as script

	# VARIABLE DEFINITIONS

	data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

	dt = get_env('TIME', default=-1, vartype=int)	# lag time for displacement

	init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
	int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of intervals of length dt considered in correlations calculations

	parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
	with open(parameters_file, 'rb') as param_file:
		parameters = pickle.load(param_file)				# parameters hash table

	r_cut = parameters['a']*get_env('R_CUT', default=_r_cut, vartype=float)	# cut-off radius for coarse graining function
	sigma = get_env('SIGMA', default=r_cut, vartype=float)			# length scale of the spatial extent of the coarse graining function

	box_size = get_env('BOX_SIZE', default=parameters['box_size'],
		vartype=float)									# size of the square box to consider
	centre = (get_env('X_ZERO', default=0, vartype=float),
		get_env('Y_ZERO', default=0, vartype=float))	# centre of the box

	prep_frames = ceil(parameters['prep_steps']/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

	Ncases = get_env('N_CASES', default=ceil(np.sqrt(parameters['N'])),
		vartype=int)	# number of boxes in each direction to compute the shear strain and displacement vorticity grid

	Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
	init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame
	Nframes = Nentries - init_frame									# number of frames available for the calculation

	dt = Nframes + dt if dt <= 0 else dt	# length of the interval of time for which displacements are calculated

	# NAMING

	attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'init_frame': init_frame, 'dt': dt,
		'int_max': int_max, 'Ncases': Ncases, 'r_cut': r_cut, 'sigma': sigma,
		'box_size': box_size, 'x_zero': centre[0], 'y_zero': centre[1]}			# attributes displayed in filenames
	naming_Css = naming.Css()													# Css naming object
	Css_filename, = naming_Css.filename(**attributes)							# Css filename
	naming_Ccc = naming.Ccc()													# Ccc naming object
	Ccc_filename, = naming_Ccc.filename(**attributes)							# Css filename

	# STANDARD OUTPUT

	if 'SLURM_JOB_ID' in envvar:	# script executed from Slurm job scheduler
		slurm_output(joinpath(data_dir, 'out'), naming_Css, attributes)

	# MODE SELECTION

	if get_env('COMPUTE', default=False, vartype=bool):	# COMPUTE mode

		startTime = datetime.now()

		# VARIABLE DEFINITIONS

		wrap_file_name = get_env('WRAPPED_FILE',
			default=joinpath(data_dir, naming.wrapped_trajectory_file))		# wrapped trajectory file (.gsd)
		unwrap_file_name = get_env('UNWRAPPED_FILE',
			default=joinpath(data_dir, naming.unwrapped_trajectory_file))	# unwrapped trajectory file (.dat)

		grid_points = np.array([(x, y) for x in\
			relative_positions(np.linspace(- box_size*(1 - 1./Ncases)/2,
			box_size*(1 - 1./Ncases)/2, Ncases, endpoint=True) + centre[0],
			0, parameters['box_size'])\
			for y in\
			relative_positions(np.linspace(- box_size*(1 - 1./Ncases)/2,
			box_size*(1 - 1./Ncases)/2, Ncases, endpoint=True) + centre[1],
			0, parameters['box_size'])\
			])	# grid points at which shear strain will be evaluated

		times = np.array(list(OrderedDict.fromkeys(map(
			lambda x: int(x),
			np.linspace(init_frame, Nentries - dt - 1, int_max)
			))))	# frames at which shear strain will be calculated

		display_grid = get_env('DISPLAY_GRID', default=0, vartype=int)	# index of map in list of variable maps to display

		# SHEAR STRAIN, DISPLACEMENT VORTICITY AND THEIR CORRELATIONS

		with open(wrap_file_name, 'rb') as wrap_file,\
			open(unwrap_file_name, 'rb') as unwrap_file:	# opens wrapped and unwrapped trajectory files

			w_traj = Gsd(wrap_file, prep_frames=prep_frames)	# wrapped trajectory object
			u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object
			Sgrid, Cgrid = tuple(np.transpose(list(map(lambda time:
				strain_vorticity_grid(parameters['box_size'], Ncases,
				grid_points, time, dt, w_traj, u_traj, sigma, r_cut)
				, times)), (1, 0, 2, 3)))						# lists of shear strain and displacement vorticity correlations

		Css2D, Ccc2D = tuple(map(corField2D_scalar_average, [Sgrid, Cgrid]))	# shear strain and displacement vorticity fields correlations

		# SAVING

		sgrid = Sgrid[display_grid]
		cgrid = Cgrid[display_grid]

		with open(joinpath(data_dir, Css_filename), 'wb') as Css_dump_file,\
			open(joinpath(data_dir, Ccc_filename), 'wb') as Ccc_dump_file:
			pickle.dump([sgrid, Css2D], Css_dump_file)
			pickle.dump([cgrid, Ccc2D], Ccc_dump_file)

		# EXECUTION TIME

		print("Execution time: %s" % (datetime.now() - startTime))

	if get_env('PLOT', default=False, vartype=bool):	# PLOT mode

		# DATA

		with open(joinpath(data_dir, Css_filename), 'rb') as Css_dump_file,\
			open(joinpath(data_dir, Ccc_filename), 'rb') as Ccc_dump_file:
			sgrid, Css2D = pickle.load(Css_dump_file)
			cgrid, Ccc2D = pickle.load(Ccc_dump_file)

	if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

		# PLOT

		r_max = get_env('R_MAX', default=_r_max, vartype=float)	# half size of the box showed for 2D correlation
		r_max = box_size/2 if r_max < 0 else r_max

		fig_Css, ax_Css, gridcircle_Css = plot(sgrid, Css2D, box_size,
			'\epsilon_{xy}', naming_Css)	# plotting shear strain map and correlation
		fig_Ccc, ax_Ccc, gridcircle_Ccc = plot(cgrid, Ccc2D, box_size,
			'\omega', naming_Ccc)			# plotting displacement vorticity map and correlation

		if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
			plt.show()
