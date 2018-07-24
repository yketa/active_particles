"""
Module css calculates or plots shear strain and displacement vorticity as well
as their correlations.

Files are saved according to active_particles.naming.Css (shear strain) and
active_particles.naming.Ccc (displacement vorticity) standards.

A brief description of the algorithm can be found at:
https://yketa.github.io/UBC_2018_Wiki/#Shear%20strain%20and%20vorticity%20fields

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
FROM_FT [COMPUTE or PLOT mode] : bool
	Calculates shear strain and displacement vorticity in Fourier space rather
	than in real space, or plots the resulting correlations in real space.
	DEFAULT: False
SUPERIMPOSE_C44 [FROM_FT and SHOW mode] : bool
	Superimpose C44 curve to Css(r, theta) curve.
	DEFAULT: False
FITTING_LINE [FROM_FT and SHOW mode] : bool
	Display adjustable fitting line on graph of projection of shear strain on
	cos(4 \\theta).
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
R_CUT [not(FROM_FT) mode] : float
	Cut-off radius for coarse graining function.
	NOTE: Value of R_CUT is then multiplied by average particle diameter from
	      simulation parameters file.
	DEFAULT: active_particles.analysis.css._r_cut
SIGMA [not(FROM_FT) mode] : float
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
DISPLAY_GRID [COMPUTE and not(FROM_FT) mode] : int
	Index of map in list of variable maps to display.
	DEFAULT : 0
FONT_SIZE : int
    Plot font size.
    DEFAULT: active_particles.analysis.css._font_size
POINTS_X_C44 [FROM_FT and PLOT mode] : int
	Number of radii at which to evaluate integrated strain correlation.
	DEFAULT: active_particles.analysis.css._points_x_c44
POINTS_THETA_C44 [FROM_FT and PLOT mode] : int
	Number of angles to evaluate integrated strain correlation.
	DEFAULT: active_particles.analysis.css._points_theta_c44
Y_MIN_C44 [FROM_FT and PLOT mode] : float
	Minimum plot value for C44.
	DEFAULT: active_particles.analysis.css._y_min_c44
Y_MAX_C44 [FROM_FT and PLOT mode] : float
	Maximum plot value for C44.
	DEFAULT: active_particles.analysis.css._y_max_c44
R_MIN_C44 [FROM_FT and PLOT mode] : float
	Minimum radius in average particle separation for C44 calculation.
	DEFAULT: active_particles.analysis.css._r_min_c44
R_MAX_C44 [FROM_FT and PLOT mode] : float
	Maximum radius in average particle separation for C44 calculation.
	DEFAULT: active_particles.analysis.css._r_max_c44
SMOOTH [FROM_FT and PLOT mode] : float
	C44 Gaussian smoothing length scale in units of average particle
	separation.
	DEFAULT: 0
R_CUT_FOURIER [FROM_FT and PLOT mode] : float
    Initial wave length Gaussian cut-off radius in units of average particle
	separation.
    DEFAULT: active_particles.analysis.css._r_cut_fourier
SLOPE_C44 [FROM_FT, PLOT and FITTING_LINE mode] : slope
	Initial slope for fitting line in C44 figure.
	DEFAULT: active_particles.analysis.css._slope0_c44
SLOPE_MIN_C44 [FROM_FT, PLOT and FITTING_LINE mode] : float
	Minimum slope for fitting line in C44 figure.
	DEFAULT: active_particles.analysis.css._slope_min_c44
SLOPE_MAX_C44 [FROM_FT, PLOT and FITTING_LINE mode] : float
	Maximum slope for fitting line in C44 figure.
	DEFAULT: active_particles.analysis.css._slope_max_c44
THETA [FROM_FT and PLOT mode] : float
	Angle at which to evaluate Css(r, theta) (in multiple of pi).
	DEFAULT:  active_particles.analysis.css._theta
POINTS_X_THETA [FROM_FT and PLOT mode] : int
	Number of radii at which to evaluate Css(r, theta).
	DEFAULT: POINTS_X_C44
Y_MIN_THETA [FROM_FT and PLOT mode] : float
	Minimum plot value for Css(r, theta).
	DEFAULT: active_particles.analysis.css._y_min_theta
Y_MAX_THETA [FROM_FT and PLOT mode] : float
	Maximum plot value for Css(r, theta).
	DEFAULT: active_particles.analysis.css._y_max_theta
R_MIN_THETA [FROM_FT and PLOT mode] : float
	Minimum radius in average particle separation for Css(r, theta)
	calculation.
	DEFAULT: active_particles.analysis.css._r_min_theta
R_MAX_THETA [FROM_FT and PLOT mode] : float
	Maximum radius in average particle separation for Css(r, theta)
	calculation.
	DEFAULT: active_particles.analysis.css._r_max_theta
SLOPE_THETA [FROM_FT, PLOT and FITTING_LINE mode] : slope
	Initial slope for fitting line in Css(r, theta) figure.
	DEFAULT: SLOPE_C44
SLOPE_MIN_THETA [FROM_FT, PLOT and FITTING_LINE mode] : float
	Minimum slope for fitting line in Css(r, theta) figure.
	DEFAULT: SLOPE_MIN_C44
SLOPE_MAX_THETA [FROM_FT, PLOT and FITTING_LINE mode] : float
	Maximum slope for fitting line in Css(r, theta) figure.
	DEFAULT: SLOPE_MAX_C44

Output
------
[COMPUTE and not(FROM_FT) mode]
> Prints neigbours grid computation time and execution time.
> Saves a computed map of shear strain and the averaged shear strain
correlation according to active_particles.naming.Css standards in
DATA_DIRECTORY.
> Saves a computed map of displacement vorticity and the averaged displacement
vorticity correlations according to active_particles.naming.Ccc standards
in DATA_DIRECTORY.
[COMPUTE and FROM_FT mode]
> Saves average square norm of shear strain Fourier transforms according to
active_particles.naming.Css standards in DATA_DIRECTORY.
> Saves average square norm of displacement vorticity Fourier transforms
according to active_particles.naming.Ccc standards in DATA_DIRECTORY.
[SHOW or PLOT mode]
> Plots data map and/or correlation for shear strain and/or displacement
vorticity.
[SAVE and not(FROM_FT) mode]
> Saves data map and correlation figures in DATA_DIRECTORY.
[SAVE and FROM_FT mode]
> Saves shear strain correlation figure in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output, linframes
from active_particles.dat import Dat, Gsd
from active_particles.maths import relative_positions, wave_vectors_2D,\
	FFT2Dfilter, gaussian_smooth_1D
from active_particles.quantities import nD0_active

from active_particles.analysis.neighbours import NeighboursGrid
from active_particles.analysis.correlations import corField2D_scalar_average
from active_particles.analysis.correlations import CorGrid
from active_particles.analysis.coarse_graining import GaussianCG,\
	CoarseGraining
from active_particles.analysis.cuu import displacement_grid, Cnn
from active_particles.analysis.number import count_particles

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

from math import ceil

import numpy as np

import pickle

from operator import itemgetter

from datetime import datetime

import matplotlib as mpl
if not(get_env('SHOW', default=False, vartype=bool)):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from active_particles.plot.mpl_tools import GridCircle, FittingLine

# DEFAULT VARIABLES

_r_cut = 2	# default cut-off radius for coarse graining function
_r_max = 20	# default half size of the box showed for 2D correlation

_c_min = -0.2	# default minimum value for correlations
_c_max = 1		# default maximum value for correlations

_font_size = 10	# default plot font size

_slope0_c44 = -2    # default initial slope for fitting line slider in C44 figure
_slope_min_c44 = -5 # default minimum slope for fitting line slider in C44 figure
_slope_max_c44 = 0  # default maximum slope for fitting line slider in C44 figure

_points_x_c44 = 100     # default number of radii at which to evaluate integrated strain correlation
_points_theta_c44 = 100 # default number of angles to evaluate integrated strain correlation

_y_min_c44 = 1e-4   # default minimum C44 value
_y_max_c44 = 2e-1   # default maximum C44 value
_r_min_c44 = 1      # default minimum radius over average particle separation value for C44 figure
_r_max_c44 = 20     # default maximum radius over average particle separation value for C44 figure

_theta = 0	# default angle at which to evaluate Css(r, theta) (in multiple of pi)

_y_min_theta = 1e-4	# default minimum Css(r, theta) value
_y_max_theta = 2e-1	# default maximum Css(r, theta) value
_r_min_theta = 1	# default minimum radius over average particle separation value for Css(r, theta) figure
_r_max_theta = 20	# default maximum radius over average particle separation value for Css(r, theta) figure

_r_cut_fourier = 0	# default initial wave length Gaussian cut-off radius

# FUNCTIONS AND CLASSES

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
		Frame at which shear strain and displacement vorticity will be
		calculated.
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

	print("Neighbours grid computation time (time = %e): %s" %
		(time, datetime.now() - startTime0))	# neighbours grid computation time

	# SHEAR STRAIN AND DISPLACEMENT VORTICITY GRIDS CALCULATION

	displacements = u_traj.displacement(time, time + dt)	# array of particle displacements

	sgrid, cgrid = tuple(np.transpose(list(map(
		lambda point: strain_vorticity(point, time, dt, positions,
		displacements, sigma, r_cut, box_size, neighbours_grid),
		grid_points))))	# shear strain and displacement vorticity lists

	correct_grid = lambda grid: np.transpose(
		np.reshape(grid, (Ncases, Ncases)))[::-1]	# get grids with the same orientation as positions
	return correct_grid(sgrid), correct_grid(cgrid)	# shear strain and displacement vorticity grids

def strain_vorticity_fftsqnorm_grid(box_size, new_box_size, centre, Ncases,
	time, dt, w_traj, u_traj):
	"""
	Calculates grids of square norm of fast Fourier transforms of (linearised)
	shear strain and displacement vorticity from fast Fourier transform of
	displacement grid.

	Parameters
	----------
	box_size : float
		Length of the system's square box.
	new_box_size : float
		Length of the considered system's square box.
	centre : float array
		Centre of the box.
	Ncases : int
		Number of boxes in each direction to compute the displacements.
	time : int
		Frame at which displacements will be calculated.
	dt : int
		Length of the interval of time for which the displacements are
		calculated.
	w_traj : active_particles.dat.Gsd
		Wrapped trajectory object.
	u_traj : active_particles.dat.Dat
		Unwrapped trajectory object.

	Returns
	-------
	FFTsgridsqnorm : 2D Numpy array
		Square norm of shear strain fast Fourier transform grid.
	FFTcgridsqnorm : 2D Numpy array
		Square norm of displacement vorticity fast Fourier transform grid.
	"""

	# DISPLACEMENT GRID FOURIER TRANSFORM

	ugrid = displacement_grid(box_size, new_box_size, centre, Ncases, time, dt,
		w_traj, u_traj)							# displacment grid
	FFTugrid = np.fft.fft2(ugrid, axes=(0, 1))	# displacement grid Fourier transform

	# SHEAR STRAIN AND DISPLACEMENT VORTICITY FOURIER TRANSFORM CALCULATION

	wave_vectors = wave_vectors_2D(Ncases, Ncases, d=new_box_size/Ncases)	# wave vectors corresponding to displacement grid Fourier transform
	eikxm1 = np.exp(1j*wave_vectors[:, :, 0]) - 1							# grid exponential of i times wave vectors' x-coordinates - 1
	eikym1 = np.exp(1j*wave_vectors[:, :, 1]) - 1							# grid exponential of i times wave vectors' y-coordinates - 1

	FFTsgrid = (eikxm1*FFTugrid[:, :, 1] + eikym1*FFTugrid[:, :, 0])/2	# shear strain Fourier transform grid
	FFTcgrid = eikxm1*FFTugrid[:, :, 1] - eikym1*FFTugrid[:, :, 0]		# displacement vorticity Fourier transform grid

	return np.conj(FFTsgrid)*FFTsgrid, np.conj(FFTcgrid)*FFTcgrid

def suptitle():
	"""
	Returns suptitles for figures.
	"""

	return (
		r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e,$'
		% (parameters['N'], parameters['density'], parameters['vzero'],
		parameters['dr'])
		+ r'$\bar{N}=%.2e, \Delta t=%.2e, nD_0 \Delta t=%.2e, L=%.2e$'
		% (Nmean, dt*parameters['period_dump']*parameters['time_step'],
		nD0*dt*parameters['period_dump']*parameters['time_step'], box_size)
		+ '\n' + r'$x_0=%.2e, y_0=%.2e,$' % centre
		+ r'$S_{init}=%.2e$' % init_frame
		+ r'$, S_{max}=%.2e, N_{cases}=%.2e$' % (int_max, Ncases)
		+ r'$, dL = %.2e a$' % (np.sqrt(Nmean)/Ncases))

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
	gridcircle : active_particles.plot.mpl_tools.GridCircle
		Grid circle object.
	"""

	# GRID AND CORRELATION FIGURE

	cmap = plt.cm.jet

	fig, ax = plt.subplots(1, 2)	# 1 x 2 figure

	fig.set_size_inches(16, 16)		# figure size
	fig.subplots_adjust(wspace=0.4)	# width space
	fig.subplots_adjust(hspace=0.3)	# height space

	fig.suptitle(suptitle()
		+ r'$, r_{cut}=%.2e, \sigma=%.2e$' % (r_cut, sigma))

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

	Cmin = np.max((np.min(corr), _c_min))
	Cmax = np.max((np.min(corr), _c_max))

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

	fig_gc.suptitle(suptitle()
		+ r'$, r_{cut}=%.2e, \sigma=%.2e$' % (r_cut, sigma))

	ax_grid.set_xlabel(r'$x$')
	ax_grid.set_ylabel(r'$y$')
	ax_grid.set_title('2D ' + r'$%s$' % C)
	cb_gc.set_label(r'$%s$' % C, labelpad=20, rotation=270)

	ax_plot.set_xlabel(r'$\theta$')
	ax_plot.set_ylabel(r'$%s(r, \theta)$' % C)

	# RETURN

	return fig, ax, gridcircle

class StrainCorrelations:
	"""
	Manipulate and plot strain correlations from its Fourier transform.
	"""

	def __init__(self, wave_vectors, FFTCss):
		"""
		Parameters
		----------
		wave_vectors : (_, _, 2) array-like
			Wave vectors at which Fourier transform was calculated.
		FFTCss : 2D Numpy array
			Strain correlations Fourrier transform.
		"""

		self.wave_vectors = wave_vectors
		self.strain_correlations_FFT = FFTCss

	def strain_correlations(self, r_cut=0):
		"""
		Computes strain correlations from inverse fast Fourier transform of
		self.strain_correlations_FFT with low wave lengths Gaussian cut at
		r_cut.

		Parameters
		----------
		r_cut : float
			Wave length Gaussian cut-off radius, equivalent to coarse-graining
			cut-off radius. (default: 0)

		Returns
		-------
		sc : Numpy array
			Strain correlations.
		"""

		sc = (FFT2Dfilter(self.strain_correlations_FFT,
			wave_vectors=self.wave_vectors)
			.gaussian_filter(np.sqrt(2)*r_cut)
			.get_signal()).real
		return sc/sc[0, 0]	# correlation normalisation

	def plot(self, box_size, r_max_css, av_p_sep,
		points_x_c44=_points_x_c44, points_theta_c44=_points_theta_c44,
		y_min_c44=_y_min_c44, y_max_c44=_y_max_c44,
		r_min_c44=_r_min_c44, r_max_c44=_r_max_c44,
		smooth=0, r_cut=0,
		theta=0, points_x_theta=None,
		y_min_theta=_y_min_theta, y_max_theta=_y_max_theta,
		r_min_theta=_r_min_theta, r_max_theta=_r_max_theta,
		superimpose_c44=False):
		"""
		Plots strain correlations with slider for cut-off radius.

		Parameters
		----------
		box_size : float
			System length.
		r_max_css : float
			Maximum radius in infinite norm for strain correlations plots.
		av_p_sep : float
			Average particle separation.
		points_x_c44 : int
			Number of radii at which to evaluate integrated strain correlation.
			(default: active_particles.plot.c44._points_x)
		points_theta_c44 : int
			Number of angles to evaluate integrated strain correlation.
			(default: active_particles.plot.c44._points_theta)
		y_min_c44 : float
			Minimum plot value for C44.
			(default: active_particles.analysis.css._y_min_c44)
		y_max_c44 : float
			Maximum plot value for C44.
			(default: active_particles.analysis.css._y_max_c44)
		r_min_c44 : float
			Minimum radius in average particle separation for C44 calculation.
			(default: active_particles.analysis.css._r_min_c44)
		r_max_c44 : float
			Maximum radius in average particle separation for C44 calculation.
			(default: active_particles.analysis.css._r_max_c44)
		smooth : float
			C44 Gaussian smoothing length scale in units of average particle
			separation. (default: 0)
		r_cut : float
			Initial wave length Gaussian cut-off radius in units of average
			particles separation. (default: 0)
		theta : float
			Angle at which to plot Css(r, theta). (default: 0)
		points_x_theta : int
			Number of radii at which to evaluate Css(r, theta). (default: None)
			NOTE: if points_x_theta == None, then it is taken equal to
			      points_x_c44.
		y_min_theta : float
			Minimum plot value for Css(r, theta).
			(default: active_particles.analysis.css._y_min_theta)
		y_max_theta : float
			Maximum plot value for Css(r, theta).
			(default: active_particles.analysis.css._y_max_theta)
		r_min_theta : float
			Minimum radius in average particle separation for Css(r, theta)
			calculation.
			(default: active_particles.analysis.css._r_min_theta)
		r_max_theta : float
			Maximum radius in average particle separation for Css(r, theta)
			calculation.
			(default: active_particles.analysis.css._r_max_theta)
		superimpose_c44 : bool
			Superimpose C44 curve to Css(r, theta) curve. (default: False)
		"""

		self.box_size = box_size
		self.av_p_sep = av_p_sep

		self.r_max_css = r_max_css

		self.r_cut = r_cut	# Gaussian cut-off radius in units of average particle separation

		self.cor_name = 'C_{\\varepsilon_{xy}\\varepsilon_{xy}}'	# name of plotted correlation

		grid = self.strain_correlations_corgrid()	# correlation grid

		# CSS FIGURE

		self.fig_css, self.ax_css = plt.subplots()
		self.ax_css.set_title(
			'2D ' + r'$%s(r_{cut}/a = %.2e)$' % (self.cor_name, self.r_cut))

		Cmin = np.max((np.min(grid.grid), _c_min))
		Cmax = np.max((np.min(grid.grid), _c_max))
		cmap = plt.cm.jet
		CvNorm = colors.Normalize(vmin=Cmin, vmax=Cmax)

		self.grid_plot = self.ax_css.imshow(grid.display_grid.grid,
			cmap=cmap, norm=CvNorm, extent=
			[-self.r_max_css, self.r_max_css, -self.r_max_css, self.r_max_css])

		self.colormap = plt.colorbar(self.grid_plot)	# color bar

		self.ax_slider = make_axes_locatable(self.ax_css).append_axes('bottom',
			size='5%', pad=0.6)	# slider Axes

		self.slider = Slider(self.ax_slider, r'$r_{cut}/a$', 0,
			self.r_max_css/self.av_p_sep, valinit=self.r_cut)	# slider
		self.slider.on_changed(self.update_r_cut)				# call update_r_cut when slider value is changed

		# GRID CIRCLE FIGURE

		self.grid_circle = GridCircle(grid.display_grid.grid, extent=
			[-self.r_max_css/self.av_p_sep, self.r_max_css/self.av_p_sep,
				-self.r_max_css/self.av_p_sep, self.r_max_css/self.av_p_sep],
			min=Cmin, max=Cmax)
		self.grid_circle.ax_grid.set_title(
			'2D ' + r'$%s(r_{cut}/a = %.2e)$' % (self.cor_name, self.r_cut))

		# C44 FIGURE

		self.points_x_c44 = points_x_c44
		self.points_theta_c44 = points_theta_c44
		self.y_min_c44 = y_min_c44
		self.y_max_c44 = y_max_c44
		self.r_min_c44 = r_min_c44
		self.r_max_c44 = r_max_c44
		self.smooth = smooth	# C44 Gaussian smoothing length scale in units of average particle separation

		self.toC44 = Css2DtoC44(self.box_size,
			self.points_x_c44, self.points_theta_c44,
			self.av_p_sep*self.r_min_c44, self.av_p_sep*self.r_max_c44)

		self.fig_c44, self.ax_c44 = plt.subplots()
		self.ax_c44.set_title(r'$r_{cut}/a = %.2e$' % self.r_cut)
		self.ax_c44.set_xlim([self.r_min_c44, self.r_max_c44])
		self.ax_c44.set_ylim([self.y_min_c44, self.y_max_c44])

		self.c44 = self.css2Dtoc44(grid.grid)				# list of [r, C44(r)]
		self.line_c44, = self.ax_c44.plot(
			self.c44[:, 0]/self.av_p_sep, self.c44[:, 1])	# C44 plotting line

		# CSS(r, theta) FIGURE

		self.theta = theta
		self.points_x_theta = (points_x_theta if points_x_theta != None
			else self.points_x_c44)
		self.y_min_theta = y_min_theta
		self.y_max_theta = y_max_theta
		self.r_min_theta = r_min_theta
		self.r_max_theta = r_max_theta
		self.superimpose_c44 = superimpose_c44

		self.toCsstheta = Css2DtoCsstheta(self.box_size, self.points_x_theta,
			self.av_p_sep*self.r_min_theta, self.av_p_sep*self.r_max_theta)

		self.fig_theta, self.ax_theta = plt.subplots()
		self.ax_theta.set_title(r'$r_{cut}/a = %.2e$' % self.r_cut)
		self.ax_theta.set_xlim([self.r_min_theta, self.r_max_theta])
		self.ax_theta.set_ylim([self.y_min_theta, self.y_max_theta])

		self.csstheta = self.css2Dtocsstheta(grid.grid)				# list of [r, Css(r, theta)]
		self.line_csstheta, = self.ax_theta.plot(
			self.csstheta[:, 0]/self.av_p_sep, self.csstheta[:, 1])	# Css(r, theta) plotting line

		if self.superimpose_c44:								# SUPERIMPOSE_C44 mode
			self.line_c44_superimposed, = self.ax_theta.plot(
				self.c44[:, 0]/self.av_p_sep, self.c44[:, 1])	# C44 plotting line superimposed on Css(r, theta)

	# METHODS CALLED BY SELF.PLOT()

	def strain_correlations_corgrid(self):
		"""
		Define strain correlations grid, with wave length cut-off radius
		self.r_cut equivalent to coarse-graining cut-off radius.

		Returns
		-------
		corgrid : active_particles.analysis.correlations.CorGrid
			Strain correlations CorGrid object.
		"""

		sc = self.strain_correlations(r_cut=self.r_cut*self.av_p_sep)	# strain correlations grid
		return CorGrid(sc, self.box_size, display_size=2*self.r_max_css)

	def css2Dtoc44(self, css2D):
		"""
		Returns strain correlations projected on cos(4 \\theta) from strain
		correlations grid, smoothed over self.smooth*self.av_p_sep.

		Parameters
		----------
		css2D : array-like
			Strain correlations grid.

		Returns
		-------
		c44 : Numpy array
			Array of position and projection of strain correlations on
			cos(4 \\theta) at this position.
		"""

		return self.toC44.get_C44(css2D, smooth=self.av_p_sep*self.smooth)	# list of [r, C44(r)]

	def css2Dtocsstheta(self, css2D):
		"""
		Returns strain correlations along half-line of origin (0, 0) and
		inclination self.theta.

		Parameters
		----------
		css2D : array-like
			Strain correlations grid.

		Returns
		-------
		csstheta : Numpy array
			Array of radii r and values of strain correlations (r, theta).
		"""

		return self.toCsstheta.get_Csstheta(css2D, self.theta)

	def update_r_cut(self, val):
		"""
		Updates cut-off radius on slider change.
		"""

		self.r_cut = self.slider.val	# new cut-off radius
		self.draw()						# updates figure

	def draw(self):
		"""
		Updates line position and strain correlations grid according to cut-off
		radius self.r_cut.
		"""

		self.ax_css.set_title(
			'2D ' + r'$%s(r_{cut}/a = %.2e)$' % (self.cor_name, self.r_cut))	# update title
		self.grid_circle.ax_grid.set_title(
			'2D ' + r'$%s(r_{cut} = %.2e)$' % (self.cor_name, self.r_cut))		# update title
		self.ax_c44.set_title(r'$r_{cut}/a = %.2e$' % self.r_cut)				# update title
		self.ax_theta.set_title(r'$r_{cut}/a = %.2e$' % self.r_cut)				# update title

		grid = self.strain_correlations_corgrid()		# new grid
		self.grid_plot.set_data(grid.display_grid.grid)	# plot grid
		self.grid_plot.figure.canvas.draw()				# update plot

		self.grid_circle.update_grid_plot(grid.display_grid.grid)	# plot grid

		self.c44 = self.css2Dtoc44(grid.grid)		# update C44 values
		self.line_c44.set_ydata(self.c44[:, 1])		# plot C44
		self.line_c44.figure.canvas.draw()			# update plot

		self.csstheta = self.css2Dtocsstheta(grid.grid)		# update Css(r, theta) values
		self.line_csstheta.set_ydata(self.csstheta[:, 1])	# plot Css(r, theta)
		self.line_csstheta.figure.canvas.draw()				# update plot
		if self.superimpose_c44: self.line_c44_superimposed.set_ydata(
			self.c44[:, 1])									# update plot

class Css2DtoC44:
    """
    Calculates C44 as projection of 2D shear strain correlations on
    cos(4\\theta).
    """

    def __init__(self, box_size, points_x, points_theta,
        r_min, r_max):
        """
        Sets parameters for C44 integration.

        Parameters
        ----------
        box_size : float
            Size of the square box to consider.
        points_x : int
            Number of radii at which to evaluate integrated strain correlation.
        points_theta : int
            Number of angles to evaluate integrated strain correlation.
        r_min : float
            Minimum radius.
        r_max : float
            Maximum radius.
        """

        self.box_size = box_size

        self.points_x = points_x
        self.points_theta = points_theta

        self.r_min = r_min
        self.r_max = r_max
        self.c44_x = np.linspace(self.r_min, self.r_max, self.points_x)

    def get_C44(self, Css2D, smooth=0):
        """
        From 2D strain correlations Css2D, returns values of C44 at
        self.points_x radii between r_max and r_min.

        Parameters
        ----------
        Css2D : 2D array-like
            Shear strain correlation grid.
		smooth : float
			C44 Gaussian smoothing length scale. (default: 0)

        Returns
        -------
        C44 : 2D Numpy array
            List of [r, C44(r)].
        """

        self.css2Dgrid = CorGrid(Css2D, self.box_size)	# shear strain 2D CorGrid object

        self.c44 = np.array(list(map(
            lambda r: self.css2Dgrid.integrate_over_angles(r,
            	projection=lambda theta: np.cos(4*theta)/np.pi,
            	points_theta=self.points_theta),
            self.c44_x)))

        return np.array(list(map(lambda x, y: [x, y],
			*(self.c44_x, gaussian_smooth_1D(self.c44_x, self.c44, smooth)))))

class Css2DtoCsstheta:
    """
    Calculates Css(r, theta) as values of strain correlations along the
	half-line of origin (0, 0) and inclination of angle theta.
    """

    def __init__(self, box_size, points_x, r_min, r_max):
        """
        Sets parameters for Css(r, theta) calculation.

        Parameters
        ----------
        box_size : float
            Size of the square box to consider.
        points_x : int
            Number of radii at which to evaluate strain correlations.
        r_min : float
            Minimum radius.
        r_max : float
            Maximum radius.
        """

        self.box_size = box_size

        self.points_x = points_x

        self.r_min = r_min
        self.r_max = r_max
        self.csstheta_x = np.linspace(self.r_min, self.r_max, self.points_x)

    def get_Csstheta(self, Css2D, theta=0):
        """
        From 2D strain correlations Css2D, returns values of Css(r, theta) at
        self.points_x radii between r_max and r_min.

        Parameters
        ----------
        Css2D : 2D array-like
            Shear strain correlation grid.
		theta : float
			Inclination of half-line in Css2D grid. (default: 0)

        Returns
        -------
        Csstheta : 2D Numpy array
            List of [r, Css(r, theta)].
        """

        self.css2Dgrid = CorGrid(Css2D, self.box_size)	# shear strain 2D CorGrid object
        self.csstheta = np.array(list(map(
			lambda r: [r, self.css2Dgrid.get_value_polar(r, theta)],
            self.csstheta_x)))

        return self.csstheta

def plot_fft():
	"""
	Plots shear strain correlations.

	Returns
	-------
	sc : active_particles.analysis.css.StrainCorrelations
		Strain correlations object.
	"""

	cor_name = r'$C_{\varepsilon_{xy}\varepsilon_{xy}}$'	# name of plotted correlation

	wave_vectors = wave_vectors_2D(*FFTsgridsqnorm.shape, d=box_size/Ncases)	# wave vectors at which Fourier transform was calculated

	sc = StrainCorrelations(wave_vectors, FFTsgridsqnorm)

	av_p_sep = box_size/np.sqrt(Nmean)
	sc.plot(parameters['box_size'], r_max, av_p_sep,
		points_x_c44=points_x_c44, points_theta_c44=points_theta_c44,
		y_min_c44=y_min_c44, y_max_c44=y_max_c44,
		r_min_c44=r_min_c44, r_max_c44=r_max_c44,
		smooth=smooth, r_cut=r_cut_fourier,
		theta=theta, points_x_theta=points_x_theta,
		y_min_theta=y_min_theta, y_max_theta=y_max_theta,
		r_min_theta=r_min_theta, r_max_theta=r_max_theta,
		superimpose_c44=superimpose_c44)

	# CSS FIGURE

	sc.fig_css.set_size_inches(16, 16)

	sc.fig_css.suptitle(suptitle())

	sc.ax_css.set_xlabel(r'$x$')
	sc.ax_css.set_ylabel(r'$y$')
	sc.colormap.set_label(cor_name, labelpad=20, rotation=270)

	if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
		sc.ax_slider.set_visible(False)	# hide slider
		image_name, = naming_Css.image().filename(
			**{**attributes, 'r_cut': r_cut_fourier})
		sc.fig_css.savefig(joinpath(data_dir, image_name))
		sc.ax_slider.set_visible(True)	# show slider

	# GRID CIRCLE FIGURE

	fig_gc, (ax_grid, ax_plot), cb_gc = sc.grid_circle.get_fig_ax_cmap()

	sc.grid_circle.fig.suptitle(suptitle())

	sc.grid_circle.fig.set_size_inches(16, 16)		# figure size
	sc.grid_circle.fig.subplots_adjust(wspace=0.4)	# width space
	sc.grid_circle.fig.subplots_adjust(hspace=0.3)	# height space

	sc.grid_circle.ax_grid.set_xlabel(r'$x/a$')
	sc.grid_circle.ax_grid.set_ylabel(r'$y/a$')
	sc.grid_circle.colormap.set_label(cor_name, labelpad=20, rotation=270)

	sc.grid_circle.ax_plot.set_xlabel(r'$\theta$')
	sc.grid_circle.ax_plot.set_ylabel(cor_name + r'$(r, \theta)$')

	# C44 FIGURE

	sc.fig_c44.suptitle(suptitle())

	sc.fig_c44.set_size_inches(16, 16)	# figure size

	sc.ax_c44.set_xlabel(r'$r/a$' + ' ' + r'$(a = L/\sqrt{\bar{N}})$')
	sc.ax_c44.set_ylabel(r'$C_4^4(r) = \frac{1}{\pi}\int_0^{2\pi}d\theta$'
		+ ' ' + cor_name + r'$(r, \theta)$'
		+ ' ' + r'$\cos4\theta$'
		+ ' ' + (r'($\sigma_{smooth}/a=%.2e$)' % smooth if smooth!=0 else ''))
	sc.ax_c44.set_xlim([r_min_c44, r_max_c44])
	sc.ax_c44.set_ylim([y_min_c44, y_max_c44])
	sc.ax_c44.set_yscale('log')
	sc.ax_c44.set_xscale('log')

	if get_env('FITTING_LINE', default=False, vartype=bool):	# FITTING_LINE mode
		sc.fl_c44 = FittingLine(sc.ax_c44, slope0_c44, slope_min_c44,
			slope_max_c44, x_fit='(r/a)', y_fit='C_4^4(r)')		# add fitting line to plot

	# CSS(r, theta) FIGURE

	sc.fig_theta.suptitle(suptitle())

	sc.fig_theta.set_size_inches(16, 16)	# figure size

	sc.ax_theta.set_xlabel(r'$r/a$' + ' ' + r'$(a = L/\sqrt{\bar{N}})$')
	sc.ax_theta.set_ylabel('%s' % cor_name + r'$(r, \theta=%.2e)$' % theta)
	sc.ax_theta.set_xlim([r_min_theta, r_max_theta])
	sc.ax_theta.set_ylim([y_min_theta, y_max_theta])
	sc.ax_theta.set_yscale('log')
	sc.ax_theta.set_xscale('log')

	if sc.superimpose_c44:	# SUPERIMPOSE_C44 mode
		sc.ax_theta.set_ylabel(r'$C(t)$')
		plt.sca(sc.ax_theta)
		sc.ax_theta.add_artist(plt.legend(handles=[
			Line2D([0], [0], color=sc.line_csstheta.get_color(),
				label='%s' % cor_name + r'$(r, \theta=%.2e)$' % theta),
			Line2D([0], [0], color=sc.line_c44_superimposed.get_color(),
				label=r'$C_4^4(r)$')]))
	else:
		sc.ax_theta.set_ylabel('%s' % cor_name + r'$(r, \theta=%.2e)$' % theta)

	if get_env('FITTING_LINE', default=False, vartype=bool):					# FITTING_LINE mode
		sc.fl_theta = FittingLine(sc.ax_theta, slope0_theta, slope_min_theta,
			slope_max_theta, x_fit='(r/a)', y_fit='%s(r, \\theta)' % cor_name)	# add fitting line to plot

	# RETURN

	return sc

# SCRIPT

if __name__ == '__main__':	# executing as script

	# VARIABLE DEFINITIONS

	data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

	wrap_file_name = get_env('WRAPPED_FILE',
		default=joinpath(data_dir, naming.wrapped_trajectory_file))	# wrapped trajectory file (.gsd)

	dt = get_env('TIME', default=-1, vartype=int)	# lag time for displacement

	init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
	int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of intervals of length dt considered in correlations calculations

	parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
	with open(parameters_file, 'rb') as param_file:
		parameters = pickle.load(param_file)				# parameters hash table

	r_cut = parameters['a']*get_env('R_CUT', default=_r_cut, vartype=float)	# cut-off radius for coarse graining function
	sigma = get_env('SIGMA', default=r_cut, vartype=float)					# length scale of the spatial extent of the coarse graining function

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

	times = linframes(init_frame, Nentries - dt, int_max)	# frames at which shear strain is calculated

	from_ft = get_env('FROM_FT', default=False, vartype=bool)	# calculation of shear strain and vorticity in real space (False) or in Fourier space (True)

	# NAMING

	attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'init_frame': init_frame, 'dt': dt,
		'int_max': int_max, 'Ncases': Ncases, 'r_cut': r_cut, 'sigma': sigma,
		'box_size': box_size, 'x_zero': centre[0], 'y_zero': centre[1]}			# attributes displayed in filenames
	naming_Css = naming.Css(from_ft=from_ft)									# Css naming object
	Css_filename, = naming_Css.filename(**attributes)							# Css filename
	naming_Ccc = naming.Ccc(from_ft=from_ft)									# Ccc naming object
	Ccc_filename, = naming_Ccc.filename(**attributes)							# Css filename

	# STANDARD OUTPUT

	if 'SLURM_JOB_ID' in envvar:	# script executed from Slurm job scheduler
		slurm_output(joinpath(data_dir, 'out'), naming_Css, attributes)

	# MODE SELECTION

	if get_env('COMPUTE', default=False, vartype=bool):	# COMPUTE mode

		startTime = datetime.now()

		# VARIABLE DEFINITIONS

		unwrap_file_name = get_env('UNWRAPPED_FILE',
			default=joinpath(data_dir, naming.unwrapped_trajectory_file))	# unwrapped trajectory file (.dat)

		if not(from_ft):	# calculation of shear strain and vorticity in real space

			grid_points = np.array([(x, y) for x in\
				relative_positions(np.linspace(- box_size*(1 - 1./Ncases)/2,
				box_size*(1 - 1./Ncases)/2, Ncases, endpoint=True) + centre[0],
				0, parameters['box_size'])\
				for y in\
				relative_positions(np.linspace(- box_size*(1 - 1./Ncases)/2,
				box_size*(1 - 1./Ncases)/2, Ncases, endpoint=True) + centre[1],
				0, parameters['box_size'])\
				])	# grid points at which shear strain will be evaluated

			display_grid = get_env('DISPLAY_GRID', default=0, vartype=int)	# index of map in list of variable maps to display

			# SHEAR STRAIN, DISPLACEMENT VORTICITY AND THEIR CORRELATIONS

			with open(wrap_file_name, 'rb') as wrap_file,\
				open(unwrap_file_name, 'rb') as unwrap_file:	# opens wrapped and unwrapped trajectory files

				w_traj = Gsd(wrap_file, prep_frames=prep_frames)	# wrapped trajectory object
				u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object
				Sgrid, Cgrid = tuple(np.transpose(list(map(lambda time:
					strain_vorticity_grid(parameters['box_size'], Ncases,
					grid_points, time, dt, w_traj, u_traj, sigma, r_cut),
					times)), (1, 0, 2, 3)))							# lists of shear strain and displacement vorticity correlations

			Css2D, Ccc2D = tuple(map(corField2D_scalar_average,
				[Sgrid, Cgrid]))	# shear strain and displacement vorticity fields correlations

			# SAVING

			sgrid = Sgrid[display_grid]
			cgrid = Cgrid[display_grid]

			with open(joinpath(data_dir, Css_filename),
				'wb') as Css_dump_file,\
				open(joinpath(data_dir, Ccc_filename),
				'wb') as Ccc_dump_file:
				pickle.dump([sgrid, Css2D], Css_dump_file)
				pickle.dump([cgrid, Ccc2D], Ccc_dump_file)

		else:	# calculation of shear strain and vorticity in Fourier space

			# SHEAR STRAIN AND DISPLACEMENT VORTICITY FAST FOURIER TRANSFORMS

			with open(wrap_file_name, 'rb') as wrap_file,\
				open(unwrap_file_name, 'rb') as unwrap_file:	# opens wrapped and unwrapped trajectory files

				w_traj = Gsd(wrap_file, prep_frames=prep_frames)	# wrapped trajectory object
				u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object
				FFTsgridsqnorm, FFTcgridsqnorm = tuple(np.mean(np.transpose(
					list(map(lambda time:
					strain_vorticity_fftsqnorm_grid(parameters['box_size'],
					box_size, centre, Ncases, time, dt, w_traj, u_traj),
					times)), (1, 0, 2, 3)), axis=1))				# average square norm of shear strain and displacement vorticity Fourier transforms

			# SAVING

			with open(joinpath(data_dir, Css_filename),
				'wb') as Css_dump_file,\
				open(joinpath(data_dir, Ccc_filename),
				'wb') as Ccc_dump_file:
				pickle.dump(FFTsgridsqnorm, Css_dump_file)
				pickle.dump(FFTcgridsqnorm, Ccc_dump_file)

		# EXECUTION TIME

		print("Execution time: %s" % (datetime.now() - startTime))

	if get_env('PLOT', default=False, vartype=bool):	# PLOT mode

		# DATA

		if not(from_ft):	# shear strain and vorticity in real space
			with open(joinpath(data_dir, Css_filename),
				'rb') as Css_dump_file,\
				open(joinpath(data_dir, Ccc_filename),
				'rb') as Ccc_dump_file:
				sgrid, Css2D = pickle.load(Css_dump_file)
				cgrid, Ccc2D = pickle.load(Ccc_dump_file)
		else:				# shear strain in Fourier space
			with open(joinpath(data_dir, Css_filename), 'rb') as Css_dump_file:
				FFTsgridsqnorm = pickle.load(Css_dump_file)

	if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

		# PARAMETERS

		with open(wrap_file_name, 'rb') as wrap_file:
			w_traj = Gsd(wrap_file, prep_frames=prep_frames)		# wrapped trajectory object
			Nmean = np.mean(count_particles(
				w_traj, *times, box_size=box_size, centre=centre))	# average number of particles in box
		nD0 = nD0_active(
			Nmean, parameters['vzero'], parameters['dr'], box_size)	# product of particle density and active diffusion constant
		av_p_sep = parameters['box_size']/np.sqrt(Nmean)			# average particle separation

		# PLOT

		font_size = get_env('FONT_SIZE', default=_font_size, vartype=int)       # plot font size
		mpl.rcParams.update({'font.size': font_size})

		r_max = get_env('R_MAX', default=_r_max, vartype=float)	# half size of the box showed for 2D correlation
		r_max = box_size/2 if r_max < 0 else r_max

		if not(from_ft):	# shear strain and vorticity in real space

			fig_Css, ax_Css, gridcircle_Css = plot(sgrid, Css2D, box_size,
				'\epsilon_{xy}', naming_Css)	# plotting shear strain map and correlation
			fig_Ccc, ax_Ccc, gridcircle_Ccc = plot(cgrid, Ccc2D, box_size,
				'\omega', naming_Ccc)			# plotting displacement vorticity map and correlation

		else:				# shear strain in Fourier space

			points_x_c44 = get_env('POINTS_X_C44', default=_points_x_c44,
				vartype=int)							# number of radii at which to evaluate integrated strain correlation
			points_theta_c44 = get_env('POINTS_THETA_C44',
				default=_points_theta_c44, vartype=int)	# number of angles to evaluate integrated strain correlation

			y_min_c44 = get_env('Y_MIN_C44', default=_y_min_c44, vartype=float)	# minimum plot value for C44
			y_max_c44 = get_env('Y_MAX_C44', default=_y_max_c44, vartype=float)	# maximum plot value for C44

			r_min_c44 = get_env('R_MIN_C44', default=_r_min_c44, vartype=float)	# minimum radius in average particle separation for C44 calculation
			r_max_c44 = get_env('R_MAX_C44', default=_r_max_c44, vartype=float)	# maximum radius in average particle separation for C44 calculation

			smooth = get_env('SMOOTH', default=0, vartype=float)	# C44 Gaussian smoothing length scale in units of average particle separation

			r_cut_fourier = get_env('R_CUT_FOURIER', default=_r_cut_fourier,
				vartype=float)	# initial wave length Gaussian cut-off radius in units of average particle separation

			theta = np.pi*get_env('THETA', default=_theta, vartype=float)	# angle at which to evaluate Css(r, theta)

			points_x_theta = get_env('POINTS_X_THETA', default=points_x_c44,
				vartype=int)	# number of radii at which to evaluate Css(r, theta)

			y_min_theta = get_env('Y_MIN_THETA', default=_y_min_theta,
				vartype=float)	# minimum plot value for Css(r, theta)
			y_max_theta = get_env('Y_MAX_THETA', default=_y_max_theta,
				vartype=float)	# maximum plot value for Css(r, theta)

			r_min_theta = get_env('R_MIN_THETA', default=_r_min_theta,
				vartype=float)	# minimum radius in average particle separation for Css(r, theta) calculation
			r_max_theta = get_env('R_MAX_THETA', default=_r_max_theta,
				vartype=float)	# maximum radius in average particle separation for Css(r, theta) calculation

			superimpose_c44 = get_env('SUPERIMPOSE_C44',
				default=False, vartype=bool)	# superimpose C44 curve to Css(r, theta) curve

			if get_env('FITTING_LINE', default=False, vartype=bool):	# FITTING_LINE mode

				slope0_c44 = get_env('SLOPE_C44',
					default=_slope0_c44, vartype=float)		# initial slope for fitting line in C44 figure
				slope_min_c44 = get_env('SLOPE_MIN_C44',
					default=_slope_min_c44, vartype=float)	# minimum slope for fitting line in C44 figure
				slope_max_c44 = get_env('SLOPE_MAX_C44',
					default=_slope_max_c44, vartype=float)	# maximum slope for fitting line in C44 figure

				slope0_theta = get_env('SLOPE_THETA',
					default=slope0_c44, vartype=float)		# initial slope for fitting line in Css(r, theta) figure
				slope_min_theta = get_env('SLOPE_MIN_THETA',
					default=slope_min_c44, vartype=float)	# minimum slope for fitting line in Css(r, theta) figure
				slope_max_theta = get_env('SLOPE_MAX_THETA',
					default=slope_max_c44, vartype=float)	# maximum slope for fitting line in Css(r, theta) figure

			plot = plot_fft()	# plotting shear strain correlations with slider for cut-off radius

		if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
			plt.show()
