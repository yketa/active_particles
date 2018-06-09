"""
Module ctt calculates or plots mean squared cross and dot products of
normalised wave vectors and displacement Fourier transform, as well as strain
correlations from these variables.

Files are saved according to active_particles.naming.Ctt (cross product) and
active_particles.naming.Cll (dot product) naming standard.

Environment modes
-----------------
COMPUTE : bool
	Compute displacement Fourier transform and mean squared cross and dot
	products with normalised wave vectors.
	DEFAULT: False
PLOT : bool
	Plot saved cylindrical averages of mean squared cross and dot products of
	normalised wave vectors and displacement Fourier transform as functions of
	wave vector norm.
	DEFAULT: False
SHOW [COMPUTE or PLOT mode] : bool
	Show graphs.
	DEFAULT: False
SAVE [COMPUTE or PLOT mode] : bool
	Save graphs.
	DEFAULT: False
FITTING_LINE [SHOW mode] : bool
	Display adjustable fitting line on graphs.
	DEFAULT: False
STRAIN_CORRELATIONS [SHOW mode] : bool
	Computes strain correlations from cross and dot products of normalised wave
	vectors and displacement Fourier transform as functions of wave vector
	norm and plots it, with adjustable cut-off radius.
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
	Maximum number of intervals of length dt considered for the calculation.
	DEFAULT: 1
N_CASES : int
	Number of boxes in each direction to compute the displacement grid.
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
R_MIN [PLOT or SHOW mode] : float
	Minimum wave length norm for plots.
	DEFAULT: active_particles.analysis.ctt._r_min
R_MAX [PLOT or SHOW mode] : float
	Maximum wave length norm for plots.
	DEFAULT: sqrt(2) times BOX_SIZE
Y_MIN [PLOT or SHOW mode] : float
	Minimum y-coordinate for plots.
	DEFAULT: fit to data
Y_MAX [PLOT or SHOW mode] : float
	Maximum y-coordinate for plots.
	DEFAULT: fit to data
SLOPE [FITTING_LINE mode] : float
	Initial slope for fitting line.
	DEFAULT: active_particles.analysis.ctt._slope0
SLOPE_MIN [FITTING_LINE mode] : float
	Minimum slope for fitting line.
	DEFAULT: active_particles.analysis.ctt._slope_min
SLOPE_MAX [FITTING_LINE mode] : float
	Maximum slope for fitting line.
	DEFAULT: active_particles.analysis.ctt._slope_max
R_MAX_CSS [STRAIN_CORRELATIONS mode] : float
	Maximum radius in infinite norm for strain correlations plot.
	DEFAULT: active_particles.analysis.css._r_max

Output
------
[COMPUTE MODE]
> Prints execution time.
> Saves wave vectors grid, 2D grid and 1D cylindrical average of mean squared
cross products of normalised wave vectors and displacement Fourier transform
according to active_particles.naming.Ctt standards in DATA_DIRECTORY.
> Saves wave vectors grid, 2D grid and 1D cylindrical average of mean squared
dot products of normalised wave vectors and displacement Fourier transform
according to active_particles.naming.Cll standards in DATA_DIRECTORY.
[SHOW or PLOT mode]
> Plots cylindrical averages as functions of wave length.
[SAVE mode]
> Saves figure in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output
from active_particles.dat import Dat, Gsd
from active_particles.maths import g2Dto1Dgrid, kFFTgrid, wave_vectors_2D,\
	FFT2Dfilter, divide_arrays

from active_particles.analysis.cuu import displacement_grid
from active_particles.analysis.css import _r_max as _r_max_css
from active_particles.analysis.correlations import CorGrid
from active_particles.plot.mpl_tools import FittingLine

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

from math import ceil

import numpy as np

import pickle

from collections import OrderedDict

from datetime import datetime

import matplotlib as mpl
if not(get_env('SHOW', default=False, vartype=bool)):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider

def plot_product(product, ax):
	"""
	Plot cross or product data on ax, on log-log axis, as a function of
	wave length.

	Parameters
	----------
	product : Numpy array
		Cross or dot product.
	ax : matplotlib axis
		Axis to plot on.

	Returns
	-------
	line : matplotlib.lines.Line2D
		Plotted product.
	"""

	line, = ax.loglog(2*np.pi/product[1:, 0], product[1:, 1])
	return line

def plot_cross(ax):
	"""
	Plot cross product data on ax.

	Parameters
	----------
	ax : matplotlib axis
		Axis to plot on.

	Returns
	-------
	cross : matplotlib.lines.Line2D
		Plotted cross product.
	"""

	cross = plot_product(k_cross_FFTugrid1D_sqnorm, ax)
	cross.set_color('blue')
	cross.set_label(
		r'$\left<||\vec{k}\wedge\tilde{\vec{u}}(\vec{k})||^2\right>/k^2$')

	return cross

def plot_dot(ax):
	"""
	Plot dot product data on ax.

	Parameters
	----------
	ax : matplotlib axis
		Axis to plot on.

	Returns
	-------
	dot : matplotlib.lines.Line2D
		Plotted dot product.
	"""

	dot = plot_product(k_dot_FFTugrid1D_sqnorm, ax)
	dot.set_color('orange')
	dot.set_label(
		r'$\left<||\vec{k}\cdot\tilde{\vec{u}}(\vec{k})||^2\right>/k^2$')

	return dot

class StrainCorrelations:
	"""
	Manipulate and plot strain correlations, computed from dot and cross
	products of normalised wave vector and displacement Fourier transform.
	"""

	def __init__(self, wave_vectors,
	    k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm):
		"""
		Parameters
		----------
		wave_vectors : (_, _, 2) array-like
			Wave vectors at which products were calculated.
		k_cross_FFTugrid2D_sqnorm : 2D array-like
			Square norm of cross product of normalised wave vector and
			displacement Fourier transform.
		k_dot_FFTugrid2D_sqnorm : 2D array-like
			Square norm of dot product of normalised wave vector and
			displacement Fourier transform.

		"""

		self.wave_vectors = wave_vectors
		self.kxsq = np.array(self.wave_vectors)[:, :, 0]**2	# squared wave vectors x-coordinates
		self.kysq = np.array(self.wave_vectors)[:, :, 1]**2	# squared wave vectors y-coordinates
		self.ksq = self.kxsq + self.kysq					# squared wave vectors norm

		self.k_cross_FFTugrid2D_sqnorm = k_cross_FFTugrid2D_sqnorm
		self.k_dot_FFTugrid2D_sqnorm = k_dot_FFTugrid2D_sqnorm

		self.strain_correlations_FFT = (
			- divide_arrays(
				(self.k_cross_FFTugrid2D_sqnorm - self.k_dot_FFTugrid2D_sqnorm)
				*self.kxsq*self.kysq,
				self.ksq)
			+ self.k_cross_FFTugrid2D_sqnorm*self.ksq/4)	# Fourier transform of strain correlations

	def strain_correlations(self, r_cut=0):
		"""
		Compute strain correlations from inverse fast Fourier transform of
		self.strain_correlations_FFT with low wave lengths cut at r_cut.

		Parameters
		----------
		r_cut : float
			Wave length cut-off radius, equivalent to coarse-graining cut-off
			radius.

		Returns
		-------
		sc : Numpy array
			Strain correlations.
		"""

		sc = (FFT2Dfilter(self.strain_correlations_FFT,
			wave_vectors=self.wave_vectors)
			.cut_low_wave_lengths(r_cut)
			.get_signal()).real
		return sc/sc[0, 0]	# correlation normalisation

	def plot(self, box_size, r_max_css):
		"""
		Plots strain correlations with slider for cut-off radius.

		Parameters
		----------
		box_size : float
			System length.
		r_max_css : float
			Maximum radius in infinite norm for strain correlations plots.
		"""

		self.box_size = box_size
		self.r_max_css = r_max_css

		self.fig, (self.ax_sc, self.ax_kFFTugrid) = plt.subplots(1, 2)

		self.r_cut = 0	# cut-off radius

		# STRAIN CORRELATIONS

		grid = self.strain_correlations_corgrid()	# correlation grid

		Cmin = np.min(grid.grid)
		Cmax = np.max(grid.grid)
		cmap = plt.cm.jet
		CvNorm = colors.Normalize(vmin=Cmin, vmax=Cmax)

		self.grid_plot = self.ax_sc.imshow(grid.display_grid.grid,
			cmap=cmap, norm=CvNorm, extent=
			[-self.r_max_css, self.r_max_css, -self.r_max_css, self.r_max_css])

		self.colormap_ax = make_axes_locatable(self.ax_sc).append_axes(
		    'right', size='5%', pad=0.05)			# color map axes
		self.colormap = mpl.colorbar.ColorbarBase(self.colormap_ax, cmap=cmap,
		    norm=CvNorm, orientation='vertical')    # color map

		# CROSS AND DOT PRODUCTS

		self.ax_kFFTugrid.set_xlim([r_min, r_max])
		self.ax_kFFTugrid.set_ylim([y_min, y_max])

		cross = plot_cross(self.ax_kFFTugrid)
		dot = plot_dot(self.ax_kFFTugrid)

		self.r_cut_line = self.ax_kFFTugrid.axvline(x=self.r_cut, color='red',
			label=r'$r_{cut}$')	# vertical line delimiting wave lengths cutting domain

		self.ax_kFFTugrid.legend()	# add legend to plot

		self.r_cut_line.figure.canvas.mpl_connect('button_press_event',
		    self.update_r_cut_line)	# call self.update_r_cut_line() on button press event

		self.slider_ax = make_axes_locatable(self.ax_kFFTugrid).append_axes(
		    'bottom', size='5%', pad=0.6)			# slider Axes
		self.slider = Slider(self.slider_ax, r'$r_{cut}$', self.r_cut,
			self.r_max_css, valinit=0)				# slider
		self.slider.on_changed(self.update_r_cut)	# call update_r_cut when slider value is changed

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

		sc = self.strain_correlations(r_cut=self.r_cut)	# strain correlations grid
		return CorGrid(sc, self.box_size, display_size=2*self.r_max_css)

	def update_r_cut(self, val):
		"""
		Updates cut-off radius on slider change.
		"""

		self.r_cut=self.slider.val	# new cut-off radius
		self.draw()					# updates figure

	def update_r_cut_line(self, event):
		"""
		Updates cut-off radius on click on figure.
		"""

		if event.inaxes != self.r_cut_line.axes: return	# if Axes instance mouse is over is different than cut-off line figure Axes

		self.r_cut = event.xdata		# new cut-off radius
		self.slider.set_val(self.r_cut)	# updates slider value

		self.draw() # updates figure

	def draw(self):
		"""
		Updates line position and strain correlations grid according to cut-off
		radius self.r_cut.
		"""

		grid = self.strain_correlations_corgrid()		# new grid
		self.grid_plot.set_data(grid.display_grid.grid)	# plot grid
		self.grid_plot.figure.canvas.draw()				# update plot

		self.r_cut_line.set_xdata([self.r_cut]*2)	# update line position
		self.r_cut_line.figure.canvas.draw()		# update plot

def plot():
	"""
	Plots mean square norm of dot and cross products of wave vector and
	displacement Fourier transform.

	Returns
	-------
	fig : matplotlib figure
		Figure.
	ax_cross : matplotlib axis
		Mean square cross product axis.
	ax_dot : matplotlib axis
		Mean square dot product axis.
	ax_super : matplotlib axis
		Superimposed mean square cross and dot products axis.
	fl_cross [FITTING_LINE mode] : active_particles.plot.mpl_tools.FittingLine
		Fitting line object for cross axis.
	fl_dot [FITTING_LINE mode] : active_particles.plot.mpl_tools.FittingLine
		Fitting line object for dot axis.
		NOTE: Only in FITTING_LINE mode.
	fl_super [FITTING_LINE mode] : active_particles.plot.mpl_tools.FittingLine
		Fitting line object for superimposition axis.
		NOTE: Only in FITTING_LINE mode.
	"""

	to_return = ()	# variables to return

	fig = plt.figure()
	to_return += (fig,)

	fig.set_size_inches(30, 30)
	fig.subplots_adjust(wspace=0.3)
	fig.subplots_adjust(hspace=0.3)

	fig.suptitle(r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
		% (parameters['N'], parameters['density'], parameters['vzero'],
		parameters['dr']) + '\n' +
		r'$S_{init}=%.2e, \Delta t=%.2e$' % (init_frame,
		dt*parameters['period_dump']*parameters['time_step']) +
		r'$, S_{max}=%.2e, N_{cases}=%.2e$' % (int_max, Ncases))

	gs = GridSpec(2, 2)

	# CROSS

	ax_cross = plt.subplot(gs[0, 0])
	to_return += (ax_cross,)
	cross_cross = plot_cross(ax_cross)	# cross product

	ax_cross.set_xlabel(r'$\lambda = 2\pi/k$')
	ax_cross.set_xlim([r_min, r_max])
	ax_cross.set_ylabel(cross_cross.get_label())
	ax_cross.set_ylim([y_min, y_max])

	if get_env('FITTING_LINE', default=False, vartype=bool):	# FITTING_LINE mode
		fl_cross = FittingLine(ax_cross, slope0, slope_min, slope_max,
			x_fit='\lambda')									# add fitting line to plot
		to_return += (fl_cross,)

	# DOT

	ax_dot = plt.subplot(gs[1, 0])
	to_return += (ax_dot,)
	dot_dot = plot_dot(ax_dot)	# dot product

	ax_dot.set_xlabel(r'$\lambda = 2\pi/k$')
	ax_dot.set_xlim([r_min, r_max])
	ax_dot.set_ylabel(dot_dot.get_label())
	ax_dot.set_ylim([y_min, y_max])

	if get_env('FITTING_LINE', default=False, vartype=bool):	# FITTING_LINE mode
		fl_dot = FittingLine(ax_dot, slope0, slope_min, slope_max,
			x_fit='\lambda')									# add fitting line to plot
		to_return += (fl_dot,)

	# SUPER

	ax_super = plt.subplot(gs[:, 1])
	to_return += (ax_super,)
	cross_super = plot_cross(ax_super)	# cross product
	dot_super = plot_dot(ax_super)		# dot product

	legend = ax_super.legend()
	ax_super.add_artist(legend)

	ax_super.set_xlabel(r'$\lambda = 2\pi/k$')
	ax_super.set_xlim([r_min, r_max])
	ax_super.set_ylabel(r'$S(k)$')
	ax_super.set_ylim([y_min, y_max])

	if get_env('FITTING_LINE', default=False, vartype=bool):	# FITTING_LINE mode
		fl_super = FittingLine(ax_super, slope0, slope_min, slope_max,
			x_fit='\lambda')									# add fitting line to plot
		to_return += (fl_super,)

	# SAVING

	if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
		image_name, = naming_Ctt.image().filename(**attributes)
		fig.savefig(joinpath(data_dir, image_name))

	# STRAIN CORRELATIONS

	if get_env('STRAIN_CORRELATIONS', default=False, vartype=bool):	# STRAIN_CORRELATIONS mode

		sc = StrainCorrelations(wave_vectors,
		    k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm)
		to_return += (sc,)
		sc.plot(parameters['box_size'], r_max_css)

		sc.fig.set_size_inches(16, 16)
		sc.fig.subplots_adjust(wspace=0.3)
		sc.fig.subplots_adjust(hspace=0.3)

		sc.fig.suptitle(r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
			% (parameters['N'], parameters['density'], parameters['vzero'],
			parameters['dr']) + '\n' +
			r'$S_{init}=%.2e, \Delta t=%.2e$' % (init_frame,
			dt*parameters['period_dump']*parameters['time_step']) +
			r'$, S_{max}=%.2e, N_{cases}=%.2e$' % (int_max, Ncases))

		sc.ax_sc.set_title('2D ' + r'$C_{\varepsilon_{xy}\varepsilon_{xy}}$')
		sc.ax_sc.set_xlabel(r'$x$')
		sc.ax_sc.set_ylabel(r'$y$')
		sc.colormap.set_label(r'$C_{\varepsilon_{xy}\varepsilon_{xy}}$',
			labelpad=20, rotation=270)

		sc.ax_kFFTugrid.set_xlabel(ax_super.get_xlabel())
		sc.ax_kFFTugrid.set_xlim(ax_super.get_xlim())
		sc.ax_kFFTugrid.set_ylabel(ax_super.get_ylabel())
		sc.ax_kFFTugrid.set_ylim(ax_super.get_ylim())

	# RETURN

	return to_return

# DEFAULT VARIABLES

_r_min = 1		# default minimum wave length for plots

_slope0 = 2		# default initial slope for fitting line
_slope_min = 0	# default minimum slope for fitting line
_slope_max = 4	# default maximum slope for fitting line

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    dt = get_env('TIME', default=-1, vartype=int)	# lag time for displacement

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of intervals of length dt considered for the calculation

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
	       parameters = pickle.load(param_file)				# parameters hash table

    box_size = get_env('BOX_SIZE', default=parameters['box_size'],
		vartype=float)									# size of the square box to consider
    centre = (get_env('X_ZERO', default=0, vartype=float),
		get_env('Y_ZERO', default=0, vartype=float))	# centre of the box

    prep_frames = ceil(parameters['prep_steps']/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

    Ncases = get_env('N_CASES', default=ceil(np.sqrt(parameters['N'])),
		vartype=int)        # number of boxes in each direction to compute the displacement grid
    dL = box_size/Ncases    # boxes separation

    Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame
    Nframes = Nentries - init_frame									# number of frames available for the calculation

    dt = Nframes + dt if dt <= 0 else dt	# length of the interval of time for which displacements are calculated

    attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'init_frame': init_frame, 'dt': dt,
		'int_max': int_max, 'Ncases': Ncases, 'box_size': box_size,
        'x_zero': centre[0], 'y_zero': centre[1]}		# attributes displayed in filenames
    naming_Ctt = naming.Ctt()                           # Ctt naming object
    Ctt_filename, = naming_Ctt.filename(**attributes)   # Ctt filename
    naming_Cll = naming.Cll()                           # Cll naming object
    Cll_filename, = naming_Cll.filename(**attributes)   # Cll filename

	# STANDARD OUTPUT

    if 'SLURM_JOB_ID' in envvar:	# script executed from Slurm job scheduler
        slurm_output(joinpath(data_dir, 'out'), naming_Ctt, attributes)

    # MODE SELECTION

    if get_env('COMPUTE', default=False, vartype=bool):	# COMPUTE mode

        startTime = datetime.now()

		# VARIABLE DEFINITIONS

        wrap_file_name = get_env('WRAPPED_FILE',
			default=joinpath(data_dir, naming.wrapped_trajectory_file))		# wrapped trajectory file (.gsd)
        unwrap_file_name = get_env('UNWRAPPED_FILE',
			default=joinpath(data_dir, naming.unwrapped_trajectory_file))	# unwrapped trajectory file (.dat)

        times = np.array(list(OrderedDict.fromkeys(map(
			lambda x: int(x),
			np.linspace(init_frame, Nentries - dt - 1, int_max)
			))))	# frames at which shear strain will be calculated

        # DISPLACEMENT CORRELATIONS

        with open(wrap_file_name, 'rb') as wrap_file,\
			open(unwrap_file_name, 'rb') as unwrap_file:	# opens wrapped and unwrapped trajectory files

            w_traj = Gsd(wrap_file, prep_frames=prep_frames)	# wrapped trajectory object
            u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object
            Ugrid = list(map(
                lambda time: displacement_grid(parameters['box_size'],
                box_size, centre, Ncases, time, dt, w_traj, u_traj, dL),
                times))											# lists of displacement variables

        wave_vectors = wave_vectors_2D(Ncases, Ncases, d=dL)			# wave vectors grid
        wave_vectors_norm = np.sqrt(np.sum(wave_vectors**2, axis=-1))	# wave vectors norm grid

        k_cross_FFTUgrid, k_dot_FFTUgrid = np.transpose(
			list(map(lambda grid: kFFTgrid(grid), Ugrid)),
			(1, 0, 2, 3))	# wave vector grids, and cross

        k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm = list(map(
			lambda Grid: np.mean(
			list(map(lambda grid: (np.conj(grid)*grid).real, Grid)),
			axis=0),
			[k_cross_FFTUgrid, k_dot_FFTUgrid]))	# grid of mean square norms of cross and dot products of normalised wave vectors with displacement grids Fourier transform

        k_cross_FFTugrid1D_sqnorm, k_dot_FFTugrid1D_sqnorm = list(map(
			lambda grid2D: g2Dto1Dgrid(grid2D, wave_vectors_norm),
			[k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm]))	# cylindrical averages of mean square norms of cross and dot products of normalised wave vectors with displacement grids Fourier transform

        # SAVING

        with open(joinpath(data_dir, Ctt_filename), 'wb') as Ctt_dump_file,\
			open(joinpath(data_dir, Cll_filename), 'wb') as Cll_dump_file:
            pickle.dump([wave_vectors, k_cross_FFTugrid2D_sqnorm,
				k_cross_FFTugrid1D_sqnorm], Ctt_dump_file)
            pickle.dump([wave_vectors, k_dot_FFTugrid2D_sqnorm,
				k_dot_FFTugrid1D_sqnorm], Cll_dump_file)

        # EXECUTION TIME

        print("Execution time: %s" % (datetime.now() - startTime))

    if get_env('PLOT', default=False, vartype=bool):	# PLOT mode

		# DATA

        with open(joinpath(data_dir, Ctt_filename), 'rb') as Ctt_dump_file,\
			open(joinpath(data_dir, Cll_filename), 'rb') as Cll_dump_file:
            (wave_vectors, k_cross_FFTugrid2D_sqnorm,
				k_cross_FFTugrid1D_sqnorm) = pickle.load(Ctt_dump_file)
            (_, k_dot_FFTugrid2D_sqnorm,
				k_dot_FFTugrid1D_sqnorm) = pickle.load(Cll_dump_file)

    if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

		# PLOT

        r_min = get_env('R_MIN', default=_r_min, vartype=float)					# minimum wave length for plots
        r_max = get_env('R_MAX', default=np.sqrt(2)*box_size, vartype=float)	# maximum wave length for plots

        y_min = get_env('Y_MIN',
			default=np.min((k_cross_FFTugrid1D_sqnorm[1:, 1],
			k_dot_FFTugrid1D_sqnorm[1:, 1])),
			vartype=float)	# minimum y-coordinate for plots
        y_max = get_env('Y_MAX',
			default=np.max((k_cross_FFTugrid1D_sqnorm[1:, 1],
			k_dot_FFTugrid1D_sqnorm[1:, 1])),
			vartype=float)	# maximum y-coordinate for plots

        slope0 = get_env('SLOPE', default=_slope0, vartype=float)			# initial slope for fitting line
        slope_min = get_env('SLOPE_MIN', default=_slope_min, vartype=float)	# minimum slope for fitting line
        slope_max = get_env('SLOPE_MAX', default=_slope_max, vartype=float)	# maximum slope for fitting line

        r_max_css = get_env('R_MAX_CSS', default=_r_max_css, vartype=float)	# maximum radius in infinite norm for strain correlations plot

        plot = plot()

        if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
            plt.show()
