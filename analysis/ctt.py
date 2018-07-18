"""
Module ctt calculates or plots mean squared cross and dot products of
normalised wave vectors and displacement Fourier transform and strain
correlations from these variables.

Files are saved according to active_particles.naming.Ctt (cross product) and
active_particles.naming.Cll (dot product) naming standards.

A brief description of the algorithm can be found at:
https://yketa.github.io/UBC_2018_Wiki/#Collective%20mean%20square%20displacements

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
R_CUT_FOURIER [PLOT or SHOW mode] : float
	Initial wave length Gaussian cut-off radius.
	DEFAULT: active_particles.analysis.css._r_cut_fourier
SMOOTH [PLOT or SHOW mode] : float
	C44 Gaussian smoothing length scale.
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
R_MAX_CSS [PLOT mode] : float
	Maximum radius in infinite norm for strain correlations plot.
	DEFAULT: active_particles.analysis.css._r_max
POINTS_X_C44 [PLOT mode] : int
	Number of radii at which to evaluate integrated strain correlation.
	DEFAULT: active_particles.analysis.css._points_x_c44
POINTS_THETA_C44 [PLOT mode] : int
	Number of angles to evaluate integrated strain correlation.
	DEFAULT: active_particles.analysis.css._points_theta_c44
Y_MIN_C44 [PLOT mode] : float
	Minimum plot value for C44.
	DEFAULT: active_particles.analysis.css._y_min_c44
Y_MAX_C44 [PLOT mode] : float
	Maximum plot value for C44.
	DEFAULT: active_particles.analysis.css._y_max_c44
R_MIN_C44 [PLOT mode] : float
	Minimum radius in average particle separation for C44 calculation.
	DEFAULT: active_particles.analysis.css._r_min_c44
R_MAX_C44 [PLOT mode] : float
	Maximum radius in average particle separation for C44 calculation.
	DEFAULT: active_particles.analysis.css._r_max_c44
SLOPE_C44 [FITTING_LINE mode] : slope
	Initial slope for fitting line.
	DEFAULT: active_particles.analysis.css._slope0_c44
SLOPE_MIN_C44 [FITTING_LINE mode] : float
	Minimum slope for fitting line.
	DEFAULT: active_particles.analysis.css._slope_min_c44
SLOPE_MAX_C44 [FITTING_LINE mode] : float
	Maximum slope for fitting line.
	DEFAULT: active_particles.analysis.css._slope_max_c44

Output
------
[COMPUTE mode]
> Prints execution time.
> Saves wave vectors grid, 2D grid and 1D cylindrical average of mean squared
cross products of normalised wave vectors and displacement Fourier transform
according to active_particles.naming.Ctt standards in DATA_DIRECTORY.
> Saves wave vectors grid, 2D grid and 1D cylindrical average of mean squared
dot products of normalised wave vectors and displacement Fourier transform
according to active_particles.naming.Cll standards in DATA_DIRECTORY.
[SHOW or PLOT mode]
> Plots cylindrical averages as functions of wave length and resulting
strain correlations.
[SAVE mode]
> Saves collective mean square displacement figure in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output, linframes
from active_particles.dat import Dat, Gsd
from active_particles.maths import g2Dto1Dgrid, kFFTgrid, wave_vectors_2D,\
	divide_arrays, FFT2Dfilter
from active_particles.quantities import nD0_active

from active_particles.analysis.cuu import displacement_grid
from active_particles.analysis.css import StrainCorrelations, Css2DtoC44,\
	_r_max as _r_max_css, _c_min, _c_max, _slope0_c44,\
	_slope_min_c44, _slope_max_c44, _points_x_c44, _points_theta_c44,\
	_y_min_c44, _y_max_c44, _r_min_c44, _r_max_c44, _r_cut_fourier
from active_particles.analysis.correlations import CorGrid
from active_particles.plot.mpl_tools import FittingLine, GridCircle
from active_particles.analysis.number import count_particles

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

from math import ceil

import numpy as np

import pickle

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

# DEFAULT VARIABLES

_r_min = 5e-1	# default minimum wave length for plots

_slope0 = 2		# default initial slope for fitting line
_slope_min = 0	# default minimum slope for fitting line
_slope_max = 4	# default maximum slope for fitting line

_cross_label =\
	r'$\left<||\vec{k}\wedge\tilde{\vec{u}}(\vec{k})||^2\right>/k^2$'	# transversal collective mean square displacement label
_dot_label =\
	r'$\left<||\vec{k}\cdot\tilde{\vec{u}}(\vec{k})||^2\right>/k^2$'	# longitudinal collective mean square displacement label

# FUNCTIONS AND CLASSES

def plot_product(product, av_p_sep, ax):
	"""
	Plot cross or product data on ax, on log-log axis, as a function of
	wave length.

	Parameters
	----------
	product : Numpy array
		Cross or dot product.
	av_p_sep : float
		Average particle separation.
	ax : matplotlib axis
		Axis to plot on.

	Returns
	-------
	line : matplotlib.lines.Line2D
		Plotted product.
	"""

	line, = ax.loglog(2*np.pi/(product[1:, 0]*av_p_sep), product[1:, 1])
	return line

def plot_cross(k_cross_FFTugrid1D_sqnorm, av_p_sep, ax):
	"""
	Plot cross product data on ax.

	Parameters
	----------
	k_cross_FFTugrid1D_sqnorm : 2D array-like
		Cylindrical average of the square norm of cross product of normalised
		wave vector and displacement Fourier transform.
	av_p_sep : float
		Average particle separation.
	ax : matplotlib axis
		Axis to plot on.

	Returns
	-------
	cross : matplotlib.lines.Line2D
		Plotted cross product.
	"""

	cross = plot_product(k_cross_FFTugrid1D_sqnorm, av_p_sep, ax)
	cross.set_color('blue')
	cross.set_label(_cross_label)

	return cross

def plot_dot(k_dot_FFTugrid1D_sqnorm, av_p_sep, ax):
	"""
	Plot dot product data on ax.

	Parameters
	----------
	k_dot_FFTugrid1D_sqnorm : 2D array-like
		Cylindrical average of the square norm of dot product of normalised
		wave vector and displacement Fourier transform.
	av_p_sep : float
		Average particle separation.
	ax : matplotlib axis
		Axis to plot on.

	Returns
	-------
	dot : matplotlib.lines.Line2D
		Plotted dot product.
	"""

	dot = plot_product(k_dot_FFTugrid1D_sqnorm, av_p_sep, ax)
	dot.set_color('orange')
	dot.set_label(_dot_label)

	return dot

class StrainCorrelationsCMSD(StrainCorrelations):
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
		self.ksq = self.kxsq + self.kysq					# squared wave vectors norms
		self.k = np.sqrt(self.ksq)							# wave vectors norms

		self.k_cross_FFTugrid2D_sqnorm = k_cross_FFTugrid2D_sqnorm
		self.k_cross_FFTugrid1D_sqnorm = g2Dto1Dgrid(
			k_cross_FFTugrid2D_sqnorm, self.k)

		self.k_dot_FFTugrid2D_sqnorm = k_dot_FFTugrid2D_sqnorm
		self.k_dot_FFTugrid1D_sqnorm = g2Dto1Dgrid(
			k_dot_FFTugrid2D_sqnorm, self.k)

		self.cross_label = _cross_label.replace('$', '')
		self.dot_label = _dot_label.replace('$', '')

	def get_strain_correlations_FFT(self,
		k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm):
		"""
		Returns Fourier transform of strain correlations from square norm of
		cross product of normalised wave vector and displacement Fourier
		transform and square norm of dot product of normalised wave vector and
		displacement Fourier transform.

		Parameters
		----------
		k_cross_FFTugrid2D_sqnorm : 2D array-like
			Square norm of cross product of normalised wave vector and
			displacement Fourier transform.
		k_dot_FFTugrid2D_sqnorm : 2D array-like
			Square norm of dot product of normalised wave vector and
			displacement Fourier transform.

		Returns
		-------
		strain_correlations_FFT : 2D array-like
			Fourier transform of strain correlations
		"""

		return (
			- divide_arrays(
				(k_cross_FFTugrid2D_sqnorm - k_dot_FFTugrid2D_sqnorm)
				*self.kxsq*self.kysq,
				self.ksq)
			+ k_cross_FFTugrid2D_sqnorm*self.ksq/4)

	def strain_correlations(self, r_cut=0):
		"""
		Computes strain correlations from inverse fast Fourier transform of
		self.strain_correlations_FFT, calculated from transversal and
		longitudinal mean square displacements with low wave lengths Gaussian
		cut at r_cut.

		Parameters
		----------
		r_cut : float
			Wave length Gaussian cut-off radius, equivalent to coarse-graining
			cut-off radius.

		Returns
		-------
		sc : Numpy array
			Strain correlations.
		"""

		self.filtered_k_cross_FFTugrid2D_sqnorm = (FFT2Dfilter(
			self.k_cross_FFTugrid2D_sqnorm, wave_vectors=self.wave_vectors)
			.gaussian_filter(r_cut)
			.signalFFT)	# square norm of cross product of normalised wave vector and displacement Fourier transform with low wave lengths Gaussian cut at r_cut
		self.filtered_k_cross_FFTugrid1D_sqnorm = g2Dto1Dgrid(
			self.filtered_k_cross_FFTugrid2D_sqnorm, self.k)

		self.filtered_k_dot_FFTugrid2D_sqnorm = (FFT2Dfilter(
			self.k_dot_FFTugrid2D_sqnorm, wave_vectors=self.wave_vectors)
			.gaussian_filter(r_cut)
			.signalFFT)	# square norm of dot product of normalised wave vector and displacement Fourier transform with low wave lengths Gaussian cut at r_cut
		self.filtered_k_dot_FFTugrid1D_sqnorm = g2Dto1Dgrid(
			self.filtered_k_dot_FFTugrid2D_sqnorm, self.k)

		self.strain_correlations_FFT = self.get_strain_correlations_FFT(
			self.filtered_k_cross_FFTugrid2D_sqnorm,
			self.filtered_k_dot_FFTugrid2D_sqnorm)

		return super().strain_correlations(r_cut=0)

	def plot(self, box_size, r_max_css, av_p_sep, r_min, r_max, y_min, y_max,
		points_x_c44=_points_x_c44, points_theta_c44=_points_theta_c44,
		y_min_c44=_y_min_c44, y_max_c44=_y_max_c44,
		r_min_c44=_r_min_c44, r_max_c44=_r_max_c44,
		r_cut=0, smooth=0):
		"""
		Plots collective mean square displacements strain correlations with
		slider for cut-off radius.

		Parameters
		----------
		box_size : float
			System length.
		r_max_css : float
			Maximum radius in infinite norm for strain correlations plots.
		av_p_sep : float
			Average particle separation.
		r_min : float
			Minimum wave length for collective mean square displacements.
		r_max : float
			Maximum wave length for collective mean square displacements.
		y_min : float
			Minimum value collective mean square displacements.
		y_max : float
			Maximum value collective mean square displacements.
		points_x_c44 : int
            Number of radii at which to evaluate integrated strain correlation.
			(default: active_particles.plot.c44._points_x)
        points_theta_c44 : int
            Number of angles to evaluate integrated strain correlation.
			(default: active_particles.plot.c44._points_theta)
		y_min_c44 : float
			Minimum plot value for C44.
			(default: active_particles.plot.c44._y_min)
		y_max_c44 : float
			Maximum plot value for C44.
			(default: active_particles.plot.c44._y_max)
        r_min_c44 : float
            Minimum radius in average particle separation for C44 calculation.
			(default: active_particles.plot.c44._r_min)
        r_max_c44 : float
            Maximum radius in average particle separation for C44 calculation.
			(default: active_particles.plot.c44._r_max)
		r_cut : float
			Initial wave length Gaussian cut-off radius. (default: 0)
		smooth : float
			C44 Gaussian smoothing length scale. (default: 0)
		"""

		self.cor_name = 'C_{\\varepsilon_{xy}\\varepsilon_{xy}}'	# name of plotted correlation

		self.box_size = box_size
		self.av_p_sep = av_p_sep

		self.r_max_css = r_max_css

		self.r_cut = r_cut

		grid = self.strain_correlations_corgrid()	# correlation grid

		# COLECTIVE MEAN SQUARE DISPLACEMENTS FIGURE

		self.fig_cmsd = plt.figure()
		gs = GridSpec(2, 2)

		self.ax_cross = plt.subplot(gs[0, 0])
		self.cross_line = plot_cross(self.filtered_k_cross_FFTugrid1D_sqnorm,
			self.av_p_sep, self.ax_cross)

		self.ax_dot = plt.subplot(gs[1, 0])
		self.dot_line = plot_dot(self.filtered_k_dot_FFTugrid1D_sqnorm,
			self.av_p_sep, self.ax_dot)

		self.ax_super = plt.subplot(gs[:, 1])
		self.super_cross_line = plot_cross(
			self.filtered_k_cross_FFTugrid1D_sqnorm, self.av_p_sep,
			self.ax_super)
		self.super_dot_line = plot_dot(
			self.filtered_k_dot_FFTugrid1D_sqnorm, self.av_p_sep,
			self.ax_super)

		self.ax_super.set_title(r'$r_{cut}/a = %.2e$' % self.r_cut)

		# STRAIN CORRELATIONS AND SUPERIMPOSED ORIGINAL CMSD FIGURE

		self.fig_sc, (self.ax_sc, self.ax_kFFTugrid) = plt.subplots(1, 2)

		Cmin = np.max((np.min(grid.grid), _c_min))
		Cmax = np.max((np.min(grid.grid), _c_max))
		cmap = plt.cm.jet
		CvNorm = colors.Normalize(vmin=Cmin, vmax=Cmax)

		self.grid_plot = self.ax_sc.imshow(grid.display_grid.grid,
			cmap=cmap, norm=CvNorm, extent=
			[-self.r_max_css, self.r_max_css, -self.r_max_css, self.r_max_css])

		self.colormap_ax = make_axes_locatable(self.ax_sc).append_axes(
		    'right', size='5%', pad=0.05)			# color map axes
		self.colormap = mpl.colorbar.ColorbarBase(self.colormap_ax, cmap=cmap,
		    norm=CvNorm, orientation='vertical')    # color map

		self.ax_kFFTugrid.set_xlim([r_min/self.av_p_sep, r_max/self.av_p_sep])
		self.ax_kFFTugrid.set_ylim([y_min, y_max])

		_ = plot_cross(self.k_cross_FFTugrid1D_sqnorm, self.av_p_sep,
			self.ax_kFFTugrid)
		_ = plot_dot(self.k_dot_FFTugrid1D_sqnorm, self.av_p_sep,
			self.ax_kFFTugrid)

		self.r_cut_line = self.ax_kFFTugrid.axvline(x=self.r_cut,
			color='red', label=r'$r_{cut}/a$')	# vertical line delimiting wave lengths cutting domain

		self.r_cut_line.figure.canvas.mpl_connect('button_press_event',
		    self.update_r_cut_line)	# call self.update_r_cut_line() on button press event

		self.slider_ax = make_axes_locatable(self.ax_kFFTugrid).append_axes(
		    'bottom', size='5%', pad=0.6)						# slider Axes
		self.slider = Slider(self.slider_ax, r'$r_{cut}/a$', 0,
			self.r_max_css/self.av_p_sep, valinit=self.r_cut)	# slider
		self.slider.on_changed(self.update_r_cut)				# call update_r_cut when slider value is changed

		# GRID CIRLCE FIGURE

		self.grid_circle = GridCircle(grid.display_grid.grid, extent=
			[-self.r_max_css, self.r_max_css, -self.r_max_css, self.r_max_css],
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
		self.smooth = smooth

		self.toC44 = Css2DtoC44(self.box_size,
			self.points_x_c44, self.points_theta_c44,
			self.av_p_sep*self.r_min_c44, self.av_p_sep*self.r_max_c44)

		self.fig_c44, self.ax_c44 = plt.subplots()
		self.ax_c44.set_title(r'$r_{cut}/a = %.2e$' % self.r_cut)
		self.ax_c44.set_xlim([self.r_min_c44, self.r_max_c44])
		self.ax_c44.set_ylim([self.y_min_c44, self.y_max_c44])
		self.ax_c44.set_yscale('log')
		self.ax_c44.set_xscale('log')

		self.c44 = self.css2Dtoc44(grid.grid)	# list of [r, C44(r)]
		self.line_c44, = self.ax_c44.plot(self.c44[:, 0]/self.av_p_sep,
			self.c44[:, 1])						# C44 plotting line

	# METHODS CALLED BY SELF.PLOT()

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

		self.ax_super.set_title(r'$r_{cut}/a = %.2e$' % self.r_cut)				# update title
		self.grid_circle.ax_grid.set_title(
			'2D ' + r'$%s(r_{cut}/a = %.2e)$' % (self.cor_name, self.r_cut))	# update title
		self.ax_c44.set_title(r'$r_{cut}/a = %.2e$' % self.r_cut)				# update title

		grid = self.strain_correlations_corgrid()		# new grid

		self.cross_line.set_ydata(
			self.filtered_k_cross_FFTugrid1D_sqnorm[1:, 1])	# update TCMSD values
		self.ax_cross.relim()								# recompute data limits
		self.ax_cross.autoscale_view()						# scale axis to data
		self.dot_line.set_ydata(
			self.filtered_k_dot_FFTugrid1D_sqnorm[1:, 1])	# update LCMSD values
		self.ax_dot.relim()									# recompute data limits
		self.ax_dot.autoscale_view()						# scale axis to data
		self.super_cross_line.set_ydata(
			self.filtered_k_cross_FFTugrid1D_sqnorm[1:, 1])	# update TCMSD values
		self.super_dot_line.set_ydata(
			self.filtered_k_dot_FFTugrid1D_sqnorm[1:, 1])	# update LCMSD values
		self.ax_super.relim()								# recompute data limits
		self.ax_super.autoscale_view()						# scale axis to data
		self.cross_line.figure.canvas.draw()				# update plot

		self.grid_plot.set_data(grid.display_grid.grid)	# plot grid
		self.grid_plot.figure.canvas.draw()				# update plot

		self.r_cut_line.set_xdata([self.r_cut]*2)	# update line position
		self.r_cut_line.figure.canvas.draw()		# update plot

		self.grid_circle.update_grid_plot(grid.display_grid.grid)	# plot grid

		self.c44 = self.css2Dtoc44(grid.grid)		# update C44 values
		self.line_c44.set_ydata(self.c44[:, 1])		# plot C44
		self.line_c44.figure.canvas.draw()			# update plot

def suptitle():
	"""
	Returns suptitles for figures.
	"""

	return (
		r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e,$'
		% (parameters['N'], parameters['density'], parameters['vzero'],
		parameters['dr'])
		+ r'$\bar{N}=%.2e, \Delta t=%.2e, nD_0 \Delta t=%.2e$'
		% (Nmean, dt*parameters['period_dump']*parameters['time_step'],
		nD0*dt*parameters['period_dump']*parameters['time_step'])
		+ '\n' + r'$L=%.2e, x_0=%.2e, y_0=%.2e,$' % (box_size, *centre)
		+ r'$S_{init}=%.2e$' % init_frame
		+ r'$, S_{max}=%.2e, N_{cases}=%.2e$' % (int_max, Ncases))

def plot():
	"""
	Plots collective mean square displacements, strain correlations and
	projection of strain correlations on cos(4\\theta).

	Returns
	-------
	sc : active_particles.analysis.ctt.StrainCorrelationsCMSD
		Strain correlations object.
	"""

	sc = StrainCorrelationsCMSD(wave_vectors,
		k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm)
	sc.plot(box_size, r_max_css,
		parameters['box_size']/np.sqrt(parameters['N']),
		r_min, r_max, y_min, y_max,
		points_x_c44=points_x_c44, points_theta_c44=points_theta_c44,
		y_min_c44=y_min_c44, y_max_c44=y_max_c44,
		r_min_c44=r_min_c44, r_max_c44=r_max_c44,
		r_cut=r_cut_fourier, smooth=smooth)

	# COLECTIVE MEAN SQUARE DISPLACEMENTS FIGURE

	sc.fig_cmsd.set_size_inches(16, 16)
	sc.fig_cmsd.subplots_adjust(wspace=0.3)
	sc.fig_cmsd.subplots_adjust(hspace=0.3)

	sc.fig_cmsd.suptitle(suptitle())

	sc.ax_cross.set_xlabel(r'$\lambda/a = 2\pi/ka$')
	sc.ax_cross.set_ylabel(sc.cross_line.get_label()
		+ r'$\times \tilde{\mathcal{G}}(\vec{k}, r_{cut})$')

	sc.ax_dot.set_xlabel(r'$\lambda/a = 2\pi/ka$')
	sc.ax_dot.set_ylabel(sc.dot_line.get_label()
		+ r'$\times \tilde{\mathcal{G}}(\vec{k}, r_{cut})$')

	sc.ax_super.set_xlabel(r'$\lambda/a = 2\pi/ka$')
	sc.ax_super.set_ylabel(
		r'$S(k) \times \tilde{\mathcal{G}}(\vec{k}, r_{cut})$')
	sc.ax_super.add_artist(sc.ax_super.legend())

	sc.ax_cross.set_title(
		r'$\tilde{\mathcal{G}}(\vec{k}, r_{cut}) \equiv$'
		+ r'$\exp(-\frac{1}{2} r_{cut}^2 \vec{k}^2)$')

	if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
		image_name, = naming_Ctt.image().filename(**attributes)
		sc.fig_cmsd.savefig(joinpath(data_dir, image_name))

	if get_env('FITTING_LINE', default=False, vartype=bool):	# FITTING_LINE mode
		sc.fl_cross = FittingLine(sc.ax_cross, slope0, slope_min, slope_max,
			x_fit='(\lambda/a)')
		sc.fl_dot = FittingLine(sc.ax_dot, slope0, slope_min, slope_max,
			x_fit='(\lambda/a)')
		sc.fl_super = FittingLine(sc.ax_super, slope0, slope_min, slope_max,
			x_fit='(\lambda/a)')

	# STRAIN CORRELATIONS AND SUPERIMPOSED ORIGINAL CMSD FIGURE

	sc.fig_sc.set_size_inches(16, 16)
	sc.fig_sc.subplots_adjust(wspace=0.5)
	sc.fig_sc.subplots_adjust(hspace=0.3)

	sc.fig_sc.suptitle(suptitle())

	sc.ax_sc.set_title(
		'2D ' + r'$%s(\vec{r})\ast\mathcal{G}(\vec{r}, \sigma)$' % sc.cor_name
		+ '\n' + r'$\mathcal{G}(\vec{r}, \sigma) \equiv$'
		+ r'$\frac{1}{2\pi\sigma^2}\exp(-\vec{r}^2/2\sigma^2)$')
	sc.ax_sc.set_xlabel(r'$x$')
	sc.ax_sc.set_ylabel(r'$y$')
	sc.colormap.set_label(r'$%s$' % sc.cor_name, labelpad=20, rotation=270)

	sc.ax_kFFTugrid.set_xlabel(sc.ax_super.get_xlabel())
	sc.ax_kFFTugrid.set_ylabel(sc.ax_super.get_ylabel())

	sc.ax_kFFTugrid.legend()

	# GRID CIRCLE FIGURE

	sc.grid_circle.fig.suptitle(suptitle())

	sc.grid_circle.fig.set_size_inches(16, 16)
	sc.grid_circle.fig.subplots_adjust(wspace=0.4)
	sc.grid_circle.fig.subplots_adjust(hspace=0.3)

	sc.grid_circle.ax_grid.set_xlabel(r'$x$')
	sc.grid_circle.ax_grid.set_ylabel(r'$y$')
	sc.grid_circle.colormap.set_label(r'$%s$' % sc.cor_name,
		labelpad=20, rotation=270)

	sc.grid_circle.ax_plot.set_xlabel(r'$\theta$')
	sc.grid_circle.ax_plot.set_ylabel(r'$%s(r, \theta)$' % sc.cor_name)

	# C44 FIGURE

	sc.fig_c44.suptitle(suptitle())

	sc.fig_c44.set_size_inches(16, 16)	# figure size

	sc.ax_c44.set_xlabel(r'$r/a$' + ' ' + r'$(a = L/\sqrt{N})$')
	sc.ax_c44.set_ylabel(r'$C_4^4(r) = \frac{1}{\pi}\int_0^{2\pi}d\theta$'
		+ ' ' + r'$%s(r, \theta)$' % sc.cor_name
		+ ' ' + r'$\cos4\theta$')

	if get_env('FITTING_LINE', default=False, vartype=bool):	# FITTING_LINE mode
		sc.fl_c44 = FittingLine(sc.ax_c44, slope0_c44, slope_min_c44,
			slope_max_c44, x_fit='(r/a)', y_fit='C_4^4(r)')

	# RETURN

	return sc

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    wrap_file_name = get_env('WRAPPED_FILE',
        default=joinpath(data_dir, naming.wrapped_trajectory_file))	# wrapped trajectory file (.gsd)

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

    Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame
    Nframes = Nentries - init_frame									# number of frames available for the calculation

    dt = Nframes + dt if dt <= 0 else dt	# length of the interval of time for which displacements are calculated

    times = linframes(init_frame, Nentries - dt, int_max)	# frames at which shear strain is calculated

	# NAMING

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

        unwrap_file_name = get_env('UNWRAPPED_FILE',
			default=joinpath(data_dir, naming.unwrapped_trajectory_file))	# unwrapped trajectory file (.dat)

        # DISPLACEMENT AND DENSITY CORRELATIONS

        with open(wrap_file_name, 'rb') as wrap_file,\
			open(unwrap_file_name, 'rb') as unwrap_file:	# opens wrapped and unwrapped trajectory files

            w_traj = Gsd(wrap_file, prep_frames=prep_frames)	# wrapped trajectory object
            u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object
            Ugrid = list(map(
                lambda time: displacement_grid(parameters['box_size'],
                box_size, centre, Ncases, time, dt, w_traj, u_traj),
                times))											# lists of displacement variables

        wave_vectors = wave_vectors_2D(Ncases, Ncases, d=box_size/Ncases)	# wave vectors grid
        wave_vectors_norm = np.sqrt(np.sum(wave_vectors**2, axis=-1))		# wave vectors norm grid

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

        with open(wrap_file_name, 'rb') as wrap_file:
            w_traj = Gsd(wrap_file, prep_frames=prep_frames)		# wrapped trajectory object
            Nmean = np.mean(count_particles(
				w_traj, *times, box_size=box_size, centre=centre))	# average number of particles in box
        nD0 = nD0_active(
			Nmean, parameters['vzero'], parameters['dr'], box_size)	# product of particle density and active diffusion constant

    if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

		# PLOT

        r_cut_fourier = get_env('R_CUT_FOURIER', default=_r_cut_fourier,
			vartype=float)	# initial wave length Gaussian cut-off radius

        smooth = get_env('SMOOTH', default=0, vartype=float)	# C44 Gaussian smoothing length scale

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

        points_x_c44 = get_env('POINTS_X_C44', default=_points_x_c44,
			vartype=int)							# number of radii at which to evaluate integrated strain correlation
        points_theta_c44 = get_env('POINTS_THETA_C44',
			default=_points_theta_c44, vartype=int)	# number of angles to evaluate integrated strain correlation

        y_min_c44 = get_env('Y_MIN_C44', default=_y_min_c44, vartype=float)	# minimum plot value for C44
        y_max_c44 = get_env('Y_MAX_C44', default=_y_max_c44, vartype=float)	# maximum plot value for C44

        r_min_c44 = get_env('R_MIN_C44', default=_r_min_c44, vartype=float)	# minimum radius in average particle separation for C44 calculation
        r_max_c44 = get_env('R_MAX_C44', default=_r_max_c44, vartype=float)	# maximum radius in average particle separation for C44 calculation

        slope0_c44 = get_env('SLOPE_C44', default=_slope0_c44,
			vartype=float)	# initial slope for fitting line
        slope_min_c44 = get_env('SLOPE_MIN_C44', default=_slope_min_c44,
			vartype=float)	# minimum slope for fitting line
        slope_max_c44 = get_env('SLOPE_MAX_C44', default=_slope_max_c44,
			vartype=float)	# maximum slope for fitting line

        plot = plot()

        if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
            plt.show()
