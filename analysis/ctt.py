"""
Module ctt calculates or plots mean squared cross and dot products of wave
vectors and displacement Fourier transform.

Files are saved according to active_particles.naming.Ctt (cross product) and active_particles.naming.Cll (dot product) naming standard.

Environment modes
-----------------
COMPUTE : bool
	Compute displacement Fourier transform and mean squared cross and dot
	products with wave vectors.
	DEFAULT: False
PLOT : bool
	Plot saved cylindrical averages of mean squared cross and dot products of
	wave vectors and displacement Fourier transform as functions of wave vector
	norm.
	DEFAULT: False
SHOW [COMPUTE or PLOT mode] : bool
	Show graphs.
	DEFAULT: False
SAVE [COMPUTE or PLOT mode] : bool
	Save graphs.
	DEFAULT: False
FITTING_LINE [SHOW mode] : bool
	Display adjustable fitting line on graphs.
	DEFAULT : false

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
K_MIN [PLOT or SHOW mode] : float
	Minimum wave vector norm for plots.
	DEFAULT: active_particles.analysis.ctt._k_min
K_MAX [PLOT or SHOW mode] : float
	Maximum wave vector norm for plots.
	DEFAULT: active_particles.analysis.ctt._k_max
Y_MIN [PLOT or SHOW mode] : float
	Minimum y-coordinate for plots.
	DEFAULT: active_particles.analysis.ctt._y_min
Y_MAX [PLOT or SHOW mode] : float
	Maximum y-coordinate for plots.
	DEFAULT: active_particles.analysis.ctt._y_max
SLOPE [FITTING_LINE mode] : float
	Initial slope for fitting line.
	DEFAULT: active_particles.analysis.ctt._slope0
SLOPE_MIN [FITTING_LINE mode] : float
	Minimum slope for fitting line.
	DEFAULT: active_particles.analysis.ctt._slope_min
SLOPE_MAX [FITTING_LINE mode] : float
	Maximum slope for fitting line.
	DEFAULT: active_particles.analysis.ctt._slope_max

Output
------
[COMPUTE MODE]
> Prints execution time.
> Saves wave vectors grid, 2D grid and 1D cylindrical average of mean squared
cross products of wave vectors and displacement Fourier transform according to
active_particles.naming.Ctt standards in DATA_DIRECTORY.
> Saves wave vectors grid, 2D grid and 1D cylindrical average of mean squared
dot products of wave vectors and displacement Fourier transform according to
active_particles.naming.Cll standards in DATA_DIRECTORY.
[SHOW or PLOT mode]
> Plots cylindrical averages as functions of wave vector norm.
[SAVE mode]
> Saves figure in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output
from active_particles.dat import Dat, Gsd
from active_particles.maths import g2Dto1Dgrid, kFFTgrid

from active_particles.analysis.cuu import displacement_grid
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

def plot():
	"""
	Plots mean square norm of dot and cross products of wave vector and
	displacement Fourier transform.

	Returns
	-------fig, ax_cross, ax_dot, ax_super, fl_cross, fl_dot, fl_super
	fig : matplotlib figure
		Figure.
	ax_cross : matplotlib axis
		Mean square cross product axis.
	ax_dot : matplotlib axis
		Mean square dot product axis.
	ax_super : matplotlib axis
		Superimposed mean square cross and dot products axis.
	fl_cross : active_particles.plot.mpl_tools.FittingLine
		Fitting line object for cross axis.
		NOTE: Only in FITTING_LINE mode.
	fl_dot : active_particles.plot.mpl_tools.FittingLine
		Fitting line object for dot axis.
		NOTE: Only in FITTING_LINE mode.
	fl_super : active_particles.plot.mpl_tools.FittingLine
		Fitting line object for superimposition axis.
		NOTE: Only in FITTING_LINE mode.
	"""

	fig = plt.figure()
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
	ax_cross.loglog(k_cross_FFTugrid1D_sqnorm[1:, 0],
		k_cross_FFTugrid1D_sqnorm[1:, 1],
		color='blue',
		label=r'$\left<||\vec{k}\times\tilde{\vec{u}}(\vec{k})||^2\right>$')	# cross product

	ax_cross.set_xlabel(r'$k$')
	ax_cross.set_xlim([k_min, k_max])
	ax_cross.set_ylabel(
		r'$\left<||\vec{k}\times\tilde{\vec{u}}(\vec{k})||^2\right>$')
	ax_cross.set_ylim([y_min, y_max])

	if get_env('FITTING_LINE', default=False, vartype=bool):
		fl_cross = FittingLine(ax_cross, slope0, slope_min, slope_max)	# add fitting line to plot

	# DOT

	ax_dot = plt.subplot(gs[1, 0])
	ax_dot.loglog(k_dot_FFTugrid1D_sqnorm[1:, 0],
		k_dot_FFTugrid1D_sqnorm[1:, 1],
		color='orange',
		label=r'$\left<||\vec{k}\cdot\tilde{\vec{u}}(\vec{k})||^2\right>$')		# dot product

	ax_dot.set_xlabel(r'$k$')
	ax_dot.set_xlim([k_min, k_max])
	ax_dot.set_ylabel(
		r'$\left<||\vec{k}\cdot\tilde{\vec{u}}(\vec{k})||^2\right>$')
	ax_dot.set_ylim([y_min, y_max])

	if get_env('FITTING_LINE', default=False, vartype=bool):
		fl_dot = FittingLine(ax_dot, slope0, slope_min, slope_max)	# add fitting line to plot

	# SUPER

	ax_super = plt.subplot(gs[:, 1])
	ax_super.loglog(k_cross_FFTugrid1D_sqnorm[1:, 0],
		k_cross_FFTugrid1D_sqnorm[1:, 1],
		color='blue',
		label=r'$\left<||\vec{k}\times\tilde{\vec{u}}(\vec{k})||^2\right>$')	# cross product
	ax_super.loglog(k_dot_FFTugrid1D_sqnorm[1:, 0],
		k_dot_FFTugrid1D_sqnorm[1:, 1],
		color='orange',
		label=r'$\left<||\vec{k}\cdot\tilde{\vec{u}}(\vec{k})||^2\right>$')		# dot product

	legend = ax_super.legend()
	ax_super.add_artist(legend)

	ax_super.set_xlabel(r'$k$')
	ax_super.set_xlim([k_min, k_max])
	ax_super.set_ylabel(r'$k^2S(k)$')
	ax_super.set_ylim([y_min, y_max])

	if get_env('FITTING_LINE', default=False, vartype=bool):
		fl_super = FittingLine(ax_super, slope0, slope_min, slope_max)	# add fitting line to plot

	# SAVING

	if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
		image_name, = naming_Ctt.image().filename(**attributes)
		fig.savefig(joinpath(data_dir, image_name))

	# RETURN

	try:
		return fig, ax_cross, ax_dot, ax_super, fl_cross, fl_dot, fl_super
	except NameError: return fig, ax_cross, ax_dot, ax_super

# DEFAULT VARIABLES

_k_min = 1e-3	# default minimum wave vector norm for plots
_k_max = 1		# default maximum wave vector norm for plots

_y_min = 1e-6	# default minimum y-coordinate for plots
_y_max = 1		# default maximum y-coordinate for plots

_slope0 = 2		# default initial slope for fitting line
_slope_min = 0	# default minimum slope for fitting line
_slope_max = 4	# default maximum slope for fitting line

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    dt = get_env('TIME', default=-1, vartype=int)	# lag time for displacement

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of intervals of length dt considered in correlations calculations

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
		vartype=int)        # number of boxes in each direction to compute the shear strain and displacement vorticity grid
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

        wave_vectors, k_cross_dot_FFTUgrid = np.transpose(
			list(map(lambda grid: kFFTgrid(grid, d=dL), Ugrid)),
			(1, 0, 2, 3, 4))								# wave vector grids, and cross and dot products of wave vectors with displacement grids Fourier transform
        k_cross_FFTUgrid = k_cross_dot_FFTUgrid[:, :, :, 0]	# cross product of wave vectors with displacement grids Fourier transform
        k_dot_FFTUgrid = k_cross_dot_FFTUgrid[:, :, :, 1]	# dot product of wave vectors with displacement grids Fourier transform

        wave_vectors = np.mean(wave_vectors, axis=0).real	# grid of wave vectors
        k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm = list(map(
			lambda Grid: np.mean(
			list(map(lambda grid: (np.conj(grid)*grid).real, Grid)),
			axis=0),
			[k_cross_FFTUgrid, k_dot_FFTUgrid]))			# grid of mean square norms of cross and dot products of wave vectors with displacement grids Fourier transform

        wave_vectors_norm = np.sqrt(np.sum(wave_vectors**2, axis=-1))	# grid of wave vectors norm
        k_cross_FFTugrid1D_sqnorm, k_dot_FFTugrid1D_sqnorm = list(map(
			lambda grid2D: g2Dto1Dgrid(grid2D, wave_vectors_norm),
			[k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm]))		# cylindrical averages of mean square norms of cross and dot products of wave vectors with displacement grids Fourier transform

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

        with open(joinpath(data_dir, Ctt_filename), 'wb') as Ctt_dump_file,\
			open(joinpath(data_dir, Cll_filename), 'wb') as Cll_dump_file:
            _, _, k_cross_FFTugrid1D_sqnorm = pickle.load(Ctt_dump_file)
            _, _, k_dot_FFTugrid1D_sqnorm = pickle.load(Cll_dump_file)

    if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

		# PLOT

        k_min = get_env('K_MIN', default=_k_min, vartype=float)	# minimum wave vector norm
        k_max = get_env('K_MAX', default=_k_max, vartype=float)	# maximum wave vector norm

        y_min = get_env('Y_MIN', default=_y_min, vartype=float)	# minimum y-coordinate for plots
        y_max = get_env('Y_MAX', default=_y_max, vartype=float)	# maximum y-coordinate for plots

        slope0 = get_env('SLOPE', default=_slope0, vartype=float)			# initial slope for fitting line
        slope_min = get_env('SLOPE_MIN', default=_slope_min, vartype=float)	# minimum slope for fitting line
        slope_max = get_env('SLOPE_MAX', default=_slope_max, vartype=float)	# maximum slope for fitting line

        plot = plot()

        if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
            plt.show()
