"""
Module msd calculates or plots mean square displacements or distribution of
square displacements.

Files are saved according to the active_particles.naming.Msd standard.

A brief description of the algorithm can be found at:
https://yketa.github.io/UBC_2018_Wiki/#Mean%20square%20displacement

Environment modes
-----------------
DISTRIBUTION : bool
	Analyse distributions of square displacements rather than mean square
	displacements.
	DEFAULT: False
COMPUTE : bool
	Compute mean square displacements.
	DEFAULT: False
PLOT : bool
	Plot saved mean square displacements.
	DEFAULT: False
SHOW [COMPUTE or PLOT mode] : bool
	Show graphs.
	DEFAULT: False
SAVE [COMPUTE or PLOT mode] : bool
	Save graphs.
	DEFAULT: False
FITTING_LINE [COMPUTE or PLOT mode] : bool
	Display fitting line on plot.
	NOTE: see active_particles.plot.mpl_tools.FittingLine
	DEFAULT: False
DIVIDE_BY_DT [COMPUTE or SHOW mode] : bool
	Divide square displacements by lag time.
	DEFAULT: True
SUPTITLE [COMPUTE or SHOW mode] : bool
	Display suptitle on figure.
	DEFAULT: True

Environment parameters
----------------------
DATA_DIRECTORY : string
	Data directory.
	DEFAULT: current working directory
PARAMETERS_FILE : string
	Simulation parameters file.
	DEFAULT: DATA_DIRECTORY/active_particles.naming.parameters_file
UNWRAPPED_FILE : string
	Unwrapped trajectory file. (.dat)
	NOTE: .dat files defined with active_particles.dat
	DEFAULT: DATA_DIRECTORY/active_particles.naming.unwrapped_trajectory_file
INITIAL_FRAME : int
	Frame to consider as initial.
	NOTE: INITIAL_FRAME < 0 will be interpreted as the initial frame being
	      the middle frame of the simulation.
	DEFAULT: -1
INTERVAL_MAXIMUM : int
	Maximum number of intervals of same length dt considered in mean square
	displacement calculations.
	DEFAULT: 1
INTERVAL_PERIOD : int
	Mean square displacement will be calculated for each INTERVAL_PERIOD dumps
	period of time.
	DEFAULT: 1
FONT_SIZE : int
	Font size for the plot.
	DEFAULT: active_particles.plot.pphiloc._font_size
SQ_DISP_MIN [DISTRIBUTION and PLOT or SHOW mode] : float
	Minimum included value of square displacement for histogram bins.
	DEFAULT: active_particles.analysis.msd._sq_disp_min
SQ_DISP_MAX [DISTRIBUTION and PLOT or SHOW mode] : float
	Maximum excluded value of square displacement for histogram bins.
	DEFAULT: active_particles.analysis.msd._sq_disp_max
NBINS [DISTRIBUTION and PLOT or SHOW mode] : int
	Number of histogram bins.
	DEFAULT: active_particles.analysis.msd._Nbins
SLOPE [FITTING_LINE and PLOT or SHOW mode] : float
	Initial slope of fitting line.
	DEFAULT: active_particles.analysis.msd._slope0
SLOPE_MIN [FITTING_LINE and PLOT or SHOW mode] : float
	Minimum slope of fitting line.
	DEFAULT: active_particles.analysis.msd._slope_min
SLOPE_MAX [FITTING_LINE and PLOT or SHOW mode] : float
	Maximum slope of fitting line.
	DEFAULT: active_particles.analysis.msd._slope_max
PSQDISP_MIN : float
	Minimum square displacement probability.
	DEFAULT: active_particles.analysis.msd._psqdispmin
PSQDISP_MAX : float
	Maximum square displacement probability.
	DEFAULT: active_particles.analysis.msd._psqdispmax
CONTOURS [DISTRIBUTION and PLOT or SHOW mode] : int
	Contour level value.
	DEFAULT: active_particles.analysis.msd._contours
COLORMAP [DISTRIBUTION and PLOT or SHOW mode] : string
	Plot colormap.
	DEFAULT: active_particles.analysis.msd._colormap
PAD [DISTRIBUTION and PLOT or SHOW mode] : float
	Separation between label and colormap.
	DEFAULT: active_particles.analysis.msd._colormap_label_pad

Output
------
[COMPUTE mode]
> Prints execution time.
[COMPUTE and not(DISTRIBUTION) mode]
> Saves mean square displacements, lag times and corresponding standard errors
according to the active_particles.naming.Msd(distribution=False) standard in
DATA_DIRECTORY.
[COMPUTE and DISTRIBUTION mode]
> Saves lag times list and corresponding square displacement lists according to
the active_particles.naming.Msd(distribution=True) standard.
[SHOW or PLOT and not(DISTRIBUTION) mode]
> Plots mean square displacements.
[SHOW or PLOT and DISTRIBUTION mode]
> Plots distribution of square displacements.
[SAVE mode]
> Saves figure in DATA_DIRECTORY.
[FITTING LINE mode]
> Adds fitting line to figure.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output
from active_particles.dat import Dat
from active_particles.maths import wo_mean, mean_sterr, Histogram

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

import numpy as np

import pickle

from datetime import datetime

from collections import OrderedDict

import matplotlib as mpl
if not(get_env('SHOW', default=False, vartype=bool)):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

from active_particles.plot.mpl_tools import FittingLine

# DEFAULT VARIABLES

_sq_disp_min = 1e-5	# default minimum included value of square displacement for histogram bins
_sq_disp_max = 1e5	# default maximum excluded value of square displacement for histogram bins
_Nbins = 100		# default number of histogram bins

_psqdispmin = 1e-4  # default minimum square displacement probability
_psqdispmax = 1e-1  # default maximum square displacement probability
_contours = 20      # default contour level value

_slope0 = 1     # default initial slope of fitting line
_slope_min = 0  # default minimum slope of fitting line
_slope_max = 3  # default maximum slope of fitting line

_font_size = 15				# default font size for the plot
_colormap = 'inferno'		# default plot colormap
_colormap_label_pad = 20	# default separation between label and colormap

# FUNCTIONS AND CLASSES

def square_displacement(u_traj, frame, dt):
    """
    Returns square displacement without mean drift between frames frame and
    frame + dt of all particles in uwrapped trajectory file.

    Parameters
    ----------
    u_traj : active_particles.dat.Dat
		Unwrapped trajectory object.
    frame : int
        Initial frame.
    dt : int
        Lag time.

    Returns
    -------
    sq_disp : Numpy array
        Array of all square displacements without mean drift.
    """

    displacements = u_traj.displacement(frame, frame + dt)  # displacements between frame and frame + dt
    return  np.sum(wo_mean(displacements)**2, axis=-1)

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of intervals of same length dt considered in mean square displacement calculations
    int_period = get_env('INTERVAL_PERIOD', default=1, vartype=int) # mean square displacement will be calculated for each int_period dumps period of time

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame

    distribution = get_env('DISTRIBUTION', default=False, vartype=bool)	# DISTRIBUTION mode

    divide_by_dt = get_env('DIVIDE_BY_DT', default=True, vartype=bool)	# DIVIDE_BY_DT mode

    # NAMING

    attributes = {'density': parameters['density'],
        'vzero': parameters['vzero'], 'dr': parameters['dr'],
        'N': parameters['N'], 'init_frame': init_frame, 'int_max': int_max,
        'int_period': int_period}                       # attributes displayed in filenames
    naming_msd = naming.Msd(distribution=distribution)	# mean square displacement naming object
    msd_filename, = naming_msd.filename(**attributes)   # mean square displacement filename

    # STANDARD OUTPUT

    if 'SLURM_JOB_ID' in envvar:	# script executed from Slurm job scheduler
        slurm_output(joinpath(data_dir, 'out'), naming_msd, attributes)

    # MODE SELECTION

    if get_env('COMPUTE', default=False, vartype=bool):	# COMPUTE mode

        startTime = datetime.now()

		# VARIABLE DEFINITIONS

        unwrap_file_name = get_env('UNWRAPPED_FILE',
			default=joinpath(data_dir, naming.unwrapped_trajectory_file))	# unwrapped trajectory file (.dat)

        Nframes = Nentries - init_frame # number of frames available for the calculation
        Ntimes = Nframes//int_period    # number of time intervals considered in the calculation

        lag_times = list(OrderedDict.fromkeys(map(
            int,
            np.exp(np.linspace(np.log(1), np.log(Nframes - 1), Ntimes))
            ))) # lag times logarithmically spaced for the calculation

        # CALCULATION

        with open(unwrap_file_name, 'rb') as unwrap_file,\
            open(joinpath(data_dir, msd_filename),
			'wb' if distribution else 'w') as msd_file:					# opens unwrapped trajectory file and square displacement output file
            if not(distribution): msd_file.write('time, MSD, sterr\n')	# output file header

            u_traj = Dat(unwrap_file, parameters['N'])  # unwrapped trajectory object

            sq_disps = []			# list of square displacements
            for dt in lag_times:    # for each lag time
                lag_time = dt*parameters['period_dump']*parameters['time_step']

                frames = list(OrderedDict.fromkeys(
                    init_frame + np.linspace(0, Nframes - dt - 1,
                    min(int_max, Nframes - dt), dtype=int)
                    ))                              # initial frames for mean square displacement at lag time dt
                sq_disp = list(map(
                    lambda frame: square_displacement(u_traj, frame, dt),
                    frames
                    ))                              # square displacements for lag time dt

                if not(distribution):					# not(DISTRIBUTION) mode
                    msd, sterr = mean_sterr(sq_disp)	# mean square displacement and corresponding standard error
                    msd_file.write('%e,%e,%e\n' % (lag_time, msd, sterr))

                else:	# DISTRIBUTION mode
                    sq_disps += [sq_disp]

            if distribution:
                pickle.dump([lag_times, sq_disps], msd_file)

        # EXECUTION TIME

        print("Execution time: %s" % (datetime.now() - startTime))

    if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

        # DATA

        if distribution:							# DISTRIBUTION mode
            with open(joinpath(data_dir, msd_filename), 'rb') as msd_file:
                lag_times, sq_disp_list = pickle.load(msd_file)
        else:										# not(DISTRIBUTION) mode
	        dt, msd, sterr = np.transpose(np.genfromtxt(
	            fname=joinpath(data_dir, msd_filename),
	            delimiter=',', skip_header=True))   # lag times, mean square displacements and corresponding standard errors

		# PLOT PARAMETERS

        font_size = get_env('FONT_SIZE', default=_font_size, vartype=float)
        mpl.rcParams.update({'font.size': font_size})	# font size for the plot

		# CALCULATION

        if distribution:	# DISTRIBUTION mode

	        sq_disp_min = get_env('SQ_DISP_MIN', default=_sq_disp_min,
				vartype=float)										# minimum included value of square displacement for histogram bins
	        sq_disp_max = get_env('SQ_DISP_MAX', default=_sq_disp_max,
				vartype=float)										# maximum excluded value of square displacement for histogram bins
	        Nbins = get_env('NBINS', default=_Nbins, vartype=int)	# number of histogram bins

	        hist = Histogram(Nbins, sq_disp_min, sq_disp_max, log=True)	# histogram maker
	        bins = np.log10(hist.bins)									# histogram bins
	        histogram3D = []											# 3D histogram

	        for lag_time, sq_disp in zip(lag_times, sq_disp_list):
	            hist.add_values(
					np.array(sq_disp)/lag_time if divide_by_dt else sq_disp,
					replace=True)
	            histogram = np.log10(hist.get_histogram())
	            var_value = np.full(Nbins, fill_value=np.log10(lag_time))
	            histogram3D += np.transpose(
                    [var_value, bins, histogram]).tolist()

	        histogram3D = np.transpose(histogram3D)

        # PLOT

        fig, ax = plt.subplots()

        if get_env('SUPTITLE', default=True, vartype=bool): fig.suptitle(
            r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
    		% (parameters['N'], parameters['density'], parameters['vzero'],
    		parameters['dr']) + '\n' +
            r'$S_{init}=%.2e, S_{max}=%.2e, S_{period}=%.2e$'
            % (init_frame, int_max, int_period))

        if distribution:	# DISTRIBUTION mode

            ax.set_xlabel(r'$\log \Delta t$')
            if divide_by_dt: ax.set_ylabel(r'$\log |\Delta r|^2/\Delta t$')
            else: ax.set_ylabel(r'$\log |\Delta r|^2$')

            pad = get_env('PAD', default=_colormap_label_pad, vartype=float)	# separation between label and colormap

            colormap = get_env('COLORMAP', default=_colormap)   # histogram colormap

            psqdispmin = np.log10(
		        get_env('PSQDISP_MIN', default=_psqdispmin, vartype=float)) 	# minimum square displacement probability
            histogram3D[-1, :][histogram3D[-1, :] < psqdispmin] = psqdispmin	# setting histogram minimum value to psqdispmin
            psqdispmax = np.log10(
		        get_env('PSQDISP_MAX', default=_psqdispmax, vartype=float)) 	# maximum square displacement probability
            histogram3D[-1, :][histogram3D[-1, :] > psqdispmax] = psqdispmax	# setting histogram maximum value to psqdispmin

            contours = get_env('CONTOURS', default=_contours, vartype=int)  # contour level value

            cmap = plt.get_cmap(colormap)
            norm = colors.Normalize(vmin=psqdispmin, vmax=psqdispmax)
            scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
		        orientation='vertical')
            cb.set_label(r'$\log P(|\Delta r(\Delta t)|^2)$',
				labelpad=pad, rotation=270)

            ax.tricontourf(*histogram3D, contours, cmap=cmap, norm=norm) # square displacement histogram

        else:	# not(DISTRIBUTION) mode

            ax.set_xlabel(r'$\Delta t$')
            ax.set_xscale('log')
            if divide_by_dt:
	            ax.set_ylabel(r'$<|\Delta r(\Delta t)|^2/Delta t>$')
            else: ax.set_ylabel(r'$<|\Delta r(\Delta t)|^2>$')
            ax.set_yscale('log')

            ax.errorbar(dt, msd/dt if divide_by_dt else msd,
				yerr=sterr/dt if divide_by_dt else sterr)

        if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
            image_name, = naming_msd.image().filename(**attributes)
            fig.savefig(joinpath(data_dir, image_name))

        if get_env('FITTING_LINE', default=False, vartype=bool):    # FITTING LINE mode

            slope0 = get_env('SLOPE', default=_slope0, vartype=float)           # initial slope of fitting line
            slope_min = get_env('SLOPE_MIN', default=_slope_min, vartype=float) # minimum slope of fitting line
            slope_max = get_env('SLOPE_MAX', default=_slope_max, vartype=float) # maximum slope of fitting line

            fitting_line = FittingLine(ax, slope0, slope_min, slope_max)    # interactive fitting line

        if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
            plt.show()
