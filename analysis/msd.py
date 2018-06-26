"""
Module msd calculates or plots mean square displacements.

Files are saved according to the active_particles.naming.Msd standard.

A brief description of the algorithm can be found at:
https://yketa.github.io/UBC_2018_Wiki/#Mean%20square%20displacement

Environment modes
-----------------
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

Output
------
[COMPUTE MODE]
> Prints execution time.
> Saves mean square displacements, lag times and corresponding standard errors
according to the active_particles.naming.Msd standard in DATA_DIRECTORY.
[SHOW or PLOT mode]
> Plots mean square displacements.
[SAVE mode]
> Saves mean square displacements figure in DATA_DIRECTORY.
[FITTING LINE mode]
> Adds fitting line to figure.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output
from active_particles.dat import Dat
from active_particles.maths import wo_mean, mean_sterr

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

from active_particles.plot.mpl_tools import FittingLine

# DEFAULT VARIABLES

_slope0 = 1     # default initial slope of fitting line
_slope_min = 0  # default minimum slope of fitting line
_slope_max = 3  # default maximum slope of fitting line

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

    # NAMING

    attributes = {'density': parameters['density'],
        'vzero': parameters['vzero'], 'dr': parameters['dr'],
        'N': parameters['N'], 'init_frame': init_frame, 'int_max': int_max,
        'int_period': int_period}                       # attributes displayed in filenames
    naming_msd = naming.Msd()                           # mean square displacement naming object
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

        # MEAN SQUARE DISPLACEMENT

        with open(unwrap_file_name, 'rb') as unwrap_file,\
            open(joinpath(data_dir, msd_filename), 'w') as msd_file:   # opens unwrapped trajectory file and mean square displacement output file
            msd_file.write('time, MSD, sterr\n')                    # output file header

            u_traj = Dat(unwrap_file, parameters['N'])  # unwrapped trajectory object

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
                msd, sterr = mean_sterr(sq_disp)    # mean square displacement and corresponding standard error

                msd_file.write('%e,%e,%e\n' % (lag_time, msd, sterr))

        # EXECUTION TIME

        print("Execution time: %s" % (datetime.now() - startTime))

    if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

        # DATA

        dt, msd, sterr = np.transpose(np.genfromtxt(
            fname=joinpath(data_dir, msd_filename),
            delimiter=',', skip_header=True))   # lag times, mean square displacements and corresponding standard errors

        # PLOT

        fig, ax = plt.subplots()

        fig.suptitle(
            r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
    		% (parameters['N'], parameters['density'], parameters['vzero'],
    		parameters['dr']) + '\n' +
            r'$S_{init}=%.2e, S_{max}=%.2e, S_{period}=%.2e$'
            % (init_frame, int_max, int_period))

        ax.set_xlabel(r'$\Delta t$')
        ax.set_ylabel(r'$<|\Delta r(\Delta t)|^2>$')

        ax.errorbar(dt, msd, yerr=sterr)

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
