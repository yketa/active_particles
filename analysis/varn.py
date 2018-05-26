"""
Module varn calculates local densities of the 2D system and plots histogram
of these local densities.

Files are saved according to the active_particles.naming.varN standard.

Environment modes
-----------------
COMPUTE : bool
	Compute local densities.
	DEFAULT: False
PLOT : bool
	Plots histogram of local densities.
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
INITIAL_FRAME : int
	Frame to consider as initial.
	NOTE: INITIAL_FRAME < 0 will be interpreted as the initial frame being
	the middle frame of the simulation.
	DEFAULT: -1
INTERVAL_MAXIMUM : int
	Maximum number of frames at which we compute local densities.
	DEFAULT: 1
BOX_SIZE : float
	Length of the square boxes in which particles are counted to compute local
	densities.
	DEFAULT: active_particles.analysis.varn._box_size
N_CASES : int
	Number of boxes in each direction to compute the shear strain and
	displacement vorticity grid.
	DEFAULT: smallest integer value greater than or equal to the square root of
		     the number of particles from the simulation parameters file.
N_BINS [PLOT or SHOW mode] : int
	Number of bins for the histogram of local densities.
	DEFAULT: active_particles.analysis.varn._Nbins
PHIMAX [PLOT or SHOW mode] : int
	Maximum local density for the histogram of local densities.
	DEFAULT: active_particles.analysis.varn._phimax

Output
------
[COMPUTE MODE]
> Prints neigbours grid computation time and execution time.
> Saves computed local densities according to the active_particles.naming.varN
standard in DATA_DIRECTORY.
[SHOW or PLOT mode]
> Plots histogram of local densities.
[SAVE mode]
> Saves local densities histogram figure in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output
from active_particles.dat import Gsd

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

import numpy as np

from math import ceil

import pickle

from collections import OrderedDict

from datetime import datetime

import matplotlib as mpl
if not(get_env('SHOW', default=False, vartype=bool)):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt

def density(w_traj, frame, Ncases, box_size):
    """
    Returns local densities in squares of length box_size around
    Ncases x Ncases nodes, uniformly distributed in the 2D system, at frame
    'frame'.

    Parameters
    ----------
    w_traj : active_particles.dat.Gsd
		Wrapped trajectory object.
    frame : int
        Frame index.
    Ncases : int
        Number of nodes in one direction.
    box_size : float
        Length of the square box in which we calculate the local density.

    Returns
    -------
    density_list : 1D Numpy array
        Array of calculated local densities.
    """

    L = w_traj.box_size()               # system box size
    dL = L/Ncases                       # distance between two consecutive nodes
    max_node_dist = ceil(box_size/dL)   # maximum distance in infinity norm in terms of nodes between particle and containing node

    area_sum = np.zeros((Ncases, Ncases))   # sum of particles' area close to each node (centre of grid box) of the system

    def node_position(node_index):
        """
        Returns node position from node index.

        Parameters
        ----------
        node_index : 2-uple of int
            Node index.

        Returns
        -------
        r : (2,) Numpy array
            Position of node.
        """

        return dL*(1/2 + np.array(node_index)) - L/2

    for position, area in zip(w_traj.position(frame),
        (np.pi/4)*(w_traj.diameter(frame)**2)):

        closest_node_index = np.array((position + L/2)//dL, dtype=int)  # index of closest node

        for dx in range(-max_node_dist, max_node_dist + 1):
            for dy in range(-max_node_dist, max_node_dist + 1):

                node_index = tuple(
                    (closest_node_index + np.array([dx, dy]))%Ncases)

                if (np.abs(position - node_position(node_index))
                    < box_size/2).all():    # particle within box of node
                    area_sum[node_index] += area

    return area_sum/(box_size**2)

def histogram(densities, Nbins, phimax):
    """
    Returns histogram and bin values from densities array.

    Parameters
    ----------
    densities : array-like
        Array of densities.
    Nbins : int
        Number of bins for histogram.
    phimax : float
        Maximum density for hsistogram.
        NOTE: Minimum density is 0.

    Returns
    -------
    bins : Numpy array
        Bins of the histogram.
    hist : Numpy array
        Values of the histogram at bins.
    """

    densities = np.array(densities)

    bins = np.linspace(0, phimax, Nbins, endpoint=False)
    bins_width = phimax/Nbins   # distance between bin values

    hist = np.zeros(Nbins)
    for density in densities.flatten():
        if density < phimax:
            hist[int(density//bins_width)] += 1
        else:
            hist[-1] += 1

    return bins, hist/densities.size

def plot(bins, hist):
    """
    Plots histogram of densities.

    Parameters
    ----------
    bins : array-like
        Bins of the histogram.
    hist : array-like
        Values of the histogram at bins.
    """

    fig, ax = plt.subplots()

    fig.suptitle(
        r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
        % (parameters['N'], parameters['density'], parameters['vzero'],
        parameters['dr']) + '\n' +
        r'$S_{init}=%.2e, S_{max}=%.2e, N_{cases}=%.2e, l=%.2e$'
        % (init_frame, int_max, Ncases, box_size))

    ax.set_xlabel(r'$\phi_{loc}$')
    ax.set_ylabel(r'$P(\phi_{loc})$')

    ax.semilogy(bins, hist)

    # SAVING

    if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
        image_name, = naming_varN.image().filename(**attributes)
        fig.savefig(joinpath(data_dir, image_name))

# DEFAULT VARIABLES

_box_size = 10  # default length of the square box in which particles are counted

_Nbins = 10 # default number of bins for the histogram
_phimax = 1 # default maximum local density for histogram

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of frames on which to calculate densities

    box_size = get_env('BOX_SIZE', default=_box_size, vartype=float)    # length of the square boxes in which particles are counted

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    prep_frames = ceil(parameters['prep_steps']/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

    Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame

    Ncases = get_env('N_CASES', default=ceil(np.sqrt(parameters['N'])),
		vartype=int)	# number of boxes in each direction to compute the shear strain and displacement vorticity grid

    # NAMING

    attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'init_frame': init_frame, 'int_max': int_max,
        'Ncases': Ncases, 'box_size': box_size}         # attributes displayed in filenames
    naming_varN = naming.VarN()                         # varN naming object
    varN_filename, = naming_varN.filename(**attributes) # varN file name

    # STANDARD OUTPUT

    if 'SLURM_JOB_ID' in envvar:	# script executed from Slurm job scheduler
        slurm_output(joinpath(data_dir, 'out'), naming_varN, attributes)

    # MODE SELECTION

    if get_env('COMPUTE', default=False, vartype=bool):	# COMPUTE mode

        startTime = datetime.now()

		# VARIABLE DEFINITIONS

        wrap_file_name = get_env('WRAPPED_FILE',
			default=joinpath(data_dir, naming.wrapped_trajectory_file))		# wrapped trajectory file (.gsd)

        frames = list(OrderedDict.fromkeys(map(
            int,
            np.linspace(init_frame, Nentries - 1, int_max)
            ))) # logarithmically spaced frames at which to calculate the densities

        with open(wrap_file_name, 'rb') as wrap_file:   # opens wrapped trajectory file

            w_traj = Gsd(wrap_file, prep_frames=prep_frames)    # wrapped trajectory object

            densities = list(map(
                lambda frame: density(w_traj, frame, Ncases, box_size),
                frames))    # density lists at frames

        # SAVING

        with open(joinpath(data_dir, varN_filename), 'wb') as varN_dump_file:
            pickle.dump(densities, varN_dump_file)

        # EXECUTION TIME

        print("Execution time: %s" % (datetime.now() - startTime))

    if get_env('PLOT', default=False, vartype=bool):	# PLOT mode

		# DATA

        with open(joinpath(data_dir, varN_filename), 'rb') as varN_dump_file:
            densities = pickle.load(varN_dump_file)

    if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

		# PLOT

        Nbins = get_env('N_BINS', default=_Nbins, vartype=int)      # number of bins for the histogram
        phimax = get_env('PHIMAX', default=_phimax, vartype=float)  # maximum local density for histogram

        bins, hist = histogram(densities, Nbins, phimax)    # histogram and bin values from densities array
        plot(bins, hist)

        if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
            plt.show()
