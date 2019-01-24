"""
Module varn calculates local densities of the 2D system and plots histogram
of these local densities.

Files are saved according to the active_particles.naming.varN standard.

Environment modes
-----------------
COMPUTE : bool
	Compute local densities.
	DEFAULT: False
CHECK : bool
	Evaluate difference between the parametrised global packing fraction and
	the measured averaged packing fraction.
	DEFAULT: False
PLOT : bool
	Plots histogram of local densities.
	DEFAULT: False
PLOT_MODE : string
	Histogram type.
	 _______________________________________________________________________
	| Mode   | Histogram                                                    |
	|________|______________________________________________________________|
	| 'mean' | Simple histogram of local densities from all computed times. |
	|________|______________________________________________________________|
	| 'time' | Histogram of local densities as function of time.            |
	|________|______________________________________________________________|
	DEFAULT: mean
SHOW : bool
	Show graphs.
	DEFAULT: False
PEAK [(COMPUTE and SHOW) or PLOT mode] : bool
	Highlight highest peak of the histogram.
	DEFAULT: True
SAVE [(COMPUTE and SHOW) or PLOT mode] : bool
	Save graphs.
	DEFAULT: False
SUPTITLE [(COMPUTE and SHOW) or PLOT mode] : bool
	Display suptitle.
	DEFAULT: True

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
PPHILOCMIN [PLOT or SHOW and 'time' mode] : float
	Minimum local density probability.
    DEFAULT: active_particles.analysis.varn._pphilocmin
PPHILOCMAX [PLOT or SHOW and 'time' mode] : float
	Maximum local density probability.
    DEFAULT: active_particles.analysis.varn._pphilocmax
CONTOURS : int
    Number of contour lines.
    DEFAULT: active_particles.analysis.varn._contours
FONT_SIZE : int
	Plot font size.
	DEFAULT: active_particles.analysis.varn._font_size

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
from active_particles.maths import Histogram

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
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable

# DEFAULT VARIABLES

_init_frame = -1	# default frame to consider as initial
_int_max = 1		# default maximum number of frames on which to calculate densities

_box_size = 10  # default length of the square box in which particles are counted

_Nbins = 10 # default number of bins for the histogram
_phimax = 1 # default maximum local density for histogram

_pphilocmin = 1e-4  # default minimum local density probability
_pphilocmax = 1e-1  # default maximum local density probability
_contours = 20      # default contour level value

_font_size = 10 # default plot font size

# FUNCTIONS AND CLASSES

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

    hist = Histogram(Nbins, 0, phimax)
    hist.add_values(densities)
    return hist.bins, hist.get_histogram()

class Plot:
	"""
	Plot mean histograms of densities.
	"""

	def __init__(self, suptitle=True):
		"""
		Set figure.

		Parameters
		----------
		suptitle : bool
			Display suptitle. (default: True)
		"""

		self.fig, self.ax = plt.subplots()

		if suptitle: self.fig.suptitle(
	        r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
	        % (parameters['N'], parameters['density'], parameters['vzero'],
	        parameters['dr']) + '\n' +
	        r'$S_{init}=%.2e, S_{max}=%.2e, N_{cases}=%.2e, l=%.2e$'
	        % (init_frame, int_max, Ncases, box_size))

		self.ax.set_xlabel(r'$\phi_{loc}$')
		self.ax.set_ylabel(r'$P(\phi_{loc})$')

	def add_hist(self, bins, hist, peak=True):
		"""
		Add histogram of densities.

		Parameters
	    ----------
	    bins : array-like
	        Bins of the histogram.
	    hist : array-like
	        Values of the histogram at bins.
		peak : bool
			Highlight tallest peak of histogram. (default: True)

		Returns
		-------
		line : matplotlib.lines.Line2D
			Plotted histogram line.
		"""

		line, = self.ax.semilogy(bins, hist)

		if peak:
			philocmax, Pphilocmax = bins[np.argmax(hist)], np.max(hist)
			self.ax.axhline(Pphilocmax, 0, 1,
				linestyle='--', color=line.get_color())
			self.ax.axvline(philocmax, 0, 1,
				linestyle='--', color=line.get_color(),
				label=r'$(\phi_{loc}^* = %1.2f, P(\phi_{loc}^*) = %.2e)$'
				% (philocmax, Pphilocmax))

		return line

class PlotTime:
	"""
	Plot histograms of densities as functions of time.
	"""

	def __init__(self, Nbins, phimax,
		pphilocmin=_pphilocmin, pphilocmax=_pphilocmax, contours=_contours,
		colormap=plt.cm.inferno, pad=20, suptitle=True):
		"""
		Set figure and histogram parameters.

		Parameters
		----------
		Nbins : int
			Number of bins for the histogram.
		phimax : float
			Maximum local density for histogram.
		pphilocmin : float
			Minimum local density probability.
			(default: active_particles.analysis.varn._pphilocmin)
		pphilocmax : float
			Maximum local density probability.
			(default: active_particles.analysis.varn._pphilocmax)
		contours : int
			Number of contour lines.
			(default: active_particles.analysis.varn._contours)
		colormap : matplotlib colormap
			Histogram colormap. (default: matplotlib.pyplot.cm.inferno)
		pad : float
			Separation between label and colormap. (default: 20)
		suptitle : bool
			Display suptitle. (default: True)
		"""

		self.Nbins = Nbins
		self.phimax = phimax
		self.pphilocmin = np.log10(pphilocmin)
		self.pphilocmax = np.log10(pphilocmax)
		self.contours = contours

		self.fig, self.ax = plt.subplots()
		self.cmap = colormap
		self.norm = colors.Normalize(
			vmin=self.pphilocmin, vmax=self.pphilocmax)
		self.colorbar = mpl.colorbar.ColorbarBase(
			make_axes_locatable(self.ax).append_axes(
				"right", size="5%", pad=0.05),
			cmap=self.cmap, norm=self.norm, orientation='vertical')

		if suptitle: self.fig.suptitle(
	        r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
	        % (parameters['N'], parameters['density'], parameters['vzero'],
	        parameters['dr']) + '\n' +
	        r'$S_{max}=%.2e, N_{cases}=%.2e, l=%.2e$'
			% (int_max, Ncases, box_size))

		self.ax.set_xlabel(r'$t$')
		self.ax.set_ylabel(r'$\phi_{loc}$')
		self.colorbar.set_label(r'$\log P(\phi_{loc})$',
			labelpad=pad, rotation=270)

	def plot(self, times, densities, peak=True):
		"""
		Plot histogram.

		Parameters
		----------
		times : array-like
			Array of times at which densities have been calculated.
		densities : array-like of array-like
			Array of densities arrays at times times.
		peak : bool
			Highlight tallest peak of histogram. (default: True)
		"""

		self.times = times

		self.histogram3D = []	# local densities histogram
		self.philocmax = []		# most probable local densities at times

		for time, density in zip(self.times, densities):

			time_value = np.full(self.Nbins, fill_value=time)

			bins, hist = histogram(density, self.Nbins, self.phimax)	# histogram of local densities with corresponding bins
			hist = np.log10(hist)

			histogram3D_time = np.transpose([time_value, bins, hist]).tolist()
			self.histogram3D += histogram3D_time
			self.philocmax += [max(histogram3D_time, key=lambda el: el[2])[1]]

		self.histogram3D = np.transpose(self.histogram3D)
		self.histogram3D[2][
			self.histogram3D[2] < self.pphilocmin] = self.pphilocmin	# set minimum histogram value as pphilocmin

		self.ax.tricontourf(*self.histogram3D, self.contours,
			cmap=self.cmap, norm=self.norm)				# local density histogram
		if peak: self.ax.plot(self.times, self.philocmax,
			linestyle='--', color='red', linewidth=4)	# most probable packing fraction line

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    init_frame = get_env('INITIAL_FRAME', default=_init_frame, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=_int_max, vartype=int)	# maximum number of frames on which to calculate densities

    box_size = get_env('BOX_SIZE', default=_box_size, vartype=float)    # length of the square boxes in which particles are counted

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    prep_frames = ceil(parameters['prep_steps']/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

    Nentries = parameters['N_steps']//parameters['period_dump']			# number of time snapshots in unwrapped trajectory file
    Nentries = get_env('FINAL_FRAME', default=Nentries, vartype=int)	# final frame to consider
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame		# initial frame

    frames = list(OrderedDict.fromkeys(map(
		int,
		np.linspace(init_frame, Nentries - 1, int_max)
		)))	# linearly spaced frames at which to calculate the densities

    Ncases = get_env('N_CASES', default=ceil(np.sqrt(parameters['N'])),
		vartype=int)	# number of boxes in each direction to compute the local density

    # NAMING

    attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'init_frame': init_frame, 'int_max': int_max,
        'fin_frame': Nentries, 'Ncases': Ncases, 'box_size': box_size}	# attributes displayed in filenames
    naming_varN = naming.VarN(final_frame='FINAL_FRAME' in envvar)		# varN naming object
    varN_filename, = naming_varN.filename(**attributes) 				# varN file name

    # STANDARD OUTPUT

    if 'SLURM_JOB_ID' in envvar:	# script executed from Slurm job scheduler
        slurm_output(joinpath(data_dir, 'out'), naming_varN, attributes)

    # MODE SELECTION

    if get_env('COMPUTE', default=False, vartype=bool):	# COMPUTE mode

        startTime = datetime.now()

		# VARIABLE DEFINITIONS

        wrap_file_name = get_env('WRAPPED_FILE',
			default=joinpath(data_dir, naming.wrapped_trajectory_file))		# wrapped trajectory file (.gsd)

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

    if get_env('CHECK', default=False, vartype=bool):	# CHECK mode

		# DATA

        with open(joinpath(data_dir, varN_filename), 'rb') as varN_dump_file:
            densities = pickle.load(varN_dump_file)

		# CHECK

        mean_density = np.mean(densities)
        difference = np.abs(mean_density - parameters['density'])
        relative_difference = difference/parameters['density']

        print('Parametrised packing fraction: %f' % parameters['density'])
        print('Measured averaged packing fraction: %f' % mean_density)
        print('Difference: %f' % difference)
        print('Relative difference: %f' % relative_difference)

    if get_env('PLOT', default=False, vartype=bool):	# PLOT mode

		# DATA

        with open(joinpath(data_dir, varN_filename), 'rb') as varN_dump_file:
            densities = pickle.load(varN_dump_file)

    if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

		# PLOT

        Nbins = get_env('N_BINS', default=_Nbins, vartype=int)      # number of bins for the histogram
        phimax = get_env('PHIMAX', default=_phimax, vartype=float)  # maximum local density for histogram

        peak = get_env('PEAK', default=True, vartype=bool)	# highlight highest peak of the histogram

        suptitle = get_env('SUPTITLE', default=True, vartype=bool)	# display suptitle

        font_size = get_env('FONT_SIZE', default=_font_size, vartype=int)	# plot font size
        mpl.rcParams.update({'font.size': font_size})

        mode = get_env('PLOT_MODE', default='mean')	# histogram plot mode

        if mode == 'mean':

	        plot = Plot(suptitle=suptitle)
	        plot.add_hist(*histogram(densities, Nbins, phimax), peak=peak)
	        if peak: plot.ax.legend()	# display legend of highlighted peaks in histograms

        elif mode == 'time':

	        pphilocmin = get_env('PPHILOCMIN',
				default=_pphilocmin, vartype=float)							# minimum local density probability
	        pphilocmax = get_env('PPHILOCMIN',
				default=_pphilocmax, vartype=float)							# maximum local density probability
	        contours = get_env('CONTOURS', default=_contours, vartype=int)	# number of contour lines

	        plot = PlotTime(Nbins, phimax,
				pphilocmin=pphilocmin, pphilocmax=pphilocmax,
				contours=contours, suptitle=suptitle)
	        plot.plot(
				parameters['period_dump']*parameters['time_step']
					*np.array(frames),
				densities, peak=peak)

        else: raise ValueError('Mode %s is not known.' % mode)	# mode is not known

		# SAVING

        if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
	        image_name, = naming_varN.image().filename(**attributes)
	        plot.fig.savefig(joinpath(data_dir, image_name))

		# SHOW

        if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
            plt.show()
