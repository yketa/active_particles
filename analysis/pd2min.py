"""
Module pd2min plots color histogram of the probability of nonaffine
squared displacements and most probable nonaffine squared displacement
as functions of the lag time.

Environment modes
-----------------
SUPTITLE : bool
    Display suptitle on figures.
    DEFAULT: True

Environment parameters
----------------------
DATA_DIRECTORY : string
    Data directory.
    DEFAULT: curent working directory
PARAMETERS_FILE : string
	Simulation parameters file.
	DEFAULT: DATA_DIRECTORY/active_particles.naming.parameters_file
WRAPPED_FILE : string
    Wrapped trajectory file. (.gsd)
    DEFAULT: DATA_DIRECTORY/active_particles.naming.wrapped_trajectory_file
INITIAL_FRAME : int
	Reference frame in D2min calculations.
	NOTE: INITIAL_FRAME < 0 will be interpreted as the reference frame
	      being the middle frame of the simulation.
	DEFAULT: -1
DT_MIN : int
	Minimum lag time.
	DEFAULT: 1
DT_MAX : int
	Maximum lag time.
	NOTE: DT_MAX < 0 will be interpreted as the maximum lag time being
	      the number of simulation frames starting from INITIAL_FRAME
		  + DT_MAX.
	DEFAULT: -1
INTERVAL_MAXIMUM : int
	Maximum number of lag times to compute D2min.
	DEFAULT: active_particles.analysis.pd2min._int_max
N_BINS : int
	Number of bins for the histogram.
	DEFAULT: active_particles.analysis.pd2min._Nbins
D2MINMIN : float
	Minimum D2min value.
	DEFAULT: active_particles.analysis.pd2min._d2minmin
D2MINMAX : float
	Maximum D2min value.
	DEFAULT: active_particles.analysis.pd2min._d2minmax
PD2MINMIN : float
	Minimum D2min probability.
	DEFAULT: active_particles.analysis.pd2min._pd2minmin
PD2MINMAX : float
	Maximum D2min probability.
	DEFAULT: active_particles.analysis.pd2min._pd2minmax
CONTOURS : int
	Number of contour lines.
	DEFAULT: active_particles.analysis.pd2min._contours
COLORMAP : string
	Plot colormap.
	DEFAULT: active_particles.analysis.pd2min._colormap
FONT_SIZE : float
	Font size.
	DEFAULT: active_particles.analysis.pd2min._font_size
PAD : float
	Separation between label and colormap.
	DEFAULT: active_particles.analysis.pd2min._colormap_label_pad
"""

import active_particles.naming as naming

from active_particles.init import get_env
from active_particles.dat import Gsd
from active_particles.maths import Histogram

from os import getcwd
from os import environ as envvar
if __name__ == '__main__': envvar['SHOW'] = 'True'
from os.path import join as joinpath

from math import ceil

import pickle

import numpy as np

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

from collections import OrderedDict

# DEFAULT VARIABLES

_int_max = 50	# default maximum number of lag times to compute D2min

_Nbins = 100		# default number of bins for the histogram
_d2minmin = 1e-10	# default minimum D2min value
_d2minmax = 1		# default maximum D2min value

_pd2minmin = 1e-4	# default minimum D2min probability
_pd2minmax = 1e-1	# default maximum D2min probability
_contours = 20      # default contour level value

_font_size = 15             # default font size for the plot
_colormap = 'inferno'       # default plot colormap
_colormap_label_pad = 20    # separation between label and colormap

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLES DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    wrap_file_name = get_env('WRAPPED_FILE',
        default=joinpath(data_dir, naming.wrapped_trajectory_file))	    # wrapped trajectory file (.gsd)

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)  # reference frame in D2min calculations

    dt_min = get_env('DT_MIN', default=1, vartype=int)	# minimum lag time
    dt_max = get_env('DT_MAX', default=-1, vartype=int)	# maximum lag time

    int_max = get_env('INTERVAL_MAXIMUM', default=_int_max, vartype=int)	# maximum number of lag times to compute D2min

    parameters_file = get_env('PARAMETERS_FILE',
        default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    prep_frames = ceil(parameters['prep_steps']/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

    Nentries = parameters['N_steps']//parameters['period_dump']	# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame
    dt_max = Nentries - init_frame + dt_max if dt_max < 0 else dt_max

    lag_times = np.array(list(OrderedDict.fromkeys(map(
        int,
		np.exp(np.linspace(np.log(dt_min), np.log(dt_max), int_max))))))	# lag times logarithmically spaced for the calculation
    pdtsdr = (parameters['period_dump']*parameters['time_step']
        *parameters['dr'])													# number of rotations corresponding to the distance between frames

    # PLOT PARAMETERS

    Nbins = get_env('N_BINS', default=_Nbins, vartype=int)				# number of bins for the histogram
    d2minmin = get_env('D2MINMIN', default=_d2minmin, vartype=float)	# minimum D2min value
    d2minmax = get_env('D2MINMAX', default=_d2minmax, vartype=float)	# maximum D2min value

    pd2minmin = get_env('PD2MINMIN', default=_pd2minmin, vartype=float)	# minimum D2min probability
    pd2minmax = get_env('PD2MINMAX', default=_pd2minmax, vartype=float)	# maximum D2min probability

    contours = get_env('CONTOURS', default=_contours, vartype=int)  # number of contour lines

	# FIGURE PARAMETERS

    colormap = get_env('COLORMAP', default=_colormap)   # plot colormap

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=float) # font size
    mpl.rcParams.update({'font.size': font_size})

    pad = get_env('PAD', default=_colormap_label_pad, vartype=float)    # separation between label and colormap

	# LEGEND SUPTITLE

    display_suptitle = get_env('SUPTITLE', default=True, vartype=bool)  # display suptitle
    suptitle = (
        r'$N=%.1e, \phi=%1.2f,$' % (parameters['N'], parameters['density'])
        + r'$\tilde{v} = %.2e,$' % parameters['vzero']
        + r'$\tilde{\nu}_r = %.2e$' % parameters['dr']
        + '\n' + r'$S_{init} = %.1e$' % init_frame)

	# CALCULATION

    hist = Histogram(Nbins, d2minmin, d2minmax, log=True)	# histogram generator
    bins = np.log10(hist.bins)

    histogram3D = []										# D2min histogram
    d2minpmax = []											# most probable D2min
    with open(wrap_file_name, 'rb') as wrap_file:
        w_traj = Gsd(wrap_file, prep_frames=prep_frames)	# wrapped trajectory object

        for dt in lag_times:
            hist.add_values(
                w_traj.d2min(init_frame, init_frame + dt),
                replace=True)

            drdt_value = np.full(Nbins, fill_value=np.log10(pdtsdr*dt))
            histogram = hist.get_histogram()
            histogram[histogram < pd2minmin] = pd2minmin

            histogram3D += np.transpose(
                [drdt_value, bins, np.log10(histogram)]).tolist()
            d2minpmax += [bins[np.argmax(histogram)]]

	# PLOT

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$\log(\tilde{\nu}_r\Delta t)$')
    ax.set_ylabel(r'$\log(D^2_{min})$')
    if display_suptitle: fig.suptitle(suptitle)

    cmap = plt.get_cmap(colormap)
    norm = colors.Normalize(
		vmin=np.log10(pd2minmin), vmax=np.log10(pd2minmax))
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
        orientation='vertical')
    cb.set_label(r'$\log P(D^2_{min})$', labelpad=pad, rotation=270)

    ax.tricontourf(*np.transpose(histogram3D),
		contours, cmap=cmap, norm=norm)				# D2min histogram
    ax.plot(np.log10(pdtsdr*lag_times), d2minpmax,
        linestyle='--', color='red', linewidth=4)	# most probable D2min line

    plt.show()
