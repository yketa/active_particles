"""
Module msd_comparison plots mean square displacements for different initial
frames.

Input files in simulation directory must follow the active_paricles.naming.Msd
standard.

Environment modes
-----------------
FITTING_LINE : bool
    Display adjustable fitting line on plot.
    DEFAULT: False
MULTIPLY_WITH_DR : bool
    Plot as function of the lag time multiplied with the rotation diffusion
    rate rather than the sole lag time.
    DEFAULT: True

Environment parameters
----------------------
DATA_DIRECTORY : string
	Data directory.
	DEFAULT: current working directory
PARAMETERS_FILE : string
	Simulation parameters file.
	DEFAULT: DATA_DIRECTORY/active_particles.naming.parameters_file
INTERVAL_MAXIMUM : int
	maximum number of time snapshots taken for the calculation of the mean
    square displacement at each time.
	DEFAULT: 1
INTERVAL_PERIOD : int
    Period of time at which mean square displacement was calculated.
    DEFAULT: 1
FONT_SIZE : int
    Font size.
    DEFAULT: active_particles.plot.msd_comparison._font_size
MARKER_SIZE : int
    Marker size.
    DEFAULT: active_particles.plot.msd_comparison._marker_size
COLORMAP : string
    Plot colormap.
    DEFAULT: active_particles.plot.msd_comparison._colormap
SLOPE [FITTING_LINE mode] : float
    Initial slope for fitting line slider.
    DEFAULT: active_particles.plot.msd_comparison._slope0
SLOPE_MIN [FITTING_LINE mode] : float
    Minimum slope for fitting line slider.
    DEFAULT: active_particles.plot.msd_comparison._slope_min
SLOPE_MAX [FITTING_LINE mode] : float
    Maximum slope for fitting line slider.
    DEFAULT: active_particles.plot.msd_comparison._slope_max
"""

from active_particles import naming

from active_particles.init import get_env

from active_particles.plot.plot import list_colormap
from active_particles.plot.mpl_tools import FittingLine

from os import getcwd
from os import environ as envvar
if __name__ == '__main__': envvar['SHOW'] = 'True'
from os.path import join as joinpath

import pickle

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

# DEFAULT VARIABLES

_slope0 = 0     # default initial slope for fitting line slider
_slope_min = -2 # minimum slope for fitting line slider
_slope_max = 2  # maximum slope for fitting line slider

_font_size = 10     # default font size
_marker_size = 20   # default marker size

_colormap = 'jet'   # default plot colormap

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)   # maximum number of time snapshots taken for the calculation of the mean square displacement at each time
    int_period = get_env('INTERVAL_PERIOD', default=1, vartype=int) # period of time at which mean square displacement was calculated

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    multiply_with_dr = get_env('MULTIPLY_WITH_DR', default=True, vartype=bool)  # multiply lag time with rotational diffusion rate

    # NAMING

    attributes = {'density': parameters['density'],
        'vzero': parameters['vzero'], 'dr': parameters['dr'],
        'N': parameters['N'], 'int_max': int_max, 'int_period': int_period} # attributes displayed in file names
    naming_Msd = naming.Msd()                                               # mean square displacement naming object

    # PLOT PARAMETERS

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=int)       # font size
    marker_size = get_env('MARKER_SIZE', default=_marker_size, vartype=int) # marker size
    mpl.rcParams.update(
        {'font.size': font_size, 'lines.markersize': marker_size})

    colormap = get_env('COLORMAP', default=_colormap)   # plot colormap

    # CALCULATION

    files = naming_Msd.get_files(directory=data_dir, **attributes)  # files corresponding to parameters
    init_list = np.array(list(map(
        lambda file: naming_Msd.get_data(file, 'init_frame'),
        files))).flatten()                                          # list of lag times corresponding to files

    dt, msd, sterr = {}, {}, {} # lag times, mean square displacements, and standard error hash tables with initial frames as keys
    for file, init in zip(files, init_list):
        dt[init], msd[init], sterr[init] = np.transpose(np.genfromtxt(
            fname=joinpath(data_dir, file), delimiter=',', skip_header=True))

    init_list.sort()

    # PLOT

    colors = list_colormap(init_list, colormap=colormap)    # hash table of line colors with initial frames as keys

    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(
        r'$\tilde{\nu}_r \Delta t$' if multiply_with_dr else r'$\Delta t$')
    ax.set_ylabel(r'$<|\Delta r(\Delta t)|^2>/\Delta t$')

    for init in init_list:
        ax.errorbar(
            dt[init]*parameters['dr'] if multiply_with_dr else dt[init],
            msd[init]/dt[init], yerr=sterr[init]/dt[init],
            color=colors[init], label=r'$S_{init} = %d$' % init)
    ax.add_artist(plt.legend())

    if get_env('FITTING_LINE', default=False, vartype=bool):    # FITTING_LINE mode

        slope0 = get_env('SLOPE', default=_slope0, vartype=float)           # initial slope for fitting line slider
        slope_min = get_env('SLOPE_MIN', default=_slope_min, vartype=float) # minimum slope for fitting line slider
        slope_max = get_env('SLOPE_MAX', default=_slope_max, vartype=float) # maximum slope for fitting line slider

        fitting_line = FittingLine(ax, slope0, slope_min, slope_max)
        fitting_line.line.set_zorder(10)

    title = r'$N=%.2e, \phi=%1.2f,$' % (parameters['N'], parameters['density'])
    title += r'$\tilde{v}=%.2e,$' % parameters['vzero']
    title += r'$\tilde{\nu}_r=%.2e$' % parameters['dr']
    title += '\n'
    title += r'$S_{max}=%.2e, S_{period}=%.2e$' % (int_max, int_period)
    fig.suptitle(title)

    plt.show()
