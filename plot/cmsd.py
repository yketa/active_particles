"""
Module cmsd plots transverse and longitudinal collective mean square
displacements for different lag times.

Input files in simulation directory must follow the active_particles.naming.Ctt
(transverse collective mean square displacements) and
active_particles.naming.Cll (longitudinal collective mean square displacements)
standards.

Environment modes
-----------------
FITTING_LINE : bool
    Display adjustable fitting line on plots.
    DEFAULT: False
DIVIDE_BY_MAX : bool
    Divide collective mean square displacement multiplied by squared wave
    vector norm by its maximum for each lag-time.
    DEFAULT: False

Environment parameters
----------------------
DATA_DIRECTORY : string
	Data directory.
	DEFAULT: current working directory
PARAMETERS_FILE : string
	Simulation parameters file.
	DEFAULT: DATA_DIRECTORY/active_particles.naming.parameters_file
INITIAL_FRAME : int
	Frame to consider as initial.
	NOTE: INITIAL_FRAME < 0 will be interpreted as the initial frame being
	      the middle frame of the simulation.
	DEFAULT: -1
INTERVAL_MAXIMUM : int
	Maximum number of intervals of length dt considered for the calculation.
	DEFAULT: 1
N_CASES : int
    Number of boxes in each direction for which displacement grids have been
    computed.
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
R_CUT : float
    Wave length Gaussian cut-off radii in units of average particle separation.
    DEFAULT: active_particles.plot.cmsd._r_cut
DT_MIN : int
    Minimum lag time.
    NOTE: DT_MIN=None does not set any lower bound.
    DEFAULT: None
DT_MAX : int
    Maximum lag time.
    NOTE: DT_MAX=None does not set any upper bound.
    DEFAULT: None
Y_MIN : float
    Minimum collective mean square displacement multiplied by squared wave
    vector norm value.
    NOTE: Y_MIN=None does not set any lower bound.
    DEFAULT: None
Y_MAX : float
    Maximum collective mean square displacement multiplied by squared wave
    vector norm value.
    NOTE: Y_MAX=None does not set any upper bound.
    DEFAULT: None
R_MIN : float
    Minimum wave length over average particle separation value.
    NOTE: R_MIN=None does not set any lower bound.
    DEFAULT: None
R_MAX : float
    Maximum wave length over average particle separation value.
    NOTE: R_MAX=None does not set any upper bound.
    DEFAULT: None
FONT_SIZE : int
    Plot font size.
    DEFAULT: active_particles.plot.cmsd._font_size
RATIO_LEGEND : float
    Width ratio between legend and figure.
    DEFAULT: active_particles.plot.cmsd._ratio_legend
NCOL_LEGEND : int
    Number of columns for the legend.
    DEFAULT: active_particles.plot.cmsd._ncol_legend
WSPACE : float
    Plots width space.
    DEFAULT: active_particles.plot.cmsd._wspace
HSPACE : float
    Plots height space.
    DEFAULT: active_particles.plot.cmsd._hspace
COLORMAP : string
    Plot colormap.
    DEFAULT: active_particles.plot.cmsd._colormap
SLOPE_CTT [FITTING_LINE mode] : float
    Initial slope for fitting line slider for transversal CMSD plot.
    DEFAULT: active_particles.plot.cmsd._slope0
SLOPE_MIN_CTT [FITTING_LINE mode] : float
    Minimum slope for fitting line slider for transversal CMSD plot.
    DEFAULT: active_particles.plot.cmsd._slope_min
SLOPE_MAX_CTT [FITTING_LINE mode] : float
    Maximum slope for fitting line slider for transversal CMSD plot.
    DEFAULT: active_particles.plot.cmsd._slope_max
SLOPE_CLL [FITTING_LINE mode] : float
    Initial slope for fitting line slider for longitudinal CMSD plot.
    DEFAULT: active_particles.plot.cmsd._slope0
SLOPE_MIN_CLL [FITTING_LINE mode] : float
    Minimum slope for fitting line slider for longitudinal CMSD plot.
    DEFAULT: active_particles.plot.cmsd._slope_min
SLOPE_MAX_CLL [FITTING_LINE mode] : float
    Maximum slope for fitting line slider for longitudinal CMSD plot.
    DEFAULT: active_particles.plot.cmsd._slope_max
FITTING_LINE_SLIDER [FITTING_LINE mode] : bool
    Display fitting line slider.
    DEFAULT: True
FITTING_LINE_LEGEND [FITTING_LINE mode] : bool
    Display fitting line legend.
    DEFAULT: True
"""

import active_particles.naming as naming

from active_particles.init import get_env

from os import getcwd
from os import environ as envvar
if __name__ == '__main__': envvar['SHOW'] = 'True'
from os.path import join as joinpath

from active_particles.quantities import nD0_active
from active_particles.maths import g2Dto1Dgrid, FFT2Dfilter

from active_particles.plot.plot import list_colormap
from active_particles.plot.mpl_tools import FittingLine

from math import ceil

import numpy as np

import pickle

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# DEFAULT VARIABLES

_r_cut = 0  # default wave length Gaussian cut-off radii in units of average particle separation

_font_size = 10 # default plot font size

_ratio_legend = 7   # default width ratio between legend and figure
_ncol_legend = 1    # default number of columns for the legend

_wspace = 0.4   # default plots width space
_hspace = 0.05  # default plots height space

_colormap = 'jet'   # default plot colormap

_slope0 = 0     # default initial slope for fitting line slider
_slope_min = -2 # default minimum slope for fitting line slider
_slope_max = 2  # maximum slope for fitting line slider

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of intervals of length dt considered in correlations calculations

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    av_p_sep = parameters['box_size']/np.sqrt(parameters['N'])  # avergae particles separation

    nD0 = nD0_active(parameters['N'], parameters['vzero'], parameters['dr'],
        parameters['box_size'])

    Ncases = get_env('N_CASES', default=ceil(np.sqrt(parameters['N'])),
		vartype=int)	# number of boxes in each direction for which shear strain and displacement vorticity grids have been computed

    Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame

    box_size = get_env('BOX_SIZE', default=parameters['box_size'],
		vartype=float)									# size of the square box to consider
    centre = (get_env('X_ZERO', default=0, vartype=float),
		get_env('Y_ZERO', default=0, vartype=float))	# centre of the box

    r_cut = get_env('R_CUT', default=_r_cut, vartype=float) # wave length Gaussian cut-off radii in units of average particle separation

    divide_by_max = get_env('DIVIDE_BY_MAX', default=False, vartype=bool)   # divide collective mean square displacement multiplied by squared wave vector norm by its maximum for each lag-time

    dt_min = get_env('DT_MIN', vartype=int) # minimum lag time
    dt_max = get_env('DT_MAX', vartype=int) # maximum lag time

    # NAMING

    attributes = {'density': parameters['density'],
        'vzero': parameters['vzero'], 'dr': parameters['dr'],
        'N': parameters['N'], 'init_frame': init_frame,
        'int_max': int_max, 'Ncases': Ncases, 'box_size': box_size,
        'x_zero': centre[0], 'y_zero': centre[1]}   # attributes displayed in filenames
    naming_Ctt = naming.Ctt()                       # Ctt naming object
    naming_Cll = naming.Cll()                       # Cll naming object

    # PLOT PARAMETERS

    y_min = get_env('Y_MIN', vartype=float) # minimum collective mean square displacement multiplied by squared wave vector norm value
    y_max = get_env('Y_MAX', vartype=float) # maximum collective mean square displacement multiplied by squared wave vector norm value

    r_min = get_env('R_MIN', vartype=float) # minimum wave length over average particle separation value
    r_max = get_env('R_MAX', vartype=float) # maximum wave legnth over average particle separation value

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=int)   # plot font size
    mpl.rcParams.update({'font.size': font_size})

    ratio_legend = get_env('RATIO_LEGEND', default=_ratio_legend,
        vartype=float)                                                      # width ratio between legend and figure
    ncol_legend = get_env('NCOL_LEGEND', default=_ncol_legend, vartype=int) # number of columns for the legend

    wspace = get_env('WSPACE', default=_wspace, vartype=float)  # plots width space
    hspace = get_env('HSPACE', default=_hspace, vartype=float)  # plots height space

    colormap = get_env('COLORMAP', default=_colormap)   # plot colormap

    # DATA

    files_Ctt = naming_Ctt.get_files(directory=data_dir, **attributes)  # Ctt files corresponding to parameters
    dt_list_Ctt = np.array(list(map(
        lambda file: naming_Ctt.get_data(file, 'dt'),
        files_Ctt))).flatten()                                          # list of lag times corresponding to files

    files_Cll = naming_Cll.get_files(directory=data_dir, **attributes)  # Cll files corresponding to parameters
    dt_list_Cll = np.array(list(map(
        lambda file: naming_Ctt.get_data(file, 'dt'),
        files_Cll))).flatten()                                          # list of lag times corresponding to files

    files, dt_list = [], []
    for file_Ctt, dt in zip(files_Ctt, dt_list_Ctt):
        if dt in dt_list_Cll\
            and (dt_min == None or dt >= dt_min)\
            and (dt_max == None or dt <= dt_max):   # looking for transversal and longitudinal collective mean square displacements files with same lag times with lag times in considered interval
            files += [
                (file_Ctt, files_Cll[dt_list_Cll.tolist().index(dt)])]
            dt_list += [dt]

    Ctt, Cll = {}, {}           # hash tables of the 2D mean square norms of cross and dot products of normalised wave vectors with displacement grids Fourier transform with lag times as keys
    wave_vectors, k = {}, {}    # hash tables of 2D wave vectors and wave vector norms with lag times as keys
    for (file_Ctt, file_Cll), dt in zip(files, dt_list):
        with open(joinpath(data_dir, file_Ctt), 'rb') as Ctt_dump_file,\
            open(joinpath(data_dir, file_Cll), 'rb') as Cll_dump_file:
            wave_vectors[dt], Ctt[dt], _ = pickle.load(Ctt_dump_file)
            _, Cll[dt], _ = pickle.load(Cll_dump_file)
            k[dt] = np.sqrt(np.sum(wave_vectors[dt]**2, axis=-1))

    dt_list.sort()  # sorted lag times

    # PLOT

    colors = list_colormap(dt_list, colormap=colormap)  # hash table of line colors with lag times as keys

    fig = plt.figure()
    fig.set_size_inches(30, 15)
    fig.subplots_adjust(wspace=wspace)
    fig.subplots_adjust(hspace=hspace)

    gs = GridSpec(1, 3, width_ratios=[1, 1, 2/ratio_legend])

    ax_Ctt = plt.subplot(gs[0]) # Ctt plot axis
    ax_Cll = plt.subplot(gs[1]) # Cll plot axis
    for dt in dt_list:

        Ctt1D = g2Dto1Dgrid(
            FFT2Dfilter(Ctt[dt], wave_vectors=wave_vectors[dt]
                ).gaussian_filter(np.sqrt(2)*av_p_sep*r_cut).signalFFT,
            k[dt])
        x_Ctt = 2*np.pi/(Ctt1D[1:, 0]*av_p_sep)
        y_Ctt = Ctt1D[1:, 1]*(Ctt1D[1:, 0]**2)

        Cll1D = g2Dto1Dgrid(
            FFT2Dfilter(Cll[dt], wave_vectors=wave_vectors[dt]
                ).gaussian_filter(np.sqrt(2)*av_p_sep*r_cut).signalFFT,
            k[dt])
        x_Cll = 2*np.pi/(Cll1D[1:, 0]*av_p_sep)
        y_Cll = Cll1D[1:, 1]*(Cll1D[1:, 0]**2)

        if divide_by_max:
            ax_Ctt.loglog(x_Ctt, y_Ctt/max(y_Ctt), color=colors[dt])
            ax_Cll.loglog(x_Cll, y_Cll/max(y_Cll), color=colors[dt])
        else:
            ax_Ctt.loglog(x_Ctt, y_Ctt, color=colors[dt])
            ax_Cll.loglog(x_Cll, y_Cll, color=colors[dt])

    ax_Ctt.set_xlabel(r'$\lambda/a \equiv 2\pi/ka$')
    ax_Ctt.set_xlim(left=r_min, right=r_max)
    ax_Ctt.set_ylabel(r'$C^{\perp}(k) \times k^2 \equiv$'
        + r'$\left<||\vec{k}\wedge\tilde{\vec{u}}(\vec{k})||^2\right>$')
    ax_Ctt.set_ylim(bottom=y_min, top=y_max)

    ax_Cll.set_xlabel(r'$\lambda/a \equiv 2\pi/ka$')
    ax_Cll.set_xlim(left=r_min, right=r_max)
    ax_Cll.set_ylabel(r'$C^{||}(k) \times k^2 \equiv$'
        + r'$\left<||\vec{k}\cdot\tilde{\vec{u}}(\vec{k})||^2\right>$')
    ax_Cll.set_ylim(bottom=y_min, top=y_max)

    if get_env('FITTING_LINE', default=False, vartype=bool):    # FITTING LINE mode

        slope0_Ctt = get_env('SLOPE_CTT',
            default=_slope0, vartype=float)     # initial slope for fitting line slider for transversal CMSD plot
        slope_min_Ctt = get_env('SLOPE_MIN_CTT',
            default=_slope_min, vartype=float)  # minimum slope for fitting line slider for transversal CMSD plot
        slope_max_Ctt = get_env('SLOPE_MAX_CTT',
            default=_slope_max, vartype=float)  # maximum slope for fitting line slider for transversal CMSD plot

        slope0_Cll = get_env('SLOPE_CLL',
            default=_slope0, vartype=float)     # initial slope for fitting line slider for longitudinal CMSD plot
        slope_min_Cll = get_env('SLOPE_MIN_CLL',
            default=_slope_min, vartype=float)  # minimum slope for fitting line slider for longitudinal CMSD plot
        slope_max_Cll = get_env('SLOPE_MAX_CLL',
            default=_slope_max, vartype=float)  # maximum slope for fitting line slider for longitudinal CMSD plot

        fitting_line_slider = get_env('FITTING_LINE_SLIDER',
            default=True, vartype=float)    # display fitting line slider
        fitting_line_legend = get_env('FITTING_LINE_LEGEND',
            default=True, vartype=float)    # display fitting line legend

        fitting_line_Ctt = FittingLine(ax_Ctt,
            slope0_Ctt, slope_min_Ctt, slope_max_Ctt,
            x_fit='(\lambda/a)', y_fit='C^{\perp}(k) \\times k^2',
            slider=fitting_line_slider, legend=fitting_line_legend)
        fitting_line_Cll = FittingLine(ax_Cll,
            slope0_Cll, slope_min_Cll, slope_max_Cll,
            x_fit='(\lambda/a)', y_fit='C^{||}(k) \\times k^2',
            slider=fitting_line_slider, legend=fitting_line_legend)

    leg = plt.subplot(gs[2])    # legend axis
    leg.axis('off')
    legend = [mpatches.Patch(color='none',
        label=r'$nD_0 = \frac{Nv_0^2}{2\nu_r L^2} = %.2e$' % nD0)]
    legend += list(map(lambda dt: Line2D([0], [0], color=colors[dt],
        label=r'$nD_0\Delta t = %.2e$'
        % (dt*parameters['time_step']*parameters['period_dump']*nD0)),
        dt_list))
    leg.legend(handles=legend, loc='center', ncol=ncol_legend)

    title = r'$N=%.2e, \phi=%1.2f,$' % (parameters['N'], parameters['density'])
    title += r'$\tilde{v}=%.2e, \tilde{\nu}_r=%.2e$' % (
            parameters['vzero'], parameters['dr'])
    title += r'$, L=%.2e, x_0=%.2e, y_0=%.2e$' % (box_size, *centre)
    title += '\n'
    title += r'$S_{init}=%.2e, S_{max}=%.2e,$' % (init_frame, int_max)
    title += (r'$N_{cases}=%.2e, dL=%.2e a, r_{cut}=%.2e a$'
        % (Ncases, np.sqrt(parameters['N'])/Ncases, r_cut))
    fig.suptitle(title)

    plt.show()
