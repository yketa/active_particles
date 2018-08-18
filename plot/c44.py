"""
Module c44 plots strain correlations projected on cos(4 \\theta) for different
lag times, which we call C44.

Input files in simulation directory must follow the active_particles.naming.Css
(strain correlations computation in real and Fourier space) or the
active_particles.naming.Ctt and active_particles.naming.Cll (respectively for
transverse and longitudinal mean square displacements for strain correlations
computation from these quantities) standards.

Environment modes
-----------------
MODE : string
    Strain correlations computation mode.
     _________________________________________________________________
    | Mode    | Computation mode                                      |
    |_________|_______________________________________________________|
    | real    | computation in real space                             |
    |_________|_______________________________________________________|
    | fourier | computation in Fourier space                          |
    |_________|_______________________________________________________|
    | cmsd    | computation from collective mean square displacements |
    |_________|_______________________________________________________|
    DEFAULT: real
TEMPERATURE : bool
    Simulation directory corresponds to a thermal Brownian simulation.
     ________________________________________________
    | Mode  | Simulation type  | D_0                 |
    |_______|__________________|_____________________|
    | True  | thermal Brownian | 2 k_B T / a \\gamma |
    |_______|__________________|_____________________|
    | False | active Brownian  | v_0^2 / 2 \\nu_r    |
    |_______|__________________|_____________________|
    NOTE: a \\gamma is the Brownian dumping coefficient.
    DEFAULT: False
FITTING_LINE : bool
    Display adjustable fitting line on plot.
    DEFAULT: False
DIVIDE_BY_MAX : bool
    Divide C44 by its maximum for each lag-time and display these maxima in
    graph inset.
    DEFAULT: False
DIVIDE_BY_CHI : bool
    Divide C44 by its susceptibility for each lag-time and display these
    susceptibilities in graph inset.
    NOTE: if DIVIDE_BY_CHI=True, DIVIDE_BY_MAX is set to False.
    DEFAULT: False
LINEAR_INTERPOLATION : bool
	Get value on grid by linear interpolation of neighbouring grid boxes.
	DEFAULT: True

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
R_CUT ['real' mode] : float
    Cut-off radius for coarse graining function in terms of particles' mean
    radius.
    DEFAULT: active_particles.analysis.css._r_cut
SIGMA ['real' mode ]: float
    Length scale of the spatial extent of the coarse graining function.
    DEFAULT: cut-off radius
R_CUT_FOURIER ['fourier' or 'cmsd' mode] : float
    Wave length cut-off radius.
    DEFAULT: active_particles.plot.c44._r_cut_fourier
N_CASES : int
    Number of boxes in each direction for which shear strain and displacement
    vorticity grids have been computed.
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
Y_MIN : float
    Minimum C44 value.
    DEFAULT: active_particles.analysis.css._y_min_c44
Y_MAX : float
    Maximum C44 value.
    DEFAULT: active_particles.analysis.css._y_max_c44
R_MIN : float
    Minimum radius over average particle separation value.
    DEFAULT: active_particles.analysis.css._r_min_c44
R_MAX : float
    Maximum radius over average particle separation value.
    DEFAULT: active_particles.analysis.css._r_max_c44
POINTS_X : int
    Number of radii at which to evaluate integrated strain correlation.
    DEFAULT: active_particles.analysis.css._points_x_c44
POINTS_THETA : int
    Number of angles to evaluate integrated strain correlation.
    DEFAULT: active_particles.analysis.css._points_theta_c44
SMOOTH ['fourier' mode] : float
	C44 smoothing length scale.
	DEFAULT: 0
R_MIN_CHI [DIVIDE_BY_CHI mode] : float
    Minimum radius for susceptibility integration.
    DEFAULT: active_particles.plot.chi_msd._r_min
R_MAX_CHI [DIVIDE_BY_CHI mode] : float
    Maximum radius for susceptibility integration.
    DEFAULT: active_particles.plot.chi_msd._r_max
DT_MIN : int
    Minimum lag time.
    NOTE: if DT_MIN=None, no minimum is taken.
    DEFAULT: None.
DT_MAX : int
    Maximum lag time.
    NOTE: if DT_MAX=None, no maximum is taken.
    DEFAULT: None.
FONT_SIZE : int
    Plot font size.
    DEFAULT: active_particles.plot.c44._font_size
MARKER_SIZE : int
    Plot marker size.
    DEFAULT: active_particles.plot.c44._marker_size
RATIO_LEGEND : float
    Width ratio between legend and figure.
    DEFAULT: active_particles.plot.c44._ratio_legend
NCOL_LEGEND : int
    Number of columns for the legend.
    DEFAULT: active_particles.plot.c44._ncol_legend
WSPACE : float
    Plots width space.
    DEFAULT: active_particles.plot.c44._wspace
HSPACE : float
    Plots height space.
    DEFAULT: active_particles.plot.c44._hspace
WIDTH_INSET [DIVIDE_BY_MAX mode] : float
    Maximum C44 inset width in percentage of graph width.
    DEFAULT: active_particles.plot.c44._width_inset
HEIGHT_INSET [DIVIDE_BY_MAX mode] : float
    Maximum C44 inset height in percentage of graph height.
    DEFAULT: active_particles.plot.c44._height_inset
COLORMAP : string
    Plot colormap.
    DEFAULT: active_particles.plot.c44._colormap
SLOPE [FITTING_LINE mode] : float
	Initial slope for fitting line.
	DEFAULT: active_particles.analysis.css._slope0_c44
SLOPE_MIN [FITTING_LINE mode] : float
	Minimum slope for fitting line.
	DEFAULT: active_particles.analysis.css._slope_min_c44
SLOPE_MAX [FITTING_LINE mode] : float
	Maximum slope for fitting line.
	DEFAULT: active_particles.analysis.css._slope_max_c44
"""

import active_particles.naming as naming

from active_particles.init import get_env

from os import getcwd
from os import environ as envvar
if __name__ == '__main__': envvar['SHOW'] = 'True'
from os.path import join as joinpath

from active_particles.maths import wave_vectors_2D
from active_particles.quantities import nD0_active, nD0_thermal

from active_particles.analysis.css import Css2DtoC44, StrainCorrelations,\
    _slope0_c44 as _slope0, _slope_min_c44 as _slope_min,\
    _slope_max_c44 as _slope_max,\
    _points_x_c44 as _points_x, _points_theta_c44 as _points_theta,\
    _y_min_c44 as _y_min, _y_max_c44 as _y_max,\
    _r_min_c44 as _r_min, _r_max_c44 as _r_max,\
    _r_cut
from active_particles.analysis.ctt import StrainCorrelationsCMSD
from active_particles.analysis.cuu import c1Dtochi
from active_particles.plot.plot import list_colormap
from active_particles.plot.mpl_tools import FittingLine
from active_particles.plot.chi_msd import _r_min as _r_min_chi,\
    _r_max as _r_max_chi

from math import ceil

import numpy as np

import pickle

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# DEFAULT VARIABLES

_r_cut_fourier = 0  # default cut-off radius for coarse graining function

_font_size = 10     # default plot font size
_marker_size = 10   # default plot marker size

_ratio_legend = 7   # default width ratio between legend and figure
_ncol_legend = 1    # default number of columns for the legend

_wspace = 0.2   # default plots width space
_hspace = 0.05  # default plots height space

_colormap = 'jet'   # default plot colormap

_width_inset = 30   # default maximum C44 inset width in percentage of graph width
_height_inset = 30  # default maximum C44 inset height in percentage of graph height

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLES DEFINITIONS

    mode = get_env('MODE', default='real')  # strain correlations computation mode

    if mode == 'real': naming_Css = naming.Css(from_ft=False)       # Css naming object
    elif mode == 'fourier': naming_Css = naming.Css(from_ft=True)   # Css naming object
    elif mode == 'cmsd':
        naming_Ctt = naming.Ctt()                                   # Ctt naming object
        naming_Cll = naming.Cll()                                   # Cll naming object
    else: raise ValueError('Mode %s is not known.' % mode)          # mode is not known

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of intervals of length dt considered in correlations calculations

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    av_p_sep = parameters['box_size']/np.sqrt(parameters['N'])  # avergae particles separation

    r_cut = parameters['a']*get_env('R_CUT', default=_r_cut, vartype=float) # cut-off radius for coarse graining function
    sigma = get_env('SIGMA', default=r_cut, vartype=float)                  # length scale of the spatial extent of the coarse graining function

    r_cut_fourier = get_env('R_CUT_FOURIER',
        default=_r_cut_fourier, vartype=float)  # wave length cut-off radius when computing strain correlations from Fourier transform components

    Ncases = get_env('N_CASES', default=ceil(np.sqrt(parameters['N'])),
		vartype=int)	# number of boxes in each direction for which shear strain and displacement vorticity grids have been computed

    Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame

    box_size = get_env('BOX_SIZE', default=parameters['box_size'],
		vartype=float)									# size of the square box to consider
    centre = (get_env('X_ZERO', default=0, vartype=float),
		get_env('Y_ZERO', default=0, vartype=float))	# centre of the box

    divide_by_max = get_env('DIVIDE_BY_MAX', default=False, vartype=bool)   # divide C44 by its maximum for each lag-time

    divide_by_chi = get_env('DIVIDE_BY_CHI', default=False, vartype=bool)   # divide C44 by its susceptibility for each lag time
    if divide_by_chi: divide_by_max = False
    r_min_chi = get_env('R_MIN_CHI', default=_r_min_chi, vartype=float)     # minimum radius for susceptibility integration
    r_max_chi = get_env('R_MAX_CHI', default=_r_max_chi, vartype=float)     # maximum radius for susceptibility integration

    dt_min = get_env('DT_MIN', vartype=int) # minimum lag time
    dt_max = get_env('DT_MAX', vartype=int) # maximum lag time

    # NAMING

    attributes = {'density': parameters['density'], 'N': parameters['N'],
        'init_frame': init_frame, 'int_max': int_max, 'Ncases': Ncases,
        'r_cut': r_cut, 'sigma': sigma, 'box_size': box_size,
        'x_zero': centre[0], 'y_zero': centre[1]}   # attributes displayed in file names

    # PLOT PARAMETERS

    y_min = get_env('Y_MIN', default=_y_min, vartype=float) # minimum C44 value
    y_max = get_env('Y_MAX', default=_y_max, vartype=float) # maximum C44 value

    r_min = get_env('R_MIN', default=_r_min, vartype=float) # minimum radius over average particle separation value
    r_max = get_env('R_MAX', default=_r_max, vartype=float) # maximum radius over average particle separation value

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=int)       # plot font size
    marker_size = get_env('MARKER_SIZE', default=_marker_size, vartype=int) # plot marker size
    mpl.rcParams.update({'font.size': font_size,
        'lines.markersize': marker_size})

    ratio_legend = get_env('RATIO_LEGEND', default=_ratio_legend,
        vartype=float)                                                      # width ratio between legend and figure
    ncol_legend = get_env('NCOL_LEGEND', default=_ncol_legend, vartype=int) # number of columns for the legend

    wspace = get_env('WSPACE', default=_wspace, vartype=float)  # plots width space
    hspace = get_env('HSPACE', default=_hspace, vartype=float)  # plots height space

    colormap = get_env('COLORMAP', default=_colormap)   # plot colormap

    width_inset = get_env('WIDTH_INSET',
        default=_width_inset, vartype=float)    # maximum C44 inset width in percentage of graph width
    height_inset = get_env('HEIGHT_INSET',
        default=_height_inset, vartype=float)   # maximum C44 inset height in percentage of graph height

    # CALCULATION PARAMETERS

    points_x = get_env('POINTS_X', default=_points_x, vartype=int)              # number of radii at which to evaluate integrated strain correlation
    points_theta = get_env('POINTS_THETA', default=_points_theta, vartype=int)  # number of angles to evaluate integrated strain correlation

    # CALCULATION

    toC44 = Css2DtoC44(box_size, points_x, points_theta,
        av_p_sep*r_min, av_p_sep*r_max,
        linear_interpolation=get_env(
            'LINEAR_INTERPOLATION', default=True, vartype=bool))    # Css2D to C44 conversion object
    C44 = {}                                                        # hash table of tables of C44 at X with lag times as keys

    if mode == 'real':

        files = naming_Css.get_files(directory=data_dir, **attributes)          # files corresponding to parameters
        dt_list = np.array(list(map(
            lambda file: naming_Css.get_data(file, 'dt'), files))).flatten()    # list of lag times corresponding to files

        for file, dt in zip(files, dt_list):
            with open(joinpath(data_dir, file), 'rb') as Css_dump_file:
                _, Css2D = pickle.load(Css_dump_file)
            C44[dt] = toC44.get_C44(Css2D)

    elif mode == 'fourier':

        smooth = get_env('SMOOTH', default=0, vartype=float)    # C44 smoothing length scale

        files = naming_Css.get_files(directory=data_dir, **attributes)          # files corresponding to parameters
        dt_list = np.array(list(map(
            lambda file: naming_Css.get_data(file, 'dt'), files))).flatten()    # list of lag times corresponding to files

        wave_vectors = wave_vectors_2D(Ncases, Ncases, d=box_size/Ncases)   # wave vectors at which Fourier transform was calculated

        for file, dt in zip(files, dt_list):
            with open(joinpath(data_dir, file), 'rb') as Css_dump_file:
                FFTsgridsqnorm = pickle.load(Css_dump_file)
            C44[dt] = toC44.get_C44(
                StrainCorrelations(wave_vectors, FFTsgridsqnorm)
                .strain_correlations(r_cut=r_cut_fourier),
                smooth=smooth)

    elif mode == 'cmsd':

        files_Ctt = naming_Ctt.get_files(directory=data_dir, **attributes)  # files corresponding to parameters
        dt_list_Ctt = np.array(list(map(
            lambda file: naming_Ctt.get_data(file, 'dt'),
            files_Ctt))).flatten()                                          # list of lag times corresponding to files

        files_Cll = naming_Cll.get_files(directory=data_dir, **attributes)  # files corresponding to parameters
        dt_list_Cll = np.array(list(map(
            lambda file: naming_Cll.get_data(file, 'dt'),
            files_Cll))).flatten()                                          # list of lag times corresponding to files

        files, dt_list = [], []
        for file_Ctt, dt in zip(files_Ctt, dt_list_Ctt):
            if dt in dt_list_Cll:   # looking for transversal and longitudinal collective mean square displacements files with same lag times
                files += [
                    (file_Ctt, files_Cll[dt_list_Cll.tolist().index(dt)])]
                dt_list += [dt]

        wave_vectors = wave_vectors_2D(Ncases, Ncases, d=box_size/Ncases)   # wave vectors at which Fourier transform was calculated

        for (file_Ctt, file_Cll), dt in zip(files, dt_list):
            with open(joinpath(data_dir, file_Ctt), 'rb') as Ctt_dump_file,\
                open(joinpath(data_dir, file_Cll), 'rb') as Cll_dump_file:
                _, k_cross_FFTugrid2D_sqnorm, _ = pickle.load(Ctt_dump_file)
                _, k_dot_FFTugrid2D_sqnorm, _ = pickle.load(Cll_dump_file)
            C44[dt] = toC44.get_C44(
                StrainCorrelationsCMSD(wave_vectors,
                    k_cross_FFTugrid2D_sqnorm, k_dot_FFTugrid2D_sqnorm)
                .strain_correlations(r_cut=av_p_sep*r_cut_fourier))

    dt_list = [dt for dt in sorted(dt_list)
        if (dt_min == None or dt >= dt_min)
        and (dt_max == None or dt <= dt_max)]

    # PLOT

    colors = list_colormap(dt_list, colormap=colormap)  # hash table of line colors with lag times as keys

    fig = plt.figure()
    fig.set_size_inches(30, 30)
    fig.subplots_adjust(wspace=wspace)
    fig.subplots_adjust(hspace=hspace)

    gs = GridSpec(1, 2, width_ratios=[1, 1/ratio_legend])
    leg = plt.subplot(gs[1])    # legend axis
    ax = plt.subplot(gs[0])     # plot axis
    leg.axis('off')

    ax.set_xlabel(r'$r/a$' + ' ' + r'$(a = L/\sqrt{N})$')
    if divide_by_max:
        ax.set_ylabel(
            r'$C_4^4(r)/max(C_4^4(r), %.2e \leq r \leq %.2e)$' % (r_min, r_max)
            + '\n' + r'$C_4^4(r) = \frac{1}{\pi}\int_0^{2\pi}d\theta$'
            + ' ' + r'$C_{\epsilon_{xy}\epsilon_{xy}}(r, \theta)$'
            + ' ' + r'$\cos4\theta$')
    elif divide_by_chi:
        ax.set_ylabel(
            r'$C_4^4(r)/\chi$'
            + '\n' + r'$C_4^4(r) = \frac{1}{\pi}\int_0^{2\pi}d\theta$'
            + ' ' + r'$C_{\epsilon_{xy}\epsilon_{xy}}(r, \theta)$'
            + ' ' + r'$\cos4\theta$')
    else:
        ax.set_ylabel(r'$C_4^4(r) = \frac{1}{\pi}\int_0^{2\pi}d\theta$'
            + ' ' + r'$C_{\epsilon_{xy}\epsilon_{xy}}(r, \theta)$'
            + ' ' + r'$\cos4\theta$')
    ax.set_xlim([r_min, r_max])
    ax.set_ylim([y_min, y_max])
    ax.set_yscale('log')
    ax.set_xscale('log')

    if get_env('FITTING_LINE', default=False, vartype=bool):    # FITTING LINE mode

        slope0 = get_env('SLOPE', default=_slope0, vartype=float)           # initial slope for fitting line slider
        slope_min = get_env('SLOPE_MIN', default=_slope_min, vartype=float) # minimum slope for fitting line slider
        slope_max = get_env('SLOPE_MAX', default=_slope_max, vartype=float) # maximum slope for fitting line slider

        fitting_line = FittingLine(ax, slope0, slope_min, slope_max,
            x_fit='(r/a)', y_fit='C_4^4(r/a)')

    if get_env('TEMPERATURE', default=False, vartype=bool): # TEMPERATURE mode
        nD0 = nD0_thermal(parameters['N'], parameters['kT'],
            parameters['damp_bro']*parameters['a'], parameters['box_size'])
        legend0 = [mpatches.Patch(color='none',
            label=r'$nD_0 = \frac{2 k_B T N}{\lambda a L^2} = %.2e$' % nD0)]
    else:
        nD0 = nD0_active(parameters['N'], parameters['vzero'],
            parameters['dr'], parameters['box_size'])
        legend0 = [mpatches.Patch(color='none', label=
            r'$nD_0 = \frac{Nv_0^2}{2\nu_r L^2} = %.2e$' % nD0)]
    legend0 += list(map(lambda dt: Line2D([0], [0], color=colors[dt],
        label=r'$nD_0\Delta t = %.2e$'
        % (dt*parameters['time_step']*parameters['period_dump']*nD0)),
        dt_list))
    legend0 += [Line2D([0], [0], lw=0, label='')]
    leg.legend(handles=legend0, loc='center', ncol=ncol_legend)

    if divide_by_max or divide_by_chi:                                          # DIVIDE_BY_MAX or DIVIDE_BY_CHI mode
        inset_ax = inset_axes(ax, loc=1,
            width='%e' % width_inset + '%', height='%e' % height_inset + '%')   # inset axes
        inset_ax.set_yscale('log')
        inset_ax.set_xscale('log')
        inset_ax.set_xlabel(r'$nD_0\Delta t$')
        if divide_by_max: inset_ax.set_ylabel(r'$max(C_4^4(r))$')               # DIVIDE_BY_MAX mode
        if divide_by_chi: inset_ax.set_ylabel(
            r'$\chi = \frac{2\pi}{L^2}\int_{r_{min}}^{r_{max}}dr$'
            + ' ' + r'$r$' + ' ' + r'$C_4^4(r)$')                               # DIVIDE_BY_CHI mode

    for dt in dt_list:
        if divide_by_max:   # DIVIDE_BY_MAX mode
            max_C44 = max(C44[dt][:, 1])
            C44[dt][:, 1] /= max_C44
            inset_ax.scatter(
                dt*parameters['time_step']*parameters['period_dump']*nD0,
                max_C44,
                color=colors[dt])
        if divide_by_chi:   # DIVIDE_BY_CHI mode
            chi_C44 = c1Dtochi(C44[dt], box_size, r_min=r_min, r_max=r_max)
            C44[dt][:, 1] /= chi_C44
            inset_ax.scatter(
                dt*parameters['time_step']*parameters['period_dump']*nD0,
                chi_C44,
                color=colors[dt])
        ax.plot(C44[dt][:, 0]/av_p_sep, C44[dt][:, 1], color=colors[dt])

    title = r'$N=%.2e, \phi=%1.2f,$' % (parameters['N'], parameters['density'])
    if get_env('TEMPERATURE', default=False, vartype=bool): # TEMPERATURE mode
        title += r'$kT=%.2e, k=%.2e$' % (parameters['kT'], parameters['k'])
    else:
        title += r'$\tilde{v}=%.2e, \tilde{\nu}_r=%.2e$' % (
            parameters['vzero'], parameters['dr'])
    title += r'$, L=%.2e, x_0=%.2e, y_0=%.2e$' % (box_size, *centre)
    if mode == 'fourier' and 'SMOOTH' in envvar:
        title += r'$, \sigma_{smooth}/a=%.2e$' % smooth
    title += '\n'
    title += r'$S_{init}=%.2e, S_{max}=%.2e,$' % (init_frame, int_max)
    title += (r'$N_{cases}=%.2e, dL=%.2e a$'
        % (Ncases, np.sqrt(parameters['N'])/Ncases))
    if mode == 'real':
        title += r'$, r_{cut}=%.2e, \sigma=%.2e$' % (r_cut, sigma)
    else:
        title += r'$, r_{cut}=%.2e$' % r_cut_fourier
    title += '\n'
    title += r'$N_r=%.2e, N_{\theta}=%.2e$' % (points_x, points_theta)
    if divide_by_chi:
        title += r'$, r_{min}=%.2e, r_{max}=%.2e$' % (r_min_chi, r_max_chi)
    fig.suptitle(title)

    plt.show()
