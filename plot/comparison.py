"""
Module comparison superimposes most probable local density, maximum
cooperativity, time of maximum cooperativity and ratio of transversal and
longitudinal correlations at time of maximum cooperativity, as functions of
the Péclet number, for different trajectories in the phase diagram (either
varying persistence time at fixed self-propelling velocity or varying
self-propelling velocity at fixed persistence time).

Simulation directories must follow the active_particles.naming.AHB2D naming
standard and input files in simulation directories must follow the
active_particles.naming.VarN standard (local densities) and either the
active_particles.naming.Cuu standard (for cooperativities from displacement
correlations), or active_particles.naming.Cww standard (for cooperativities
from displacement relative to centre of mass displacement correlations), or
active_particles.naming.Cdd standard (for cooperativities from displacement
norm correlations), or active_particles.naming.Cee (for cooperativities from
displacement direction correlations).

Environment modes
-----------------
CORRELATION : string
    Correlations from which to calculate cooperativities and extract
    longitudinal and transversal components.
     _____________________________________________________________
    | Mode | Correlations                                         |
    |______|______________________________________________________|
    | Cuu  | displacement                                         |
    |______|______________________________________________________|
    | Cww  | displacement relative to centre of mass displacement |
    |______|______________________________________________________|
    | Cdd  | displacement norm                                    |
    |______|______________________________________________________|
    | Cee  | displacement direction                               |
    |______|______________________________________________________|
    DEFAULT: Cuu
DRDT : bool
    Use the product of the rotation diffusion constant and lag time rather
    than the bare lag time.
    DEFAULT: True
SHOW : bool
    Show comparison plot.
    DEFAULT: True

Environment parameters
----------------------
DATA_DIRECTORY : string
    Data directory.
    DEFAULT: active_particles.naming.sim_directory
EXCLUDE : string
    Simulation directories in DATA_DIRECTORY to exclude from the plots.
    DEFAULT:
PARAMETERS_FILE : string
    Simulation parameters file name.
    DEFAULT: active_particles.naming.parameters_file
IMAGE_NAME : string
    Default image file name.
    DEFAULT: active_particles.plot.comparison._image_name
N : int
    Number of particles.
    DEFAULT: active_particles.plot.comparison._N
VARIABLES (mandatory) : string separated by ':'
    Trajectory variable in order of trajectory.
     _______________________________________
    | Mode    | Variable                    |
    |_________|_____________________________|
    | 'dr'    | Rotation diffusion constant |
    |_________|_____________________________|
    | 'vzero' | self-propelling velocity    |
    |_________|_____________________________|
VAR_MIN (mandatory) : float separated by ':'
    Minimum value of trajectory variable in order of trajectory.
VAR_MAX (mandatory) : float separated by ':'
    Maximum value of trajectory variable in order of trajectory.
VAR_C (mandatory) : float separated by ':'
    Trajectory variable transition value in order of trajectory.
    NOTE: A vertical line is displayed at the corresponding Péclet number.
FIXED_VAR (mandatory) : float separated by ':'
    Fixed variable value in order of trajectory.
DENSITIES (mandatory) : float separated by ':'
    Packing fractions of particles in order of trajectory.
INITIAL_FRAME_PHILOC : int
    Frame to consider as initial in local densities calculations.
    DEFAULT: active_particles.plot.pphiloc._init_frame
INTERVAL_MAXIMUM_PHILOC : int
    Maximum number of frames on which densities are calculated.
    DEFAULT: active_particles.plot.pphiloc._int_max
BOX_SIZE_PHILOC : float
    Length of the square boxes in which particles are counted.
    DEFAULT: active_particles.plot.pphiloc._box_size
N_CASES_PHILOC : int
    Number of boxes in each direction in which local densities are computed.
    DEFAULT: active_particles.plot.pphiloc._Ncases
N_BINS : int
    Number of bins for the local densities histogram.
    DEFAULT: active_particles.plot.pphiloc._Nbins
PHIMAX : float
    Maximum local density for the local densities histogram.
    DEFAULT: active_particles.plot.pphiloc._phimax
BOX_SIZE_COR : float
    Size of the square box which was considered for correlations.
    DEFAULT: simulation box size
X_ZERO : float
    Centre of the of the box x-coordinate for correlations.
    DEFAULT: 0
Y_ZERO : float
    Centre of the of the box y-coordinate for correlations.
    DEFAULT: 0
INITIAL_FRAME_COR : int
    Frame to consider as initial for correlations.
    DEFAULT: active_particles.plot.chi_msd._init_frame_cor
INTERVAL_MAXIMUM_COR : int
    Maximum number of intervals of length dt considered for correlations.
    DEFAULT: active_particles.plot.chi_msd._int_max_cor
N_CASES_COR : int
    Number of boxes in each direction with which the displacement grid is
    computed.
    DEFAULT: active_particles.plot.chi_msd._Ncases_cor
R_MIN : float
    Minimum radius for correlations integration.
    DEFAULT: active_particles.plot.chi_msd._r_min
R_MAX : float
    Maximum radius for correlations integration.
    DEFAULT: active_particles.plot.chi_msd._r_max
FONT_SIZE : int
    Plot font size.
    DEFAULT: active_particles.plot.comparison._font_size
MARKER_SIZE : int
    Plot marker size.
    DEFAULT: active_particles.plot.comparison._marker_size
COLORMAPS : string separated by ':'
    Plot colormaps to choose from.
    DEFAULT: active_particles.plot.comparison._colormaps
RATIO_LEGEND : float
    Width ratio between legend and figure.
    DEFAULT: active_particles.plot.comparison._ratio_legend
WSPACE : float
    Plots width space.
    DEFAULT: active_particles.plot.comparison._wspace
HSPACE : float
    Plots height space.
    DEFAULT: active_particles.plot.comparison._hspace
X_SCALE : string
    Plots x-scale.
    DEFAULT: active_particles.plot.comparison._x_scale
PHILOC_YS : string
    Most probable local density y-scale.
    DEFAULT: active_particles.plot.comparison._philoc_ys
CHI_YS : string
    Maximum cooperativity y-scale.
    DEFAULT: active_particles.plot.comparison._chi_ys
DT_YS : string
    Time of maximum cooperativity y-scale.
    DEFAULT: active_particles.plot.comparison._dt_ys
RATIOTL_YS : string
    Ratio of transversal and longitudinal correlations y-scale.
    DEFAULT: active_particles.plot.comparison._ratioTL_ys

Output
------
> Saves figure to IMAGE_NAME.
[SHOW mode]
> Displays plot.
"""

import active_particles.naming as naming

from active_particles.init import get_env, get_env_list

from os import environ as envvar
if __name__ == '__main__': envvar['SHOW'] = 'True'
from os.path import join as joinpath

from active_particles.plot.plot import list_colormap, list_markers,\
    list_linestyles
from active_particles.plot.pphiloc import Philoc,\
    _init_frame as _init_frame_philoc, _int_max as _int_max_philoc,\
    _box_size as _box_size_philoc, _Ncases as _Ncases_philoc,\
    _Nbins, _phimax
from active_particles.plot.chi_msd import ChiMsd,\
    _init_frame_cor, _int_max_cor, _Ncases_cor, _r_min, _r_max
from active_particles.plot.corlcort import CLCT

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

# DEFAULT VARIABLES

_N = int(1e5)   # default number of particles

_font_size = 25     # default plot font size
_marker_size = 20   # default plot marker size

_colormaps = ('cool', 'hot')    # default plot colormaps

_wspace = 0.4   # default plot width space
_hspace = 0.05  # default plot height space

_ratio_legend = 2   # default width ratio between graphs and legends

_x_scale = 'log'        # default plots x-scale
_philoc_ys = 'linear'   # default most probable local density y-scale
_chi_ys = 'log'         # default maximum cooperativity y-scale
_dt_ys = 'log'          # default time of maximum cooperativity y-scale
_ratioTL_ys = 'linear'  # default ratio of transversal and longitudinal correlations y-scale

_image_name = joinpath(get_env('HOME'), 'comparison.eps')   # default image file name

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLES DEFINITIONS

    var = get_env_list('VARIABLES')                         # plot variables
    var_min = get_env_list('VAR_MIN', vartype=float)        # minimum values of plot variable
    var_max = get_env_list('VAR_MAX', vartype=float)        # maximum values of plot variable
    var_c = get_env_list('VAR_C', vartype=float)            # plot variable transition value
    fixed_var = get_env_list('FIXED_VAR', vartype=float)    # values of fixed variables
    densities = get_env_list('DENSITIES', vartype=float)    # packing fractions of particles
    if not(len(var) == len(var_min) == len(var_max) == len(fixed_var)
        == len(var_c) == len(densities)):
        raise IndexError(
            'VARIABLES, VAR_MIN, VAR_MAX, VAR_C, FIXED_VAR and DENSITIES \
            must have equal lengths.')
    comparisons = len(var)                                  # number of trajectories to compare

    var_label = []      # variables labels
    fix_label = []      # fixed variable labels
    var_attribute = []  # variables attribute to be displayed in file names
    pe_func = []        # Péclet number as function of plot variable
    for index in range(comparisons):

        if var[index] == 'dr':
            var_label += ['\\tilde{\\nu}_r']
            fix_label += ['\\tilde{v}']
            var_attribute += [{'vzero': fixed_var[index]}]
            pe_func += [(lambda index: lambda x: fixed_var[index]/x)(index)]

        elif var[index] == 'vzero':
            var_label += ['\\tilde{v}']
            fix_label += ['\\tilde{\\nu}_r']
            var_attribute += [{'dr': fixed_var[index]}]
            pe_func += [(lambda index: lambda x: x/fixed_var[index])(index)]

        else: raise ValueError('Variable %s is not known.' % var[index])

    cor = get_env('CORRELATION', default='Cuu') # correlation variable

    if cor == 'Cuu':    # cooperativities from

        naming_cor = naming.Cuu()   # correlations naming object
        cor_name = 'C_{uu}'         # correlations name

    elif cor == 'Cww':

        naming_cor = naming.Cww()               # correlations naming object
        cor_name = 'C_{\\delta u \\delta u}'    # correlations name

    elif cor == 'Cdd':

        naming_cor = naming.Cdd()   # correlations naming object
        cor_name = 'C_{|u||u|}'     # correlations name

    elif cor == 'Cee':

        naming_cor = naming.Cee()           # correlations naming object
        cor_name = 'C_{\\hat{u}\\hat{u}}'   # correlations name

    else: raise ValueError('Correlation %s is not known.' % cor)    # correlation is not known

    data_dir = get_env('DATA_DIRECTORY', default=naming.sim_directory)  # data directory
    excluded_directories = get_env('EXCLUDE', default='')               # directories to exclude

    parameters_file = get_env('PARAMETERS_FILE',
        default=naming.parameters_file) # simulations parameters file name

    N = get_env('N', default=_N, vartype=int)                       # number of particles

    init_frame_philoc = get_env('INITIAL_FRAME_PHILOC',
        default=_init_frame_philoc, vartype=int)    # frame to consider as initial in local densities calculations
    int_max_philoc = get_env('INTERVAL_MAXIMUM_PHILOC',
        default=_int_max_philoc, vartype=int)       # maximum number of frames on which densities are calculated

    box_size_philoc = get_env('BOX_SIZE_PHILOC',
        default=_box_size_philoc, vartype=float)    # length of the square boxes in which particles are counted

    Ncases_philoc = get_env('N_CASES_PHILOC',
        default=_Ncases_philoc, vartype=int)    # number of boxes in each direction in which local densities are computed

    Nbins = get_env('N_BINS', default=_Nbins, vartype=int)      # number of bins for the local densities histogram
    phimax = get_env('PHIMAX', default=_phimax, vartype=float)  # maximum local density for the local densities histogram

    box_size_cor = get_env('BOX_SIZE_COR', vartype=float)   # size of the square box which was considered for correlations
    centre_cor = (get_env('X_ZERO', default=0, vartype=float),
	   get_env('Y_ZERO', default=0, vartype=float))         # centre of the box for correlations

    init_frame_cor = get_env('INITIAL_FRAME_COR',
        default=_init_frame_cor, vartype=int)                               # frame to consider as initial for correlations
    int_max_cor = get_env('INTERVAL_MAXIMUM_COR',
        default=_int_max_cor, vartype=int)                                  # maximum number of intervals of length dt considered for correlations
    Ncases_cor = get_env('N_CASES_COR', default=_Ncases_cor, vartype=int)   # number of boxes in each direction with which the displacement grid is computed

    r_min = get_env('R_MIN', default=_r_min, vartype=float) # minimum radius for correlations integration
    r_max = get_env('R_MAX', default=_r_max, vartype=float) # maximum radius for correlations integration

    # NAMING

    attributes_philoc = {'N': N, 'init_frame': init_frame_philoc,
        'int_max': int_max_philoc, 'Ncases': Ncases_philoc,
        'box_size': box_size_philoc}                        # attributes displayed in local densities file names
    attributes_cor = {'N': N, 'init_frame': init_frame_cor,
        'int_max': int_max_cor, 'Ncases': Ncases_cor, 'box_size': box_size_cor,
        'x_zero': centre_cor[0], 'y_zero': centre_cor[1]}   # attributes displayed in correlations file names
    naming_varN = naming.VarN()                     # varN naming object
    naming_simdir = naming.AHB2D()                  # simulation directory naming object

    # PLOT PARAMETERS

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=float)     # plot font size
    marker_size = get_env('MARKER_SIZE', default=_marker_size, vartype=int) # plot marker size
    mpl.rcParams.update({'font.size': font_size,
        'lines.markersize': marker_size})

    colormaps = get_env_list('COLORMAPS')                       # plot colormaps
    if colormaps == []: colormaps = _colormaps                  # no plot colormaps provided, use default
    while len(colormaps) < comparisons: colormaps += colormaps  # at least as much colormaps as trajectories to compare

    ratio_legend = get_env('RATIO_LEGEND',
        default=_ratio_legend, vartype=float)   # width ratio between graphs and legends

    wspace = get_env('WSPACE', default=_wspace, vartype=float)  # plots width space
    hspace = get_env('HSPACE', default=_hspace, vartype=float)  # plots height space

    x_scale = get_env('X_SCALE', default=_x_scale)          # plots x-scale
    philoc_ys = get_env('PHILOC_YS', default=_philoc_ys)    # most probable local density y-scale
    chi_ys = get_env('CHI_YS', default=_chi_ys)             # maximum cooperativity y-scale
    dt_ys = get_env('DT_YS', default=_dt_ys)                # time of maximum cooperativity y-scale
    ratioTL_ys = get_env('RATIOTL_YS', default=_ratioTL_ys) # ratio of transversal and longitudinal correlations y-scale

    multiply_with_dr = get_env('DRDT', default=True, vartype=bool)  # plot dr*dt rather than dt
    if multiply_with_dr:
        dt_label = r'$\tilde{\nu}_r \Delta t^*$'    # dt label
    else: dt_label = r'$\Delta t^*$'                # dt label

    # CALCULATION

    philoc = [] # local densities histogram calculators
    chimsd = [] # cooperativities calculators
    clct = []   # longtidunal and transversal correlations calculators
    for v, vmin, vmax, phi, vattribute\
        in zip(var, var_min, var_max, densities, var_attribute):

        philoc += [Philoc(data_dir, naming_simdir,
            {'density': phi, **vattribute, **attributes_philoc},
            parameters_file, v, vmin, vmax, excluded_dir=excluded_directories)]
        philoc[-1].calculate(naming_varN,
            {'density': phi, **vattribute, **attributes_philoc}, Nbins, phimax)

        chimsd += [ChiMsd(data_dir, naming_simdir,
            {'density': phi, **vattribute, **attributes_cor},
            parameters_file, v, vmin, vmax, excluded_dir=excluded_directories)]
        chimsd[-1].calculate(naming_cor,
            {'density': phi, **vattribute, **attributes_cor}, r_min, r_max,
            box_size=box_size_cor, multiply_with_dr=multiply_with_dr)

        clct += [CLCT(data_dir, naming_simdir,
            {'density': phi, **vattribute, **attributes_cor},
            parameters_file, v, vmin, vmax, excluded_dir=excluded_directories)]
        clct[-1].calculate(naming_cor,
            {'density': phi, **vattribute, **attributes_cor},
            multiply_with_dr=multiply_with_dr)
        clct[-1].calculate_max(chimsd[-1].dtmax)

    # PLOT

    colors = list(map(
        lambda philoc_traj, chimsd_traj, clct_traj, cmap:
            list_colormap(
                sorted(philoc_traj.var_list + chimsd_traj.var_list +
                clct_traj.var_list),
                colormap=cmap),
        *(philoc, chimsd, clct, colormaps)))    # plot colors hash tables
    markers = list(map(
        lambda philoc_traj, chimsd_traj, clct_traj:
            list_markers(
                sorted(philoc_traj.time_step_list + chimsd_traj.time_step_list
                + clct_traj.time_step_list)),
        *(philoc, chimsd, clct)))               # plot markers hash tables

    linestyles = list_linestyles(range(comparisons))    # Péclet transition values vertical lines

    fig = plt.figure()
    fig.set_size_inches(30, 30)
    fig.subplots_adjust(wspace=wspace)
    fig.subplots_adjust(hspace=hspace)
    fig.suptitle(
        r'$N=%.2e$' % N + '\n' + r'$[%s]:$' % cor_name
        + r'$S_{init}=%.2e, S_{max}=%.2e,$' % (init_frame_cor, int_max_cor)
        + r'$N_{cases}=%.2e, r_{min}=%.2e,$' % (Ncases_cor, r_min)
        + r'$r_{max}=%.2e$' % r_max + '\n'
        + r'$[\phi^*_{loc}]: S_{init}=%.2e,$' % init_frame_philoc
        + r'$S_{max}=%.2e, N_{cases}=%.2e,$' % (int_max_philoc, Ncases_philoc)
        + r'$r_{max}=%.2e$' % box_size_philoc)

    gs = GridSpec(4, 1 + comparisons,
        width_ratios=[1] + comparisons*[1/(comparisons*ratio_legend)])

    ax_philoc = plt.subplot(gs[0, 0])
    ax_philoc.set_xscale(x_scale)
    ax_philoc.set_yscale(philoc_ys)
    ax_philoc.set_ylabel(r'$\phi^*_{loc}$')

    ax_chi = plt.subplot(gs[1, 0])
    ax_chi.set_xscale(x_scale)
    ax_chi.set_yscale(chi_ys)
    ax_chi.set_ylabel(r'$\chi(\Delta t^*) = \frac{1}{L^2}$'
        + r'$\int_{r=r_{min}}^{r=r_{max}} dr 2 \pi r %s(r, \Delta t^*)$'
        % cor_name)

    ax_dt = plt.subplot(gs[2, 0])
    ax_dt.set_xscale(x_scale)
    ax_dt.set_yscale(dt_ys)
    ax_dt.set_ylabel(dt_label)

    ax_ratioTL = plt.subplot(gs[3, 0])
    ax_ratioTL.set_xscale(x_scale)
    ax_ratioTL.set_xlabel(r'$Pe$')
    ax_ratioTL.set_yscale(ratioTL_ys)
    ax_ratioTL.set_ylabel(r'$%s^T/%s^L(\Delta t^*)$' % (cor_name, cor_name))

    plt.setp(
        [ax.get_xticklabels() for ax in [ax_philoc, ax_chi, ax_dt]],
        visible=False)

    axes_legend = [plt.subplot(gs[:, 1 + traj]) for traj in range(comparisons)]

    for (phi, f_label, f_var, c_label, c_var, x_func, philoc_traj, chimsd_traj,
        clct_traj, colors_traj, markers_traj, linestyle, ax_legend)\
        in zip(densities, fix_label, fixed_var, var_label, var_c, pe_func,
        philoc, chimsd, clct, colors, markers, linestyles.values(),
        axes_legend):

        x_func_var_c = x_func(c_var)

        ax_philoc.axvline(x_func_var_c, color='black', linestyle=linestyle)
        for dir in philoc_traj.philocmax:
            var_value = philoc_traj.var_hash[dir]
            ax_philoc.scatter(
                x_func(var_value), philoc_traj.philocmax[dir],
                color=colors_traj[var_value],
                marker=markers_traj[philoc_traj.time_step[dir]])

        ax_chi.axvline(x_func_var_c, color='black', linestyle=linestyle)
        for dir in chimsd_traj.chimax:
            if not(chimsd_traj.isinvarinterval[dir]): continue
            var_value = chimsd_traj.var_hash[dir]
            ax_chi.scatter(
                x_func(var_value), chimsd_traj.chimax[dir],
                color=colors_traj[var_value],
                marker=markers_traj[chimsd_traj.time_step[dir]])

        ax_dt.axvline(x_func_var_c, color='black', linestyle=linestyle)
        for dir in chimsd_traj.dtmax:
            if not(chimsd_traj.isinvarinterval[dir]): continue
            var_value = chimsd_traj.var_hash[dir]
            ax_dt.scatter(
                x_func(var_value), chimsd_traj.dtmax[dir],
                color=colors_traj[var_value],
                marker=markers_traj[chimsd_traj.time_step[dir]])

        ax_ratioTL.axvline(x_func_var_c, color='black', linestyle=linestyle)
        for dir in clct_traj.ratioTL_max:
            var_value = clct_traj.var_hash[dir]
            ax_ratioTL.scatter(
                x_func(var_value), clct_traj.ratioTL_max[dir],
                color=colors_traj[var_value],
                marker=markers_traj[clct_traj.time_step[dir]])

        legend = [
            Line2D([0], [0], lw=0,
                label=r'$\phi=%1.2f, %s=%.1e$' % (phi, f_label, f_var)),
            Line2D([0], [0], lw=0),
            Line2D([0], [0], linestyle=linestyle, color='black',
                label=r'$%s = %.1e$' % (c_label, c_var)),
            Line2D([0], [0], lw=0)]
        legend += list(map(
            lambda var_value: Line2D([0], [0],
                color=colors_traj[var_value],
                label=r'$%s = %.1e$' % (c_label, var_value)),
            colors_traj))
        legend += [Line2D([0], [0], lw=0)]
        legend += list(map(
            lambda time_step: Line2D([0], [0], lw=0, color='black',
                marker=markers_traj[time_step],
                label=r'$dt = %.1e$' % time_step),
            markers_traj))
        ax_legend.axis('off')
        ax_legend.legend(handles=legend, loc='center')

    # SAVING

    fig.savefig(get_env('IMAGE_NAME', default=_image_name))

    # SHOW

    if get_env('SHOW', default=True, vartype=bool): plt.show()
