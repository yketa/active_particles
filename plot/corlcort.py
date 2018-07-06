"""
Module corlcort plots longitudinal, transverse and ratio of transverse and
longitudinal correlations, as well as single out values of these variables at
time of maximum cooperativity, for either different persistence times at fixed
self-propelling velocity or different self-propelling velocities at fixed
persistence time.

Simulation directories must follow the active_particles.naming.AHB2D naming
standard and input files in simulation directories must follow either the
active_particles.naming.Cuu standard (for displacement correlations), or
active_particles.naming.Cww standard (for displacement relative to centre of
mass displacement correlations), or active_particles.naming.Cdd standard (for
displacement norm correlations), or active_particles.naming.Cee (for direction
correlations).

Environment modes
-----------------
VARIABLE : string
    Plot of transversal and longitudinal correlations x-coordinate variable.
     _____________________________________________________________________
    | Mode    | Variable                    | x-coordinate if not(PECLET) |
    |_________|_____________________________|_____________________________|
    | 'dr'    | Rotation diffusion constant | \\tau = 1/dr                |
    |_________|_____________________________|_____________________________|
    | 'vzero' | self-propelling velocity    | vzero                       |
    |_________|_____________________________|_____________________________|
    DEFAULT: dr
PECLET : bool
    Plot transversal and longitudinal correlations at times of maximum
    cooperativity as functions of the Péclet number Pe = vzero/dr.
    DEFAULT: True
CORRELATION : string
    Correlations from which to extract longitudinal and transversal components.
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
FIT : bool
    Fit ratio of transversal and longitudinal correlations as logarithmic law
    of their x-coordinates for Péclet number lesser than the value
    corresponding to VAR_C.
    DEFAULT: False

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
DENSITY : float
    Packing fraction of particles.
    DEFAULT: active_particles.plot.corlcort._density
N : int
    Number of particles.
    DEFAULT: active_particles.plot.corlcort._N
VZERO ['dr' mode] : float
    Self-propulsion velocity.
    DEFAULT: active_particles.plot.corlcort._vzero
DR_MIN ['dr' mode] : float
    Minimum rotation diffusion constant.
    DEFAULT: active_particles.plot.corlcort._dr_min
DR_MAX ['dr' mode] : float
    Maximum rotation diffusion constant.
    DEFAULT: active_particles.plot.corlcort._dr_max
DR_C ['dr' and FIT mode] : float
    Transition rotation diffusion constant.
    DEFAULT: active_particles.plot.corlcort._dr_c
DR ['vzero' mode] : float
    Rotation diffusion constant.
    DEFAULT: active_particles.plot.corlcort._dr
VZERO_MIN ['vzero' mode] : float
    Minimum self-propulsion velocity.
    DEFAULT: active_particles.plot.corlcort._vzero_min
VZERO_MAX ['vzero' mode] : float
    Maximum self-propulsion velocity.
    DEFAULT: active_particles.plot.corlcort._vzero_max
VZERO_C ['vzero' and FIT mode] : float
    Transition self-propelling velocity.
    DEFAULT: active_particles.plot.corlcort._vzero_c
BOX_SIZE : float
	Size of the square box which was considered.
	DEFAULT: simulation box size
X_ZERO : float
	1st coordinate of the centre of the square box to consider.
	DEFAULT: 0
Y_ZERO : float
	2nd coordinate of the centre of the square box to consider.
	DEFAULT: 0
INITIAL_FRAME : int
    Frame to consider as initial for correlations.
    DEFAULT: active_particles.plot.corlcort._init_frame
INTERVAL_MAXIMUM : int
    Maximum number of intervals of length dt considered for correlations.
    DEFAULT: active_particles.plot.corlcort._int_max
N_CASES : int
    Number of boxes in each direction with which the displacement grid was
    computed.
    DEFAULT: active_particles.plot.corlcort._Ncases
R_MIN : float
    Minimum radius for correlations integration to get cooperativities.
    DEFAULT: active_particles.plot.corlcort._r_min
R_MAX : float
    Maximum radius for correlations integration to get cooperativities.
    DEFAULT: active_particles.plot.corlcort._r_max
FONT_SIZE : int
    Plot font size.
    DEFAULT: active_particles.plot.corlcort._font_size
MARKER_SIZE : int
    Plot marker size.
    DEFAULT: active_particles.plot.corlcort._marker_size
COLORMAP : string
    Plot colormap.
    DEFAULT: active_particles.plot.corlcort._colormap
RATIO_LEGEND : float
    Width ratio between legend and figure.
    DEFAULT: active_particles.plot.corlcort._ratio_legend
NCOL_LEGEND : int
    Number of columns for the legend.
    DEFAULT: active_particles.plot.corlcort._ncol_legend
WSPACE : float
    Plots width space.
    DEFAULT: active_particles.plot.corlcort._wspace
HSPACE : float
    Plots height space.
    DEFAULT: active_particles.plot.corlcort._hspace
COR_YS : string
    Correlations y-scale.
    DEFAULT: active_particles.plot.corlcort._cor_ys
DT_XS : string
    Lag time x-scale.
    DEFAULT: active_particles.plot.corlcort._dt_xs
X_SCALE : string
    Plot variable x-scale.
    DEFAULT: active_particles.plot.corlcort._x_scale
"""

import active_particles.naming as naming

from active_particles.init import get_env, dir_list, isnumber

from os import environ as envvar
if __name__ == '__main__': envvar['SHOW'] = 'True'
from os.path import join as joinpath

from active_particles.plot.plot import list_colormap, list_markers,\
    list_linestyles
from active_particles.plot.chi_msd import ChiMsd

from collections import OrderedDict

import pickle

import numpy as np

from scipy import stats as st

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

# DEFAULT VARIABLES

_dr = 3e-4      # default rotation diffusion constant

_dr_min = 1e-5  # default minimum diffusion rotation constant
_dr_max = 1e-2  # default maximum diffusion rotation constant
_dr_c = 3e-4    # default transition diffusion rotation constant

_vzero = 1e-2   # default self-propelling velocity

_vzero_min = 1e-2   # default minimum self-propelling velocity
_vzero_max = 1e-1   # default maximum self-propelling velocity
_vzero_c = 5e-1     # default transition self-propelling velocity

_density = 0.8  # default packing fraction of particles
_N = int(1e5)   # default number of particles

_init_frame = 0 # default frame to consider as initial for correlations
_int_max = 1    # default maximum number of intervals of length dt considered for correlations
_Ncases = 500   # default number of boxes in each direction with which the displacement grid is computed

_r_min = 1  # default minimum radius for correlations integration
_r_max = 20 # default maximum radius for correlations integration

_font_size = 15         # default font size for the plot
_marker_size = 5        # default plot marker size
_colormap = 'jet'       # default plot colormap

_ncol_legend = 2    # default number of legend columns
_ratio_legend = 10  # default width ratio between graph and legend

_wspace = 0.4   # default plots width space
_hspace = 0.4   # default plots height space

_cor_ys = 'linear'  # default correlations y-scale
_dt_xs = 'log'      # default lag time x-scale
_x_scale = 'log'    # default plot variable x-scale

# FUNCTIONS AND CLASSES

class CLCT:
    """
    Search and read correlation files, extract longtidunal and transversal
    correlations.
    """

    def __init__(self, data_dir, dir_standard, dir_attributes, parameters_file,
        var, var_min, var_max, excluded_dir=''):
        """
        Create list of directories to consider and compute plot variable values
        associated to them.

        Parameters
        ----------
        data_dir : string
            Data directory.
        dir_standard : active_particles.naming._File standard
            Simulation directory naming object.
        dir_attributes : hash table
            Attributes to be displayed in directory names.
        parameters_file : string
            Simulations parameters file name.
        var : string
            Plot variable name.
        var_min : float
            Minimum plot variable value.
        var_max : float
            Maximum plot variable value.
        excluded_dir : string
            Names of directories to be ignored. (default: '')
        """

        self.data_dir = data_dir
        self.dir_standard = dir_standard
        self.dir_attributes = dir_attributes
        self.excluded_dir = excluded_dir

        self.parameters_file = parameters_file

        self.var = var
        self.var_min = var_min
        self.var_max = var_max

        self.dirs, self.var_hash, self.var_list, _, _ = dir_list(
            self.data_dir, self.dir_standard, self.dir_attributes,
            self.var, self.var_min, self.var_max,
            self.parameters_file, excluded_dir=self.excluded_dir,
            include_out=False)

    def calculate(self, cor_standard, cor_attributes, multiply_with_dr=True):
        """
        Extract longtidunal and transversal correlations.

        Parameters
        ----------
        cor_standard : active_particles.naming._File standard
            Correlation files naming object.
        cor_attributes : hash table
            Attributes to be displayed in correlation file names.
        multiply_with_dr : bool
            Consider the product of rotation diffusion constant and lag time
            rather than just lag time. (default: True)
        """

        self.cor_standard = cor_standard
        self.cor_attributes = cor_attributes

        self.multiply_with_dr = multiply_with_dr

        self.time_step = {} # hash table of directories' simulation time step
        self.corL = {}      # hash table of longitudinal correlation
        self.corT = {}      # hash table of transversal correlation
        self.ratioTL = {}   # hash table of ratio of transversal and longitudinal correlations

        for dir in self.dirs:

            with open(
                joinpath(self.data_dir, dir, self.parameters_file), 'rb')\
                as param_file:
                parameters = pickle.load(param_file)    # simulation parameters hash table

            pdts = parameters['period_dump']*parameters['time_step']    # time corresponding to one dump length of time
            if self.multiply_with_dr: pdts *= parameters['dr']          # plot dr*dt rather than dt

            corL_dir, corT_dir, ratioTL_dir = [], [], []    # list of longitudinal, transversal and ration of transversal and longitudinal correlations for current directory
            for cor_filename in self.cor_standard.get_files(
                directory=joinpath(self.data_dir, dir), **self.cor_attributes): # loop over correlations files in directory

                with open(joinpath(self.data_dir, dir, cor_filename), 'rb')\
                    as cor_file:
                    corL, corT = pickle.load(cor_file)[-2:]
                    if not(isnumber(corL)) or not(isnumber(corT)): continue
                dt = pdts*self.cor_standard.get_data(cor_filename, 'dt')[0]
                corL_dir += [[dt, corL]]
                corT_dir += [[dt, corT]]
                ratioTL_dir += [[dt, corT/corL]]

            if corL_dir and corT_dir and ratioTL_dir:
                self.corL[dir] = np.transpose(
                    sorted(corL_dir, key=lambda el: el[0]))
                self.corT[dir] = np.transpose(
                    sorted(corT_dir, key=lambda el: el[0]))
                self.ratioTL[dir] = np.transpose(
                    sorted(ratioTL_dir, key=lambda el: el[0]))
                self.time_step[dir] = parameters['time_step']

        self.time_step_list = sorted(OrderedDict.fromkeys(
            self.time_step.values()))   # list of time steps

    def calculate_max(self, dtmax):
        """
        Extract longitudinal, transversal and ratio of transversal and
        longitudinal correlations at times of maximum cooperativity.

        Parameters
        ----------
        dtmax : hash table
            Hash table of time of maximum cooperativities with directory names
            as keys.
        """

        self.corL_max = {}      # hash table of longitudinal correlations at time of maximum cooperativity
        self.corT_max = {}      # hash table of transversal correlations at time of maximum cooperativity
        self.ratioTL_max = {}   # hash table of ratio of transversal and longitudinal correlations at time of maximum cooperativity

        for dir in dtmax:

            if dir in self.corL:
                self.corL_max[dir] = self.corL[dir][1,
                    np.argmin(np.abs(self.corL[dir][0] - dtmax[dir]))]

            if dir in self.corT:
                self.corT_max[dir] = self.corT[dir][1,
                    np.argmin(np.abs(self.corT[dir][0] - dtmax[dir]))]

            if dir in self.ratioTL:
                self.ratioTL_max[dir] = self.ratioTL[dir][1,
                    np.argmin(np.abs(self.ratioTL[dir][0] - dtmax[dir]))]

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLES DEFINITIONS

    mode = get_env('VARIABLE', default='dr')                # plotting variable
    peclet = get_env('PECLET', default=True, vartype=bool)  # display Péclet number rather than mode variable

    if mode == 'dr':

        vzero = get_env('VZERO', default=_vzero, vartype=float) # self-propelling velocity
        attributes = {'vzero': vzero}                           # attributes displayed in filenames

        var = 'dr'                                                  # plot variable
        var_min = get_env('DR_MIN', default=_dr_min, vartype=float) # minimum rotation diffusion constant
        var_max = get_env('DR_MAX', default=_dr_max, vartype=float) # maximum rotation diffusion constant
        var_c = get_env('DR_C', default=_dr_c, vartype=float)       # transition rotation diffusion constant

        var_label = r'$\tilde{\nu}_r$'  # variable label

        if peclet:
            x_func = lambda x: vzero/x                      # x-coordinate as function of plot variable
        else:
            x_label = r'$\tau_r \equiv \tilde{\nu}_r^{-1}$' # x-coordinate label
            x_func = lambda x: 1/x                          # x-coordinate as function of plot variable

    elif mode == 'vzero':

        dr = get_env('DR', default=_dr, vartype=float)  # rotation diffusion constant
        attributes = {'dr': dr}                         # attributes displayed in filenames

        var = 'vzero'                                                       # plot variable
        var_min = get_env('VZERO_MIN', default=_vzero_min, vartype=float)   # minimum self-propelling velocity
        var_max = get_env('VZERO_MAX', default=_vzero_max, vartype=float)   # maximum self-propelling velocity
        var_c = get_env('VZERO_C', default=_vzero_c, vartype=float)         # transition self-propelling velocity

        var_label = r'$\tilde{v}$'  # variable label

        if peclet:
            x_func = lambda x: x/dr     # x-coordinate as function of plot variable
        else:
            x_label = r'$\tilde{v}$'    # x-coordinate label
            x_func = lambda x: x        # x-coordinate as function of plot variable

    else: raise ValueError('Mode %s is not known.' % mode)  # mode is not known

    if peclet: x_label = r'$Pe$'    # x-coordinate label

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

    density = get_env('DENSITY', default=_density, vartype=float)   # packing fraction of particles
    N = get_env('N', default=_N, vartype=int)                       # number of particles

    box_size = get_env('BOX_SIZE', vartype=float)     # size of the square box which was considered
    centre = (get_env('X_ZERO', default=0, vartype=float),
		get_env('Y_ZERO', default=0, vartype=float))  # centre of the box

    init_frame = get_env('INITIAL_FRAME', default=_init_frame, vartype=int) # frame to consider as initial for correlations
    int_max = get_env('INTERVAL_MAXIMUM', default=_int_max, vartype=int)    # maximum number of intervals of length dt considered for correlations
    Ncases = get_env('N_CASES', default=_Ncases, vartype=int)               # number of boxes in each direction with which the displacement grid is computed

    r_min = get_env('R_MIN', default=_r_min, vartype=float) # minimum radius for correlations integration
    r_max = get_env('R_MAX', default=_r_max, vartype=float) # maximum radius for correlations integration

    # NAMING

    attributes = {**attributes, 'density': density, 'N': N,
        'init_frame': init_frame, 'int_max': int_max, 'Ncases': Ncases,
        'box_size': box_size, 'x_zero': centre[0], 'y_zero': centre[1]} # attributes displayed in filenames
    naming_simdir = naming.AHB2D()                                      # simulation directory naming object

    # PLOT PARAMETERS

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=float)     # plot font size
    marker_size = get_env('MARKER_SIZE', default=_marker_size, vartype=int) # plot marker size
    mpl.rcParams.update({'font.size': font_size,
        'lines.markersize': marker_size})

    colormap = get_env('COLORMAP', default=_colormap)       # plot colormap

    ratio_legend = get_env('RATIO_LEGEND',
        default=_ratio_legend, vartype=float)                               # width ratio between graph and legend
    ncol_legend = get_env('NCOL_LEGEND', default=_ncol_legend, vartype=int) # number of legend columns

    wspace = get_env('WSPACE', default=_wspace, vartype=float)  # plots width space
    hspace = get_env('HSPACE', default=_hspace, vartype=float)  # plots height space

    cor_ys = get_env('COR_YS', default=_cor_ys)     # correlations y-scale
    dt_xs = get_env('DT_XS', default=_dt_xs)        # lag time x-scale
    x_scale = get_env('X_SCALE', default=_x_scale)  # plot variable x-scale

    multiply_with_dr = get_env('DRDT', default=True, vartype=bool)  # plot dr*dt rather than dt
    if multiply_with_dr:
        dt_label = '\\tilde{\\nu}_r \\Delta t'                      # dt label
    else: dt_label = '\\Delta t'                                    # dt label

    fit_max = get_env('FIT', default=False, vartype=bool)   # fit ratio of transversal and longitudinal correlations as logarithmic law of their x-coordinates for Péclet number lesser than the value corresponding to variable transition value

    # CALCULATION

    chimsd = ChiMsd(data_dir, naming_simdir, attributes, parameters_file,
        var, var_min, var_max, excluded_dir=excluded_directories)   # cooperativities calculator
    chimsd.calculate(naming_cor, attributes, r_min, r_max, box_size=box_size,
        multiply_with_dr=multiply_with_dr)                          # calculate cooperativities

    clct = CLCT(data_dir, naming_simdir, attributes, parameters_file,
        var, var_min, var_max, excluded_dir=excluded_directories)               # longtidunal and transversal correlations calculator
    clct.calculate(naming_cor, attributes, multiply_with_dr=multiply_with_dr)   # calculate longtidunal and transversal correlations
    clct.calculate_max(chimsd.dtmax)                                            # calculate longitudinal, transversal and ratio of transversal and longitudinal correlations at times of maximum cooperativity

    # PLOT

    colors = list_colormap(clct.var_list, colormap=colormap)    # plot colors hash table
    markers = list_markers(clct.time_step_list)                 # plot markers hash table

    # RATIO

    fig_ratio = plt.figure()
    fig_ratio.subplots_adjust(wspace=wspace, hspace=hspace)
    suptitle = (
        r'$\phi=%1.2f, N=%.2e$' % (density, N)
        + (r'$, \tilde{v}=%.2e$' % vzero if mode == 'dr' else
        r'$, \tilde{\nu}_r = %.2e$' % dr)
        + '\n'
        + r'$S_{init}=%.2e, S_{max}=%.2e$' % (init_frame, Ncases)
        + r'$, N_{cases}=%.2e$' % int_max
        + r'$, r_{min}=%.2e, r_{max}=%.2e$' % (r_min, r_max))
    fig_ratio.suptitle(suptitle)

    gs = GridSpec(2, 2, width_ratios=[1, 1/ratio_legend])

    ax_ratio = plt.subplot(gs[0, 0])
    ax_ratio.set_xlabel(r'$%s$' % dt_label)
    ax_ratio.set_xscale(dt_xs)
    ax_ratio.set_ylabel(r'$%s^T/%s^L$' % (cor_name, cor_name))
    ax_ratio.set_yscale(cor_ys)
    lines_dt = [Line2D([0], [0], lw=0, marker='s', color='black',
        label=r'$%s^*$' % dt_label)]
    ax_ratio.legend(handles=lines_dt)

    ax_ratio_max = plt.subplot(gs[1, 0])
    ax_ratio_max.set_xlabel(x_label)
    ax_ratio_max.set_xscale(x_scale)
    ax_ratio_max.set_ylabel(r'$%s^T/%s^L(\Delta t^*)$' % (cor_name, cor_name))
    ax_ratio_max.set_yscale(cor_ys)
    ax_ratio_max.axhline(1, linestyle='--', color='black')
    lines_time_step = list(map(
        lambda time_step: Line2D([0], [0], lw=0, marker=markers[time_step],
            color='black', label=r'$dt=%.0e$' % time_step),
        clct.time_step_list))
    lines_ax_ratio_max = lines_time_step + [Line2D([0], [0], lw=0, label=''),
        Line2D([0], [0], linestyle='--', color='black', label='1')]

    for dir in clct.ratioTL:
        ax_ratio.plot(*clct.ratioTL[dir], color=colors[clct.var_hash[dir]])

    for dir in clct.ratioTL_max:
        var_value = clct.var_hash[dir]
        time_step_value = clct.time_step[dir]
        ax_ratio.scatter(chimsd.dtmax[dir], clct.ratioTL_max[dir],
            marker='s', color=colors[var_value])
        ax_ratio_max.scatter(x_func(var_value), clct.ratioTL_max[dir],
            marker=markers[time_step_value], color=colors[clct.var_hash[dir]])

    if fit_max: # fit ratio of transversal and longitudinal correlations as logarithmic law of their x-coordinates

        # VALUES

        x_func_var_c = x_func(var_c)
        ratioTL_max_lt_pec = []
        for dir in clct.ratioTL_max:

            x_func_var_value = x_func(clct.var_hash[dir])
            if x_func_var_value < x_func_var_c:
                ratioTL_max_lt_pec += [[np.log(x_func_var_value),
                    clct.ratioTL_max[dir]]]

        ratioTL_max_lt_pec = np.transpose(ratioTL_max_lt_pec)

        # FIT

        slope, intercept, _, _, stderr = st.linregress(*ratioTL_max_lt_pec)

        # SORT FOR PLOTS

        ratioTL_max_lt_pec[0].sort()

        # PLOT

        ax_ratio_max.plot(np.exp(ratioTL_max_lt_pec[0]),
            slope*ratioTL_max_lt_pec[0] + intercept,
            color='black', linestyle='-.')
        lines_ax_ratio_max += [Line2D([0], [0], linestyle='-.', color='black',
            label='slope = ' + r'$%.1e\pm%.0e$' % (slope, stderr))]

    ax_ratio_max.legend(handles=lines_ax_ratio_max)

    ax_legend_ratio = plt.subplot(gs[:, 1])
    ax_legend_ratio.axis('off')
    lines_var = list(map(
        lambda var_value: Line2D([0], [0], color=colors[var_value], lw=2,
            label='%s = %.0e' % (var_label, var_value)),
        clct.var_list))
    ax_legend_ratio.legend(handles=lines_var, loc='center', ncol=ncol_legend)

    # SEPARATE

    fig_separate = plt.figure()
    fig_separate.subplots_adjust(wspace=wspace, hspace=hspace)
    fig_separate.suptitle(suptitle)

    gs = GridSpec(2, 3, width_ratios=[1, 1, 2/ratio_legend])

    ax_corL = plt.subplot(gs[0, 0])
    ax_corL.set_xlabel(r'$%s$' % dt_label)
    ax_corL.set_xscale(dt_xs)
    ax_corL.set_ylabel(r'$%s^L$' % cor_name)
    ax_corL.set_yscale(cor_ys)
    ax_corL.legend(handles=lines_dt)

    ax_corL_max = plt.subplot(gs[0, 1])
    ax_corL_max.set_xlabel(x_label)
    ax_corL_max.set_xscale(x_scale)
    ax_corL_max.set_ylabel(r'$%s^L(\Delta t^*)$' % cor_name)
    ax_corL_max.set_yscale(cor_ys)
    ax_corL_max.legend(handles=lines_time_step)

    for dir in clct.corL:
        ax_corL.plot(*clct.corL[dir], color=colors[clct.var_hash[dir]])

    for dir in clct.corL_max:
        var_value = clct.var_hash[dir]
        time_step_value = clct.time_step[dir]
        ax_corL.scatter(chimsd.dtmax[dir], clct.corL_max[dir],
            marker='s', color=colors[var_value])
        ax_corL_max.scatter(x_func(var_value), clct.corL_max[dir],
            marker=markers[time_step_value], color=colors[clct.var_hash[dir]])

    ax_corT = plt.subplot(gs[1, 0])
    ax_corT.set_xlabel(r'$%s$' % dt_label)
    ax_corT.set_xscale(dt_xs)
    ax_corT.set_ylabel(r'$%s^T$' % cor_name)
    ax_corT.set_yscale(cor_ys)
    ax_corT.legend(handles=lines_dt)

    ax_corT_max = plt.subplot(gs[1, 1])
    ax_corT_max.set_xlabel(x_label)
    ax_corT_max.set_xscale(x_scale)
    ax_corT_max.set_ylabel(r'$%s^T(\Delta t^*)$' % cor_name)
    ax_corT_max.set_yscale(cor_ys)
    ax_corT_max.legend(handles=lines_time_step)

    for dir in clct.corT:
        ax_corT.plot(*clct.corT[dir], color=colors[clct.var_hash[dir]])

    for dir in clct.corT_max:
        var_value = clct.var_hash[dir]
        time_step_value = clct.time_step[dir]
        ax_corT.scatter(chimsd.dtmax[dir], clct.corT_max[dir],
            marker='s', color=colors[var_value])
        ax_corT_max.scatter(x_func(var_value), clct.corT_max[dir],
            marker=markers[time_step_value], color=colors[clct.var_hash[dir]])

    ax_legend_separate = plt.subplot(gs[:, 2])
    ax_legend_separate.axis('off')
    ax_legend_separate.legend(handles=lines_var,
        loc='center', ncol=ncol_legend)

    # SHOW

    plt.show()
