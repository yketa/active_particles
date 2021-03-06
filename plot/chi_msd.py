"""
Module chi_msd plots cooperativities, maximum cooperativities, times of
maximum cooperativities and mean square displacements for either different
persistence times at fixed self-propelling velocity or different
self-propelling velocities at fixed persistence time.

Simulation directories must follow the active_particles.naming.AHB2D naming
standard and input files in simulation directories must follow either the
active_particles.naming.Cuu standard (for cooperativities from displacement
correlations), or active_particles.naming.Cww standard (for cooperativities
from displacement relative to centre of mass displacement correlations), or
active_particles.naming.Cdd standard (for cooperativities from displacement
norm correlations), or active_particles.naming.Cee (for cooperativities from
displacement direction correlations), and the active_particles.naming.Msd (mean
square displacements).

Environment modes
-----------------
VARIABLE : string
    Plot of maximum cooperativities and times of maximum cooperativity
    x-coordinate variable.
     _____________________________________________________________________
    | Mode    | Variable                    | x-coordinate if not(PECLET) |
    |_________|_____________________________|_____________________________|
    | 'dr'    | Rotation diffusion constant | \\tau = 1/dr                |
    |_________|_____________________________|_____________________________|
    | 'vzero' | self-propelling velocity    | vzero                       |
    |_________|_____________________________|_____________________________|
    DEFAULT: dr
PECLET : bool
    Plot maximum cooperativities and times of maximum cooperativity as
    functions of the Péclet number Pe = vzero/dr.
    DEFAULT: True
CORRELATION : string
    Correlations from which to calculate cooperativities.
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
MSDDT : bool
    Divide mean square displacement by lag time.
    DEFAULT: True
FIT : bool
    Fit maximum cooperativity and time of maximum cooperativity as power law of
    their x-coordinates, on both sides of the variabe transition value VAR_C.
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
    DEFAULT: active_particles.plot.chi_msd._density
N : int
    Number of particles.
    DEFAULT: active_particles.plot.chi_msd._N
VZERO ['dr' mode] : float
    Self-propulsion velocity.
    DEFAULT: active_particles.plot.chi_msd._vzero
DR_MIN ['dr' mode] : float
    Minimum rotation diffusion constant.
    NOTE: Rotation diffusion constants lesser than DR_MIN will appear in the
    mean square displacement figure but not in the maximum cooperativity
    figure.
    DEFAULT: active_particles.plot.chi_msd._dr_min
DR_MAX ['dr' mode] : float
    Maximum rotation diffusion constant.
    NOTE: Rotation diffusion constants greater than DR_MAX will appear in the
    mean square displacement figure but not in the maximum cooperativity
    figure.
    DEFAULT: active_particles.plot.chi_msd._dr_max
DR_C ['dr' and FIT mode] : float
    Transition rotation diffusion constant.
    DEFAULT: active_particles.plot.chi_msd._dr_c
DR ['vzero' mode] : float
    Rotation diffusion constant.
    DEFAULT: active_particles.plot.chi_msd._dr
VZERO_MIN ['vzero' mode] : float
    Minimum self-propulsion velocity.
    NOTE: Self-propelling velocities lesser than VZERO_MIN will appear in the
    mean square displacement figure but not in the maximum cooperativity
    figure.
    DEFAULT: active_particles.plot.chi_msd._vzero_min
VZERO_MAX ['vzero' mode] : float
    Maximum self-propulsion velocity.
    NOTE: Self-propelling velocities greater than VZERO_MIN will appear in the
    mean square displacement figure but not in the maximum cooperativity
    figure.
    DEFAULT: active_particles.plot.chi_msd._vzero_max
VZERO_C ['vzero' and FIT mode] : float
    Transition self-propelling velocity.
    DEFAULT: active_particles.plot.chi_msd._vzero_c
BOX_SIZE : float
	Size of the square box which was considered.
	DEFAULT: simulation box size
X_ZERO : float
	1st coordinate of the centre of the square box to consider.
	DEFAULT: 0
Y_ZERO : float
	2nd coordinate of the centre of the square box to consider.
	DEFAULT: 0
INITIAL_FRAME_COR : int
    Frame to consider as initial for correlations.
    DEFAULT: active_particles.plot.chi_msd._init_frame_cor
INTERVAL_MAXIMUM_COR : int
    Maximum number of intervals of length dt considered for correlations.
    DEFAULT: active_particles.plot.chi_msd._int_max_cor
N_CASES : int
    Number of boxes in each direction with which the displacement grid was
    computed.
    DEFAULT: active_particles.plot.chi_msd._Ncases_cor
INITIAL_FRAME_MSD : int separated by ':'
    Display only mean square displacements calculated with initial frame from
    this list except if this list is empty.
    DEFAULT: []
INTERVAL_MAXIMUM_MSD : int
    Maximum number of intervals of length dt considered for mean square
    displacements calculation.
    DEFAULT: active_particles.plot.chi_msd._int_max_msd
INTERVAL_PERIOD : int
    Period of dumps for which mean square displacements were calculated.
    DEFAULT: active_particles.plot.chi_msd._int_period_msd
R_MIN : float
    Minimum radius for correlations integration.
    DEFAULT: active_particles.plot.chi_msd._r_min
R_MAX : float
    Maximum radius for correlations integration.
    DEFAULT: active_particles.plot.chi_msd._r_max
FONT_SIZE : int
    Plot font size.
    DEFAULT: active_particles.plot.chi_msd._font_size
MARKER_SIZE : int
    Plot marker size.
    DEFAULT: active_particles.plot.chi_msd._marker_size
COLORMAP : string
    Plot colormap.
    DEFAULT: active_particles.plot.chi_msd._colormap
COLORMAP0 : string
    Plot colormap for variables out of variable window.
    DEFAULT: active_particles.plot.chi_msd._colormap0
RATIO_LEGEND : float
    Width ratio between legend and figure.
    DEFAULT: active_particles.plot.chi_msd._ratio_legend
NCOL_LEGEND : int
    Number of columns for the legend.
    DEFAULT: active_particles.plot.chi_msd._ncol_legend
WSPACE : float
    Plots width space.
    DEFAULT: active_particles.plot.chi_msd._wspace
HSPACE : float
    Plots height space.
    DEFAULT: active_particles.plot.chi_msd._hspace
CHI_XS : string
    Cooperativity plot x-scale.
    DEFAULT: active_particles.plot.chi_msd._chi_xs
CHI_YS : string
    Cooperativity plot y-scale.
    DEFAULT: active_particles.plot.chi_msd._chi_ys
X_SCALE : string
    Maximum cooperativity and time of maximum cooperativity plots x-scale.
    DEFAULT: active_particles.plot.chi_msd._x_scale
DTMAX_YS : string
    Time of maximum cooperativity plot y-scale.
    DEFAULT: active_particles.plot.chi_msd._dtmax_ys
CHIMAX_YS : string
    Maximum cooperativity plot y-scale.
    DEFAULT: active_particles.plot.chi_msd._chimax_ys
MSD_XS : string
    Mean square displacement plot x-scale.
    DEFAULT: active_particles.plot.chi_msd._msd_xs
MSD_YS : string
    Mean square displacement plot y-scale.
    DEFAULT: active_particles.plot.chi_msd._msd_ys
"""

import active_particles.naming as naming

from active_particles.init import get_env, get_env_list, dir_list

from os import environ as envvar
if __name__ == '__main__': envvar['SHOW'] = 'True'
from os.path import join as joinpath

from active_particles.plot.plot import list_colormap, list_markers,\
    list_linestyles
from active_particles.analysis.cuu import c1Dtochi

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

_init_frame_cor = 0 # default frame to consider as initial for correlations
_int_max_cor = 1    # default maximum number of intervals of length dt considered for correlations
_Ncases_cor = 500   # default number of boxes in each direction with which the displacement grid is computed

_int_max_msd = 1    # default maximum number of intervals of length dt considered for correlations
_int_period_msd = 1 # default period of dumps for which mean square displacements were calculated

_r_min = 1  # default minimum radius for correlations integration
_r_max = 20 # default maximum radius for correlations integration

_font_size = 15         # default font size for the plot
_marker_size = 5        # default plot marker size
_colormap = 'jet'       # default plot colormap
_colormap0 = 'Greys'    # default plot colormap for variables out of variable window

_ncol_legend = 2    # default number of legend columns
_ratio_legend = 10  # default width ratio between graph and legend

_wspace = 0.4   # default plots width space
_hspace = 0.05  # default plots height space

_chi_xs = 'log'     # default cooperativity plot x-scale
_chi_ys = 'log'     # default cooperativity plot y-scale
_x_scale = 'log'    # default x-scale
_dtmax_ys = 'log'   # default time of maximum cooperativity plot y-scale
_chimax_ys = 'log'  # default maximum cooperativity plot y-scale
_msd_xs = 'log'     # default mean square displacement plot x-scale
_msd_ys = 'log'     # default mean square displacement plot y-scale

# FUNCTIONS AND CLASSES

class ChiMsd:
    """
    Search and read correlations files, calculate cooperativities.
    Also search and read mean square displacement files.
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

        (self.dirs, self.var_hash, self.var_list, self.var0_list,
            self.isinvarinterval) = dir_list(
                self.data_dir, self.dir_standard, self.dir_attributes,
                self.var, self.var_min, self.var_max,
                self.parameters_file, excluded_dir=self.excluded_dir,
                include_out=True)

        self.calculate_msd = False  # calculate mean square displacements with cooperativities

    def calculate(self, cor_standard, cor_attributes, r_min, r_max,
        box_size=None, multiply_with_dr=True):
        """
        Calculate cooperativities, maximum cooperativities and times of
        maximum cooperativites.
        Also calculates mean square displacements.

        Parameters
        ----------
        cor_standard : active_particles.naming._File standard
            Correlation files naming object.
        cor_attributes : hash table
            Attributes to be displayed in correlation file names.
        r_min : float
            Minimum radius for correlations integration.
        r_max : float
            Maximum radius for correlations integration.
        box_size : float or None
            Size of the square box which was considered. (default: None)
            NOTE: if None, then the size is taken as the simulation box size.
        multiply_with_dr : bool
            Consider the product of rotation diffusion constant and lag time
            rather than just lag time. (default: True)
        """

        self.cor_standard = cor_standard
        self.cor_attributes = cor_attributes

        self.r_min = r_min
        self.r_max = r_max

        self.box_size = box_size
        self.multiply_with_dr = multiply_with_dr

        self.time_step = {}                     # hash table of directories' simulation time step
        self.chi = {}                           # hash table of list of lag times and corresponding cooperativities
        self.dtmax = {}                         # hash table of time of maximum cooperativities
        self.chimax = {}                        # hash table of maximum cooperativities
        self.islocalmax = {}                    # hash table of booleans indicating if the time of maximum cooperativity is a local maximum

        for dir in self.dirs:

            # COOPERATIVITY

            with open(
                joinpath(self.data_dir, dir, self.parameters_file), 'rb')\
                as param_file:
                parameters = pickle.load(param_file)    # simulation parameters hash table

            if self.box_size == None: L = parameters['box_size']
            else: L = self.box_size # box size

            pdts = parameters['period_dump']*parameters['time_step']    # time corresponding to one dump length of time
            if self.multiply_with_dr: pdts *= parameters['dr']          # plot dr*dt rather than dt

            chidir = []                                                         # list of lag times and corresponding cooperativities for current directory
            for cor_filename in self.cor_standard.get_files(
                directory=joinpath(self.data_dir, dir), **self.cor_attributes): # loop over correlations files in directory

                with open(joinpath(self.data_dir, dir, cor_filename), 'rb')\
                    as cor_file:
                    c1D = pickle.load(cor_file)[1]                          # 1D correlation
                chidir += [[
                    pdts*self.cor_standard.get_data(cor_filename, 'dt')[0],
                    c1Dtochi(c1D, L, r_min=self.r_min, r_max=self.r_max)]]  # cooperativity

            if not(chidir): continue    # no cooperativities files with requested attributes

            self.time_step[dir] = parameters['time_step']       # simulation time step
            self.dtmax[dir], self.chimax[dir] = max(chidir,
                key=lambda el: el[1])                           # time of maximum cooperativity and maximum cooperativity
            self.chi[dir] = np.transpose(sorted(chidir,
                key=lambda el: el[0]))                          # cooperativity
            self.islocalmax[dir] =\
                self.dtmax[dir] > min(self.chi[dir][:, 0])\
                and self.dtmax[dir] < max(self.chi[dir][:, 0])  # is dtmax a local maximum

            if not(self.calculate_msd): continue    # do not calculate mean square displacements

            # MEAN SQUARE DISPLACEMENT

            for msd_filename in self.msd_standard.get_files(
                directory=joinpath(self.data_dir, dir), **self.msd_attributes): # loop over mean square displacements files in directory
                init_frame = self.msd_standard.get_data(
                    msd_filename, 'init_frame')[0]                              # initial frame

                if not(self.init_frame_msd)\
                    or init_frame in self.init_frame_msd:

                    self.msd[(dir, init_frame)] = np.genfromtxt(
                        fname=joinpath(self.data_dir, dir, msd_filename),
                        delimiter=',', skip_header=True)    # mean square displacement

                    if self.divide_by_dt:
                        self.msd[(dir, init_frame)][:, 1] /=\
                            self.msd[(dir, init_frame)][:, 0]
                        self.msd[(dir, init_frame)][:, 2] /=\
                            self.msd[(dir, init_frame)][:, 0]

                    if self.multiply_with_dr:
                        self.msd[(dir, init_frame)][:, 0] *= parameters['dr']

        self.time_step_list = sorted(OrderedDict.fromkeys(
            self.time_step.values()))   # list of time steps

    def calculate_with_msd(self, cor_standard, cor_attributes, r_min, r_max,
        msd_standard, msd_attributes, init_frame_msd, box_size=None,
        multiply_with_dr=True, divide_by_dt=True):
        """
        Calculate cooperativities, maximum cooperativities and times of
        maximum cooperativites.
        Also calculates mean square displacements.

        Parameters
        ----------
        cor_standard : active_particles.naming._File standard
            Correlation files naming object.
        cor_attributes : hash table
            Attributes to be displayed in correlation file names.
        r_min : float
            Minimum radius for correlations integration.
        r_max : float
            Maximum radius for correlations integration.
        msd_standard : active_particles.naming._File standard
            Mean square displacement files naming object.
        msd_attributes : hash table
            Attributes to be displayed in mean square displacement file names.
        init_frame_msd : list of int
            Calculate only mean square displacements calculated with initial
            frame from this list except if this list is empty.
        box_size : float or None
            Size of the square box which was considered. (default: None)
            NOTE: if None, then the size is taken as the simulation box size.
        multiply_with_dr : bool
            Consider the product of rotation diffusion constant and lag time
            rather than just lag time. (default: True)
        divide_by_dt : bool
            Divide mean square displacement by lag time. (default: True)
        """

        self.calculate_msd = True

        self.msd_standard = msd_standard
        self.msd_attributes = msd_attributes
        self.init_frame_msd = init_frame_msd

        self.divide_by_dt = divide_by_dt

        self.msd = {}   # hash table of lists of lag times and corresponding mean square displacement and associated standard error

        self.calculate(cor_standard, cor_attributes, r_min, r_max,
            box_size=box_size, multiply_with_dr=multiply_with_dr)

        self.init_frame_list = sorted(OrderedDict.fromkeys(
            [init_frame for dir, init_frame in self.msd]))  # list of mean square displacements initial frames

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

    init_frame_cor = get_env('INITIAL_FRAME_COR',
        default=_init_frame_cor, vartype=int)                           # frame to consider as initial for correlations
    int_max_cor = get_env('INTERVAL_MAXIMUM_COR',
        default=_int_max_cor, vartype=int)                              # maximum number of intervals of length dt considered for correlations
    Ncases_cor = get_env('N_CASES', default=_Ncases_cor, vartype=int)   # number of boxes in each direction with which the displacement grid is computed

    init_frame_msd = get_env_list('INITIAL_FRAME_MSD',
        delimiter=':', vartype=int)             # display only mean square displacements calculated with initial frame from this list except if this list is empty
    int_max_msd = get_env('INTERVAL_MAXIMUM_MSD',
        default=_int_max_msd, vartype=int)      # maximum number of intervals of length dt considered for mean square displacements calculation
    int_period_msd = get_env('INTERVAL_PERIOD',
        default=_int_period_msd, vartype=int)   # period of dumps for which mean square displacements were calculated

    r_min = get_env('R_MIN', default=_r_min, vartype=float) # minimum radius for correlations integration
    r_max = get_env('R_MAX', default=_r_max, vartype=float) # maximum radius for correlations integration

    # NAMING

    common_attributes = {**attributes, 'density': density, 'N': N}  # attributes to be displayed in file names
    attributes_cor = {**common_attributes, 'init_frame': init_frame_cor,
        'int_max': int_max_cor, 'Ncases': Ncases_cor, 'box_size': box_size,
        'x_zero': centre[0], 'y_zero': centre[1]}                   # attributes displayed in file names specifically for correlations
    attributes_msd = {**common_attributes, 'int_max': int_max_msd,
        'int_period': int_period_msd}                               # attributes displayed in file names specifically for mean square displacements
    naming_msd = naming.Msd()                                       # mean square displacements naming object
    naming_simdir = naming.AHB2D()                                  # simulation directory naming object

    # PLOT PARAMETERS

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=float)     # plot font size
    marker_size = get_env('MARKER_SIZE', default=_marker_size, vartype=int) # plot marker size
    mpl.rcParams.update({'font.size': font_size,
        'lines.markersize': marker_size})

    colormap = get_env('COLORMAP', default=_colormap)       # plot colormap
    colormap0 = get_env('COLORMAP0', default=_colormap0)    # plot colormap for variables out of variable window

    ratio_legend = get_env('RATIO_LEGEND',
        default=_ratio_legend, vartype=float)                               # width ratio between graph and legend
    ncol_legend = get_env('NCOL_LEGEND', default=_ncol_legend, vartype=int) # number of legend columns

    wspace = get_env('WSPACE', default=_wspace, vartype=float)  # plots width space
    hspace = get_env('HSPACE', default=_hspace, vartype=float)  # plots height space

    chi_xs = get_env('CHI_XS', default=_chi_xs)             # cooperativity plot x-scale
    chi_ys = get_env('CHI_YS', default=_chi_ys)             # cooperativity plot y-scale
    x_scale = get_env('X_SCALE', default=_x_scale)          # maximum cooperativity and time of maximum cooperativity plots x-scale
    dtmax_ys = get_env('DTMAX_YS', default=_dtmax_ys)       # time of maximum cooperativity plot y-scale
    chimax_ys = get_env('CHIMAX_YS', default=_chimax_ys)    # maximum cooperativity plot y-scale
    msd_xs = get_env('MSD_XS', default=_msd_xs)             # mean square displacement plot x-scale
    msd_ys = get_env('MSD_YS', default=_msd_ys)             # mean square displacement plot y-scale

    multiply_with_dr = get_env('DRDT', default=True, vartype=bool)  # plot dr*dt rather than dt
    if multiply_with_dr:
        dt_label = '\\tilde{\\nu}_r \\Delta t'                      # dt label
    else: dt_label = '\\Delta t'                                    # dt label

    divide_by_dt = get_env('MSDDT', default=True, vartype=bool) # divide mean square displacement by lag time
    msd_label = r'$<|\Delta \vec{r}(t)|^2>$'                    # mean square displacement plot label
    if divide_by_dt: msd_label += r'$/\Delta t$'

    fit_max = get_env('FIT', default=False, vartype=bool)   # fit maximum cooperativity and time of maximum cooperativity as power law of their x-coordinates on both sides of the variable transition value

    # CALCULATION

    chimsd = ChiMsd(data_dir, naming_simdir, common_attributes,
        parameters_file, var, var_min, var_max,
        excluded_dir=excluded_directories)                              # cooperativities and mean square displacements calculator
    chimsd.calculate_with_msd(naming_cor, attributes_cor, r_min, r_max,
        naming_msd, attributes_msd, init_frame_msd, box_size=box_size,
        multiply_with_dr=multiply_with_dr, divide_by_dt=divide_by_dt)   # calculate cooperativities and mean square displacements

    # PLOT

    colors = {**list_colormap(chimsd.var_list, colormap=colormap),
        **list_colormap(chimsd.var0_list, colormap=colormap0)}  # plot colors hash table
    markers = list_markers(chimsd.time_step_list)               # plot markers hash table
    linestyles = list_linestyles(chimsd.init_frame_list)        # plot linestyles hash table

    # CHI, CHIMAX, DTMAX

    fig_chi = plt.figure()
    fig_chi.subplots_adjust(wspace=wspace, hspace=hspace)
    fig_chi.suptitle(
        r'$N=%.1e, \phi=%1.2f,$' % (N, density)
        + (r'$\tilde{v} = %.2e$' % vzero if mode == 'dr' else
        r'$\tilde{\nu}_r = %.2e$' % dr)
        + r'$, r_{min} = %.2e, r_{max} = %.2e$' % (r_min, r_max))

    gs = GridSpec(2, 3, width_ratios=[1, 1, 2/ratio_legend])

    ax_chi = plt.subplot(gs[:, 0])
    ax_chi.set_xlabel(r'$%s$' % dt_label)
    ax_chi.set_xscale(chi_xs)
    ax_chi.set_ylabel(r'$\chi(\Delta t) = \frac{1}{L^2}$'
        + r'$\int_{r=r_{min}}^{r=r_{max}} dr 2 \pi r %s(r, \Delta t)$'
        % cor_name)
    ax_chi.set_yscale(chi_ys)
    ax_chi.set_title(
        r'$S_{init} = %.1e, S_{max} = %.1e,$' % (init_frame_cor, int_max_cor)
        + r'$N_{cases} = %.1e$' % Ncases_cor)

    ax_dtmax = plt.subplot(gs[0, 1])
    ax_dtmax.set_xscale(x_scale)
    ax_dtmax.set_ylabel(r'$%s^*$' % dt_label)
    ax_dtmax.set_yscale(dtmax_ys)

    ax_chimax = plt.subplot(gs[1, 1], sharex=ax_dtmax)
    plt.setp(ax_dtmax.get_xticklabels(), visible=False)
    ax_chimax.set_xlabel(x_label)
    ax_chimax.set_xscale(x_scale)
    ax_chimax.set_ylabel(r'$\chi(\Delta t^*) = \frac{1}{L^2}$'
        + r'$\int_{r=r_{min}}^{r=r_{max}} dr 2 \pi r %s(r, \Delta t^*)$'
        % cor_name)
    ax_chimax.set_yscale(chimax_ys)

    ax_legend = plt.subplot(gs[:, 2])
    ax_legend.axis('off')
    lines_fig_chi = list(map(
        lambda var_value: Line2D([0], [0], color=colors[var_value], lw=2,
            label='%s = %.0e' % (var_label, var_value)),
        chimsd.var_list))
    lines_fig_chi += [Line2D([0], [0], lw=0, label='')]
    lines_fig_chi += list(map(
        lambda time_step: Line2D([0], [0], marker=markers[time_step],
            color='black', label=r'$dt=%.0e$' % time_step),
        chimsd.time_step_list))
    ax_legend.legend(handles=lines_fig_chi, loc='center', ncol=ncol_legend)

    for dir in chimsd.chi:
        if chimsd.isinvarinterval[dir]:

            time_step_value = chimsd.time_step[dir]
            var_value = chimsd.var_hash[dir]
            plot_parameters = {'color': colors[var_value],
                'marker': markers[time_step_value]}

            ax_chi.plot(*chimsd.chi[dir], **plot_parameters)
            ax_dtmax.plot(x_func(var_value), chimsd.dtmax[dir],
                **plot_parameters)
            ax_chimax.plot(x_func(var_value), chimsd.chimax[dir],
                **plot_parameters)

    if fit_max: # fit maximum cooperativity and time of maximum cooperativity as power law of dt or dr*dt

        # VALUES

        dtmax_lt_varc, dtmax_gt_varc = [], []
        chimax_lt_varc, chimax_gt_varc = [], []
        for dir in chimsd.chi:
            if chimsd.isinvarinterval[dir]:

                var_value = chimsd.var_hash[dir]
                if var_value < var_c:
                    dtmax_lt_varc += [[np.log(x_func(var_value)),
                        np.log(chimsd.dtmax[dir])]]
                    chimax_lt_varc += [[np.log(x_func(var_value)),
                        np.log(chimsd.chimax[dir])]]
                elif var_value > var_c:
                    dtmax_gt_varc += [[np.log(x_func(var_value)),
                        np.log(chimsd.dtmax[dir])]]
                    chimax_gt_varc += [[np.log(x_func(var_value)),
                        np.log(chimsd.chimax[dir])]]

        dtmax_lt_varc = np.transpose(dtmax_lt_varc)
        dtmax_gt_varc = np.transpose(dtmax_gt_varc)
        chimax_lt_varc = np.transpose(chimax_lt_varc)
        chimax_gt_varc = np.transpose(chimax_gt_varc)

        # FIT

        dtmax_lt_varc_slo, dtmax_lt_varc_int, _, _, dtmax_lt_varc_std =\
            st.linregress(*dtmax_lt_varc)
        dtmax_gt_varc_slo, dtmax_gt_varc_int, _, _, dtmax_gt_varc_std =\
            st.linregress(*dtmax_gt_varc)
        chimax_lt_varc_slo, chimax_lt_varc_int, _, _, chimax_lt_varc_std =\
            st.linregress(*chimax_lt_varc)
        chimax_gt_varc_slo, chimax_gt_varc_int, _, _, chimax_gt_varc_std =\
            st.linregress(*chimax_gt_varc)

        # SORT FOR PLOTS

        dtmax_lt_varc[0].sort()
        dtmax_gt_varc[0].sort()
        chimax_lt_varc[0].sort()
        chimax_gt_varc[0].sort()

        # PLOT

        ax_dtmax.plot(np.exp(dtmax_lt_varc[0]),
            np.exp(dtmax_lt_varc_int)
                *(np.exp(dtmax_lt_varc[0])**dtmax_lt_varc_slo),
            color='black', linestyle='-.',
            label=r'$%s \propto (%s)^{%.1e \pm %.1e}$'
                % (ax_dtmax.get_ylabel().replace('$', ''),
                    x_label.replace('$', ''),
                    dtmax_lt_varc_slo, dtmax_lt_varc_std))
        ax_dtmax.plot(np.exp(dtmax_gt_varc[0]),
            np.exp(dtmax_gt_varc_int)
                *(np.exp(dtmax_gt_varc[0])**dtmax_gt_varc_slo),
            color='black', linestyle='--',
            label=r'$%s \propto (%s)^{%.1e \pm %.1e}$'
                % (ax_dtmax.get_ylabel().replace('$', ''),
                    x_label.replace('$', ''),
                    dtmax_gt_varc_slo, dtmax_gt_varc_std))
        ax_chimax.plot(np.exp(chimax_lt_varc[0]),
            np.exp(chimax_lt_varc_int)
                *(np.exp(chimax_lt_varc[0])**chimax_lt_varc_slo),
            color='black', linestyle='-.',
            label=r'$%s \propto (%s)^{%.1e \pm %.1e}$'
                % ('\\chi(\\Delta t^*)', x_label.replace('$', ''),
                    chimax_lt_varc_slo, chimax_lt_varc_std))
        ax_chimax.plot(np.exp(chimax_gt_varc[0]),
            np.exp(chimax_gt_varc_int)
                *(np.exp(chimax_gt_varc[0])**chimax_gt_varc_slo),
            color='black', linestyle='--',
            label=r'$%s \propto (%s)^{%.1e \pm %.1e}$'
                % ('\\chi(\\Delta t^*)', x_label.replace('$', ''),
                    chimax_gt_varc_slo, chimax_gt_varc_std))

        # LABEL

        ax_dtmax.legend()
        ax_chimax.legend()

    # CHI, MSD

    if chimsd.msd:  # if there are mean square displacements to plot

        fig_msd = plt.figure()
        fig_msd.subplots_adjust(wspace=wspace, hspace=hspace)
        fig_msd.suptitle(
            r'$N=%.1e, \phi=%1.2f,$' % (N, density)
            + (r'$\tilde{v} = %.2e$' % vzero if mode == 'dr' else
            r'$\tilde{\nu}_r = %.2e$' % dr)
            + r'$, r_{min} = %.2e, r_{max} = %.2e$' % (r_min, r_max))

        gs = GridSpec(1, 3, width_ratios=[1, 1, 2/ratio_legend])

        ax_msd_chi = plt.subplot(gs[0])
        ax_msd_chi.set_xlabel(r'$%s$' % dt_label)
        ax_msd_chi.set_xscale(chi_xs)
        ax_msd_chi.set_ylabel(r'$\chi(\Delta t) = \frac{1}{L^2}$'
            + r'$\int_{r=r_{min}}^{r=r_{max}} dr 2 \pi r %s(r, \Delta t)$'
            % cor_name)
        ax_msd_chi.set_yscale(chi_ys)
        ax_msd_chi.set_title(
            r'$S_{init} = %.1e$' % init_frame_cor
            + r'$, S_{max} = %.1e,$' % int_max_cor
            + r'$N_{cases} = %.1e$' % Ncases_cor)
        lines_ax_msd_chi = list(map(
            lambda time_step: Line2D([0], [0], marker=markers[time_step],
                color='black', label=r'$dt=%.0e$' % time_step),
            chimsd.time_step_list))
        ax_msd_chi.legend(handles=lines_ax_msd_chi, loc='best')

        ax_msd = plt.subplot(gs[1])
        ax_msd.set_xlabel(r'$%s$' % dt_label)
        ax_msd.set_xscale(msd_xs)
        ax_msd.set_ylabel(msd_label)
        ax_msd.set_yscale(msd_ys)
        ax_msd.set_title(
            r'$S_{max} = %.1e$' % int_max_msd
            + r'$, S_{period} = %.1e$' % int_period_msd)
        lines_ax_msd = list(map(
            lambda init_frame: Line2D([0], [0],
                linestyle=linestyles[init_frame],
                color='black', label=r'$S_{init}=%.1e$' % init_frame),
            chimsd.init_frame_list))
        ax_msd.legend(handles=lines_ax_msd, loc='best')

        ax_legend_msd = plt.subplot(gs[2])
        ax_legend_msd.axis('off')
        lines_fig_msd = list(map(
            lambda var_value: Line2D([0], [0], color=colors[var_value], lw=2,
                label='%s = %.0e' % (var_label, var_value)),
            sorted(chimsd.var_list + chimsd.var0_list)))
        ax_legend_msd.legend(handles=lines_fig_msd,
            loc='center', ncol=ncol_legend)

        for dir in chimsd.chi:
            ax_msd_chi.plot(*chimsd.chi[dir],
                color=colors[chimsd.var_hash[dir]],
                marker=markers[chimsd.time_step[dir]])
        for dir, init_frame in chimsd.msd:
            ax_msd.errorbar(
                chimsd.msd[(dir, init_frame)][:, 0],
                chimsd.msd[(dir, init_frame)][:, 1],
                yerr=chimsd.msd[(dir, init_frame)][:, 2],
                color=colors[chimsd.var_hash[dir]],
                linestyle=linestyles[init_frame])

    # SHOW

    plt.show()
