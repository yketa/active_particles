"""
Module pphiloc plots color histogram of the probability of local density as
function of the local density and either the persistence time at fixed
self-propelling velocity or the self-propelling velocity at fixed persistence
time.

Simulation directories must follow the active_particles.naming.AHB2D naming
standard and input files in simulation directories must follow the
active_particles.naming.VarN naming standard.

Environment modes
-----------------
VARIABLE : string
    Plot x-coordinate variable.
     ______________________________________________________
    | Mode    | Variable                    | x-coordinate |
    |_________|_____________________________|______________|
    | 'dr'    | Rotation diffusion constant | \\tau = 1/dr |
    |_________|_____________________________|______________|
    | 'vzero' | self-propelling velocity    | vzero        |
    |_________|_____________________________|______________|
    DEFAULT: dr
PECLET : bool
    Plot as function of the Péclet number Pe = vzero/dr.
    DEFAULT: True
PHILOCMAX : bool
    Plot most probable packing fraction instead of global system packing
    fraction.
    DEFAULT: False
TITLE : bool
    Display title on figure.
    DEFAULT: True

Environment parameters
----------------------
DATA_DIRECTORY : string
    Data directory.
    DEFAULT: active_particles.naming.sim_directory
EXCLUDE : string
    Simulation directories in DATA_DIRECTORY to exclude from the plot.
    DEFAULT:
PARAMETERS_FILE : string
    Simulation parameters file name.
    DEFAULT: active_particles.naming.parameters_file
DENSITY : float
    Packing fraction of particles.
    DEFAULT: active_particles.plot.pphiloc._density
N : int
    Number of particles.
    DEFAULT: active_particles.plot.pphiloc._N
VZERO ['dr' mode] : float
    Self-propulsion velocity.
    DEFAULT: active_particles.plot.pphiloc._vzero
DR_MIN ['dr' mode] : float
    Minimum rotation diffusion constant.
    DEFAULT: active_particles.plot.pphiloc._dr_min
DR_MAX ['dr' mode] : float
    Maximum rotation diffusion constant.
    DEFAULT: active_particles.plot.pphiloc._dr_max
DR ['vzero' mode] : float
    Rotation diffusion constant.
    DEFAULT: active_particles.plot.pphiloc._dr
VZERO_MIN ['vzero' mode] : float
    Minimum self-propulsion velocity.
    DEFAULT: active_particles.plot.pphiloc._vzero_min
VZERO_MAX ['vzero' mode] : float
    Maximum self-propulsion velocity.
    DEFAULT: active_particles.plot.pphiloc._vzero_max
INITIAL_FRAME : int
    Frame to consider as initial.
    DEFAULT: active_particles.plot.pphiloc._init_frame
INTERVAL_MAXIMUM : int
    Maximum number of frames on which densities are calculated.
    DEFAULT: active_particles.analysis.varn._int_max
BOX_SIZE : float
    Length of the square boxes in which particles are counted.
    DEFAULT: active_particles.analysis.varn._box_size
N_CASES : int
    Number of boxes in each direction in which local density is computed.
    DEFAULT: active_particles.plot.pphiloc._Ncases
N_BINS : int
    Number of bins for the histogram.
    DEFAULT: active_particles.analysis.varn._Nbins
PHIMAX : float
    Maximum local density for histogram.
    DEFAULT: active_particles.plot.pphiloc._phimax
PPHILOC_MIN : float
    Minimum local density probability.
    DEFAULT: active_particles.plot.pphiloc._pphilocmin
PPHILOC_MAX : float
    Maximum local density probability.
    DEFAULT: active_particles.plot.pphiloc._pphilocmin
CONTOURS : int
    Contour level value.
    DEFAULT: active_particles.plot.pphiloc._contours
FONT_SIZE : int
    Font size for the plot.
    DEFAULT: active_particles.plot.pphiloc._font_size
COLORMAP : string
    Plot colormap.
    DEFAULT: active_particles.plot.pphiloc._colormap
PAD : float
    Separation between label and colormap.
    DEFAULT: active_particles.plot.pphiloc._colormap_label_pad
"""

import active_particles.naming as naming

from active_particles.init import get_env, dir_list

from os import environ as envvar
if __name__ == '__main__': envvar['SHOW'] = 'True'
from os.path import join as joinpath

from active_particles.analysis.varn import _int_max, _box_size, _Nbins,\
    _phimax, histogram as get_histogram

import numpy as np
np.seterr(divide='ignore')

import pickle

from collections import OrderedDict

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

# DEFAULT VARIABLES

_dr = 3e-4      # default rotation diffusion constant

_dr_min = 1e-5  # default minimum diffusion rotation constant
_dr_max = 1e-2  # default maximum diffusion rotation constant

_vzero = 1e-2   # default self-propulsion velocity

_vzero_min = 1e-2   # default minimum self-propulsion velocity
_vzero_max = 1e-1   # default maximum self-propulsion velocity

_density = 0.8  # default packing fraction of particles
_N = int(1e5)   # default number of particles

_init_frame = 0 # default frame to consider as initial

_Ncases = 500   # default number of boxes in each direction to compute the local density

_pphilocmin = 1e-4  # default minimum local density probability
_pphilocmax = 1e-1  # default maximum local density probability
_contours = 20      # default contour level value

_font_size = 15             # default font size for the plot
_colormap = 'inferno'       # default plot colormap
_colormap_label_pad = 20    # separation between label and colormap

# FUNCTIONS AND CLASSES

class Philoc:
    """
    Search and read local densities file, compute local densities histogram.
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

    def calculate(self, varN_standard, varN_attributes, Nbins, phimax):
        """
        Calculate local densities histogram.

        Parameters
        ----------
        varN_standard : active_particles.naming._File standard
            Local densities files naming object.
        varN_attributes : hash table
            Attributes to be displayed in local densities file names.
        Nbins : int
            Number of bins for the histogram.
        phimax : float
            Maximum local density for histogram.
        """

        self.varN_standard = varN_standard
        self.varN_attributes = varN_attributes
        self.Nbins = Nbins
        self.phimax = phimax

        self.time_step = {}     # hash table of directories' simulation time step
        self.histogram3D = []   # local densities histogram
        self.philocmax = {}     # hash table of most probable local density with directory name as keys

        for dir in sorted(self.dirs):
            try:
                varN_filename, = self.varN_standard.get_files(
                    directory=joinpath(self.data_dir, dir),
                    **self.varN_attributes)
            except ValueError: continue

            with open(
                joinpath(self.data_dir, dir, self.parameters_file), 'rb')\
                as param_file:
                self.time_step[dir] = pickle.load(param_file)['time_step']  # simulation parameters hash table

            with open(joinpath(self.data_dir, dir, varN_filename),
                'rb') as varN_file:

                var_value = np.full(self.Nbins,
                    fill_value=self.var_hash[dir])  # plot variable value

                densities = pickle.load(varN_file)  # list of local densities
                bins, histogram = get_histogram(densities,
                    self.Nbins, self.phimax)        # histogram of local densities with corresponding bins

                histogram = np.log10(histogram)

                histogram3D_dir = np.transpose(
                    [var_value, bins, histogram]).tolist()
                self.histogram3D += histogram3D_dir
                _, self.philocmax[dir], _ = max(histogram3D_dir,
                    key=lambda el: el[2])

        self.histogram3D = np.transpose(self.histogram3D)
        self.time_step_list = sorted(OrderedDict.fromkeys(
            self.time_step.values()))   # list of time steps

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLES DEFINITIONS

    mode = get_env('VARIABLE', default='dr')                # plotting variable
    peclet = get_env('PECLET', default=True, vartype=bool)  # display Péclet number rather than mode variable

    if mode == 'dr':

        vzero = get_env('VZERO', default=_vzero, vartype=float) # self-propulsion velocity
        attributes = {'vzero': vzero}                           # attributes displayed in filenames

        var = 'dr'                                                  # plot variable
        var_min = get_env('DR_MIN', default=_dr_min, vartype=float) # minimum rotation diffusion constant
        var_max = get_env('DR_MAX', default=_dr_max, vartype=float) # maximum rotation diffusion constant

        if peclet:
            x_func = lambda x: np.log10(vzero/x)    # x-coordinate as function of plot variable
        else:
            x_func = lambda x: np.log10(1/x)        # x-coordinate as function of plot variable
            x_label = r'$\log(1/\tilde{\nu}_r)$'    # x-coordinate label

    elif mode == 'vzero':

        dr = get_env('DR', default=_dr, vartype=float)  # rotation diffusion constant
        attributes = {'dr': dr}                         # attributes displayed in filenames

        var = 'vzero'                                                       # plot variable
        var_min = get_env('VZERO_MIN', default=_vzero_min, vartype=float)   # minimum self-propulsion velocity
        var_max = get_env('VZERO_MAX', default=_vzero_max, vartype=float)   # maximum self-propulsion velocity

        if peclet:
            x_func = lambda x: np.log10(x/dr)   # x-coordinate as function of plot variable
        else:
            x_func = lambda x: np.log10(x)      # x-coordinate as function of plot variable
            x_label = r'$\log(\tilde{v})$'      # x-coordinate label

    else: raise ValueError('Mode %s is not known.' % mode)  # mode is not known

    if peclet: x_label = r'$\log(Pe)$'  # x-coordinate label

    data_dir = get_env('DATA_DIRECTORY', default=naming.sim_directory)  # data directory
    excluded_directories = get_env('EXCLUDE', default='')               # directories to exclude

    parameters_file = get_env('PARAMETERS_FILE',
        default=naming.parameters_file) # simulations parameters file name

    density = get_env('DENSITY', default=_density, vartype=float)   # packing fraction of particles
    N = get_env('N', default=_N, vartype=int)                       # number of particles

    init_frame = get_env('INITIAL_FRAME', default=_init_frame, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=_int_max, vartype=int)	# maximum number of frames on which densities are calculated

    box_size = get_env('BOX_SIZE', default=_box_size, vartype=float)    # length of the square boxes in which particles are counted

    Ncases = get_env('N_CASES', default=_Ncases, vartype=int)   # number of boxes in each direction in which local density is computed

    # NAMING

    attributes = {**attributes, 'density': density, 'N': N,
        'init_frame': init_frame, 'int_max': int_max, 'Ncases': Ncases,
        'box_size': box_size}
    naming_varN = naming.VarN()     # varN naming object
    naming_simdir = naming.AHB2D()  # simulation directory naming object

    # PLOT PARAMETERS

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=float)
    mp.rcParams.update({'font.size': font_size})    # font size for the plot

    pad = get_env('PAD', default=_colormap_label_pad, vartype=float)    # separation between label and colormap

    colormap = get_env('COLORMAP', default=_colormap)   # plot colormap

    Nbins = get_env('N_BINS', default=_Nbins, vartype=int)      # number of bins for the histogram
    phimax = get_env('PHIMAX', default=_phimax, vartype=float)  # maximum local density for histogram

    pphilocmin = np.log10(
        get_env('PPHILOC_MIN', default=_pphilocmin, vartype=float)) # minimum local density probability
    pphilocmax = np.log10(
        get_env('PPHILOC_MAX', default=_pphilocmax, vartype=float)) # maximum local density probability
    contours = get_env('CONTOURS', default=_contours, vartype=int)  # contour level value

    philocmax = get_env('PHILOCMAX', default=False, vartype=bool)   # plot most probable packing fraction instead of global system packing fraction

    # CALCULATION

    philoc = Philoc(data_dir, naming_simdir, attributes, parameters_file,
        var, var_min, var_max, excluded_dir=excluded_directories)   # local densities histogram calculator
    philoc.calculate(naming_varN, attributes, Nbins, phimax)        # calculate local densities histogram

    philoc.histogram3D[0] = x_func(philoc.histogram3D[0])                   # apply x_func to plot variable
    philoc.histogram3D[2][philoc.histogram3D[2] < pphilocmin] = pphilocmin  # set minimum histogram value as pphilocmin

    # PLOT

    fig, ax = plt.subplots()

    cmap = plt.get_cmap(colormap)
    norm = colors.Normalize(vmin=pphilocmin, vmax=pphilocmax)
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = mp.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
        orientation='vertical')
    cb.set_label(r'$\log P(\phi_{loc})$', labelpad=pad, rotation=270)

    ax.tricontourf(*philoc.histogram3D, contours, cmap=cmap, norm=norm) # local density histogram
    if philocmax:
        philocmax_line = [[x_func(philoc.var_hash[dir]), philoc.philocmax[dir]]
            for dir in philoc.dirs
            if dir in philoc.var_hash and dir in philoc.philocmax]
        philocmax_line.sort(key=lambda el: el[0])
        philocmax_line = np.array(philocmax_line)
        ax.plot(philocmax_line[:, 0], philocmax_line[:, 1],
            linestyle='--', color='red', linewidth=4)   # most probable packing fraction line
    else:
        ax.plot(
            philoc.histogram3D[0],
            np.full(philoc.histogram3D.shape[1], fill_value=density),
            linestyle='--', color='red', linewidth=4)   # global system packing fraction line

    title = (
        r'$N=%.1e, \phi=%1.2f,$' % (N, density)
        + (r'$\tilde{v} = %.2e$' % vzero if mode == 'dr' else
        r'$\tilde{\nu}_r = %.2e$' % dr)
        + '\n' + r'$S_{init} = %.1e, S_{max} = %.1e,$' % (init_frame, int_max)
        + r'$N_{cases} = %.1e, r_{max} = %.1e$' % (Ncases, box_size))
    if get_env('TITLE', default=True, vartype=bool): ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(r'$\phi_{loc}$')

    plt.show()
