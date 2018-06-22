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
X_VARIABLE : string
    Plot x-coordinate variable.
     _____________________________________________________________
    | Mode    | Variable                    | x-coordinate        |
    |_________|_____________________________|_____________________|
    | 'dr'    | Rotation diffusion constant | \\tau = \\log(1/dr) |
    |_________|_____________________________|_____________________|
    | 'vzero' | self-propelling velocity    | vzero               |
    |_________|_____________________________|_____________________|
    DEFAULT: dr

Environment parameters
----------------------
DATA_DIRECTORY : string
    Data directory.
    DEFAULT: active_particles.naming.sim_directory
EXCLUDE : string
    Simulation directories in DATA_DIRECTORY to exclude from the plot.
    DEFAULT:
PARAMETERS_FILE : string
    Simulation parameters file.
    DEFAULT: DATA_DIRECTORY/active_particles.naming.parameters_file
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
"""

import active_particles.naming as naming

from active_particles.init import get_env

from os import getcwd
from os import environ as envvar
envvar['SHOW'] = 'True'
from os.path import join as joinpath

from active_particles.analysis.varn import _int_max, _box_size, _Nbins,\
    _phimax, histogram as get_histogram

import numpy as np
np.seterr(divide='ignore')

import pickle

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

# DEFAULT VARIABLES

_dr = 3e-4      # default rotation diffusion constant

_dr_min = 1e-5  # default minimum diffusion rotation constant
_dr_max = 1e-2  # default maximum diffusion rotation constant

_vzero = 1e-2  # default self-propulsion velocity

_vzero_min = 1e-2   # default minimum self-propulsion velocity
_vzero_max = 1e-1   # default maximum self-propulsion velocity

_density = 0.8  # default packing fraction of particles
_N = int(1e5)   # default number of particles

_init_frame = 0 # default frame to consider as initial

_Ncases = 500   # default number of boxes in each direction to compute the local density

_pphilocmin = 1e-4  # default minimum local density probability
_pphilocmax = 1e-1  # default maximum local density probability
_contours = 20      # default contour level value

_font_size = 15         # default font size for the plot
_colormap = 'inferno'   # default plot colormap

if __name__ == '__main__':  # executing as script

    # VARIABLES DEFINITIONS

    mode = get_env('X_VARIABLE', default='dr')  # plotting mode

    if mode == 'dr':

        vzero = get_env('VZERO', default=_vzero, vartype=float) # self-propulsion velocity
        attributes = {'vzero': vzero}                           # attributes displayed in filenames

        var = 'dr'                                                  # plot variable
        var_min = get_env('DR_MIN', default=_dr_min, vartype=float) # minimum rotation diffusion constant
        var_max = get_env('DR_MAX', default=_dr_max, vartype=float) # maximum rotation diffusion constant

        x_func = lambda x: np.log10(1/x)    # x-coordinate as function of plot variable

        x_label = r'$\log( \tau_r \equiv \tilde{\nu}_r^{-1})$'  # x-coordinate label

    elif mode == 'vzero':

        dr = get_env('DR', default=_dr, vartype=float)  # rotation diffusion constant
        attributes = {'dr': dr}                         # attributes displayed in filenames

        var = 'vzero'                                                       # plot variable
        var_min = get_env('VZERO_MIN', default=_vzero_min, vartype=float)   # minimum self-propulsion velocity
        var_max = get_env('VZERO_MAX', default=_vzero_max, vartype=float)   # maximum self-propulsion velocity

        x_func = lambda x: x    # x-coordinate as function of plot variable

        x_label = r'$\tilde{v}$'  # x-coordinate label

    else: raise ValueError('Mode %s is not known.' % mode)  # mode is not known

    data_dir = get_env('DATA_DIRECTORY', default=naming.sim_directory)  # data directory
    excluded_directories = get_env('EXCLUDE', default='')               # directories to exclude

    parameters_file = get_env('PARAMETERS_FILE',
        default= naming.parameters_file)    # simulations parameters file

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

    colormap = get_env('COLORMAP', default=_colormap)   # plot colormap

    Nbins = get_env('N_BINS', default=_Nbins, vartype=int)      # number of bins for the histogram
    phimax = get_env('PHIMAX', default=_phimax, vartype=float)  # maximum local density for histogram

    pphilocmin = np.log10(
        get_env('PPHILOC_MIN', default=_pphilocmin, vartype=float)) # minimum local density probability
    pphilocmax = np.log10(
        get_env('PPHILOC_MAX', default=_pphilocmax, vartype=float)) # maximum local density probability
    contours = get_env('CONTOURS', default=_contours, vartype=int)  # contour level value

    # CALCULATION

    dirs = [dir
        for dir in naming_simdir.get_files(directory=data_dir, **attributes)
        if not(dir in excluded_directories) and
        naming_simdir.get_data(dir, var)[0] >= var_min and
        naming_simdir.get_data(dir, var)[0] <= var_max] # directories to consider

    histogram3D = [] # plot histogram
    for dir in sorted(dirs):
        try:
            varN_filename, = naming_varN.get_files(
                directory=joinpath(data_dir, dir), **attributes)
        except ValueError: continue

        with open(joinpath(data_dir, dir, parameters_file),
            'rb') as param_file, open(joinpath(data_dir, dir, varN_filename),
            'rb') as varN_file:

            parameters = pickle.load(param_file)                            # simulation parameters
            var_value = np.full(Nbins, fill_value=x_func(parameters[var]))  # plot variable value

            densities = pickle.load(varN_file)                          # list of local densities
            bins, histogram = get_histogram(densities, Nbins, phimax)   # histogram of local densities with corresponding bins

            histogram = np.log10(histogram)
            histogram[histogram < pphilocmin] = pphilocmin

            histogram3D += np.transpose([var_value, bins, histogram]).tolist()

    histogram3D = np.array(histogram3D)

    # PLOT

    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20)

    cmap = plt.get_cmap(colormap)
    norm = colors.Normalize(vmin=pphilocmin, vmax=pphilocmax)
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = mp.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
        orientation='vertical')
    cb.set_label(r'$\log P(\phi_{loc})$', labelpad=20, rotation=270)

    ax.tricontourf(histogram3D[:, 0], histogram3D[:, 1], histogram3D[:, 2],
        contours, cmap=cmap, norm=norm)             # local density histogram
    ax.plot(
        histogram3D[:, 0], np.full(histogram3D.shape[0], fill_value=density),
        linestyle='--', color='red', linewidth=4)   # average density line

    ax.set_title(
        r'$N=%.1e, \phi=%1.2f,$' % (N, density)
        + (r'$\tilde{v} = %.2e$' % vzero if mode == 'dr' else
        r'$\tilde{\nu}_r = %.2e$' % dr)
        + '\n' + r'$S_{init} = %.1e, S_{max} = %.1e,$' % (init_frame, int_max)
        + r'$N_{cases} = %.1e, r_{max} = %.1e$' % (Ncases, box_size))
    ax.set_xlabel(x_label)
    ax.set_ylabel(r'$\phi_{loc}$')

    plt.show()
