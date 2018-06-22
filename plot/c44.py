"""
Module c44 plots strain correlations projected on cos(4 \\theta) for different
lag times, which we call C44.

Input files in simulation directory must follow the active_particles.naming.Css
naming standard.

Environment modes
-----------------
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
R_CUT : float
    Cut-off radius for coarse graining function in terms of particles' mean
    radius.
    DEFAULT: active_particles.analysis.css._r_cut
SIGMA : float
    Length scale of the spatial extent of the coarse graining function.
    DEFAULT: cut-off radius
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
    DEFAULT: active_particles.plot.c44._y_min
Y_MAX : float
    Maximum C44 value.
    DEFAULT: active_particles.plot.c44._y_max
R_MIN : float
    Minimum radius over average particle separation value.
    DEFAULT: active_particles.plot.c44._r_min
R_MAX : float
    Maximum radius over average particle separation value.
    DEFAULT: active_particles.plot.c44._r_max
POINTS_X : int
    Number of radii at which to evaluate integrated strain correlation.
    DEFAULT: active_particles.plot.c44._points_x
POINTS_THETA : int
    Number of angles to evaluate integrated strain correlation.
    DEFAULT: active_particles.plot.c44._points_theta
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
COLORMAP : string
    Plot colormap.
    DEFAULT: active_particles.plot.c44._colormap
SLOPE [FITTING_LINE mode] : float
	Initial slope for fitting line.
	DEFAULT: active_particles.plot.c44._slope0
SLOPE_MIN [FITTING_LINE mode] : float
	Minimum slope for fitting line.
	DEFAULT: active_particles.plot.c44._slope_min
SLOPE_MAX [FITTING_LINE mode] : float
	Maximum slope for fitting line.
	DEFAULT: active_particles.plot.c44._slope_max
"""

import active_particles.naming as naming

from active_particles.init import get_env

from os import getcwd
from os import environ as envvar
envvar['SHOW'] = 'True'
from os.path import join as joinpath

from active_particles.analysis.correlations import CorGrid
from active_particles.plot.plot import list_colormap
from active_particles.plot.mpl_tools import FittingLine

from active_particles.analysis.css import _r_cut

from math import ceil

import numpy as np

import pickle

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

class Css2DtoC44:
    """
    Calculates C44 as projection of 2D shear strain correlations on
    cos(4\\theta).
    """

    def __init__(self, box_size, points_x, points_theta,
        r_min, r_max):
        """
        Sets parameters for C44 integration.

        Parameters
        ----------
        box_size : float
            Size of the square box to consider.
        points_x : int
            Number of radii at which to evaluate integrated strain correlation.
        points_theta : int
            Number of angles to evaluate integrated strain correlation.
        r_min : float
            Minimum radius.
        r_max : float
            Maximum radius.
        """

        self.box_size = box_size

        self.points_x = points_x
        self.points_theta = points_theta

        self.r_min = r_min
        self.r_max = r_max
        self.c44_x = np.linspace(self.r_min, self.r_max, self.points_x)

    def get_C44(self, Css2D):
        """
        From 2D strain correlations Css2D, returns values of C44 at
        self.points_x radii between r_max and r_min.

        Parameters
        ----------
        Css2D : 2D array-like
            Shear strain correlation grid.

        Returns
        -------
        C44 : 2D Numpy array
            List of [r, C44(r)].
        """

        self.css2Dgrid = CorGrid(Css2D, self.box_size)  # shear strain 2D CorGrid object
        self.c44 = np.array(list(map(
            lambda r: [r, self.css2Dgrid.integrate_over_angles(r,
            projection=lambda theta: np.cos(4*theta)/np.pi,
            points_theta=self.points_theta)],
            self.c44_x)))

        return self.c44

# DEFAULT VARIABLES

_font_size = 10     # default plot font size
_marker_size = 20   # default plot marker size

_ratio_legend = 7   # default width ratio between legend and figure
_ncol_legend = 1    # default number of columns for the legend

_wspace = 0.2   # default plots width space
_hspace = 0.05  # default plots height space

_colormap = 'jet'   # default plot colormap

_slope0 = -2    # default initial slope for fitting line slider
_slope_min = -5 # default minimum slope for fitting line slider
_slope_max = 0  # default maximum slope for fitting line slider

_points_x = 100     # default number of radii at which to evaluate integrated strain correlation
_points_theta = 100 # default number of angles to evaluate integrated strain correlation

_y_min = 1e-4   # default minimum C44 value
_y_max = 2e-1   # default maximum C44 value
_r_min = 1      # default minimum radius over average particle separation value
_r_max = 20     # default maximum radius over average particle separation value

if __name__ == '__main__':  # executing as script

    # VARIABLES DEFINITIONS

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

    Ncases = get_env('N_CASES', default=ceil(np.sqrt(parameters['N'])),
		vartype=int)	# number of boxes in each direction for which shear strain and displacement vorticity grids have been computed

    Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame

    box_size = get_env('BOX_SIZE', default=parameters['box_size'],
		vartype=float)									# size of the square box to consider
    centre = (get_env('X_ZERO', default=0, vartype=float),
		get_env('Y_ZERO', default=0, vartype=float))	# centre of the box

    # NAMING

    attributes = {'density': parameters['density'], 'N': parameters['N'],
        'init_frame': init_frame, 'int_max': int_max, 'Ncases': Ncases,
        'r_cut': r_cut, 'sigma': sigma, 'box_size': box_size,
        'x_zero': centre[0], 'y_zero': centre[1]}   # attributes displayed in filenames
    naming_Css = naming.Css()                       # Css naming object

    # PLOT PARAMETERS

    y_min = get_env('Y_MIN', default=_y_min, vartype=float) # minimum C44 value
    y_max = get_env('Y_MAX', default=_y_max, vartype=float) # maximum C44 value

    r_min = get_env('R_MIN', default=_r_min, vartype=float) # minimum radius over average particle separation value
    r_max = get_env('R_MAX', default=_r_max, vartype=float) # maximum radius over average particle separation value

    font_size = get_env('FONT_SIZE', default=_font_size, vartype=int)       # plot font size
    marker_size = get_env('MARKER_SIZE', default=_marker_size, vartype=int) # plot marker size
    mpl.rcParams.update({'font.size': font_size, 'lines.markersize': marker_size})

    ratio_legend = get_env('RATIO_LEGEND', default=_ratio_legend,
        vartype=float)                                                      # width ratio between legend and figure
    ncol_legend = get_env('NCOL_LEGEND', default=_ncol_legend, vartype=int) # number of columns for the legend

    wspace = get_env('WSPACE', default=_wspace, vartype=float)  # plots width space
    hspace = get_env('HSPACE', default=_hspace, vartype=float)  # plots height space

    colormap = get_env('COLORMAP', default=_colormap)   # plot colormap

    # CALCULATION PARAMETERS

    points_x = get_env('POINTS_X', default=_points_x, vartype=int)              # number of radii at which to evaluate integrated strain correlation
    points_theta = get_env('POINTS_THETA', default=_points_theta, vartype=int)  # number of angles to evaluate integrated strain correlation

    # CALCULATION

    files = sorted(naming_Css.get_files(directory=data_dir, **attributes))  # files corresponding to parameters
    dt_list = np.array(list(map(
        lambda file: naming_Css.get_data(file, 'dt'), files))).flatten()    # list of lag times corresponding to files

    toC44 = Css2DtoC44(box_size, points_x, points_theta,
        av_p_sep*r_min, av_p_sep*r_max) # Css2D to C44 conversion object
    C44 = {}                            # hash table of tables of C44 at X with lag times as keys
    for file, dt in zip(files, dt_list):
        with open(joinpath(data_dir, file), 'rb') as Css_dump_file:
            _, Css2D = pickle.load(Css_dump_file)
        C44[dt] = toC44.get_C44(Css2D)

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

    for dt in dt_list:
        ax.loglog(C44[dt][:, 0]/av_p_sep, C44[dt][:, 1],
            color=colors[dt], label=r'$\Delta t = %.0e$'
            % (parameters['period_dump']*parameters['time_step']*dt))

    if get_env('TEMPERATURE', default=False, vartype=bool): # TEMPERATURE mode
        nD0 = (2*parameters['kT']*parameters['N'])/(
            parameters['damp_bro']*parameters['a']*(parameters['box_size']**2))
        legend0 = [mpatches.Patch(color='none',
            label=r'$nD_0 = \frac{2 k_B T N}{\lambda a L^2} = %.2e$' % nD0)]
        legend0 += list(map(lambda dt: Line2D([0], [0], color=colors[dt],
            label=r'$nD_0\Delta t = %.2e$'
            % (dt*parameters['time_step']*parameters['period_dump']*nD0)),
            dt_list))
    else:
        nD0 = (parameters['N']*(parameters['vzero']**2))/(2*parameters['dr']
            *(parameters['box_size']**2))
        legend0 = [mpatches.Patch(color='none', label=
            r'$nD_0 = \frac{Nv_0^2}{2\nu_r L^2} = %.2e$' % nD0)]
        legend0 += list(map(lambda dt: Line2D([0], [0], color=colors[dt],
            label=r'$nD_0\Delta t = %.2e$'
            % (dt*parameters['time_step']*parameters['period_dump']*nD0)),
            dt_list))
    legend0 += [Line2D([0], [0], lw=0, label='')]
    leg.legend(handles=legend0, loc='center', ncol=ncol_legend)

    title = r'$N=%.2e, \phi=%1.2f,$' % (parameters['N'], parameters['density'])
    if get_env('TEMPERATURE', default=False, vartype=bool): # TEMPERATURE mode
        title += r'$kT=%.2e, k=%.2e$' % (parameters['kT'], parameters['k'])
    else:
        title += r'$\tilde{v}=%.2e, \tilde{\nu}_r=%.2e$' % (
            parameters['vzero'], parameters['dr'])
    title += '\n'
    title += r'$S_{init}=%.2e, S_{max}=%.2e,$' % (init_frame, int_max)
    title += r'$N_{cases}=%.2e, r_{cut}=%.2e, \sigma=%.2e$' % (Ncases, r_cut,
        sigma)
    title += '\n'
    title += r'$N_r=%.2e, N_{\theta}=%.2e$' % (points_x, points_theta)
    if 'BOX_SIZE' in envvar:
        title += r'$, L = %.3e$' % box_size
    if 'X_ZERO' in envvar or 'Y_ZERO' in envvar:
        title += r'$, x_0=%.3e, y_0=%.3e$' % centre
    fig.suptitle(title)

    plt.show()
