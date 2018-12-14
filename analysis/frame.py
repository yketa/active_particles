"""
Module frame renders images of the 2D system, with particles rendered as
circles and arrows indicating either displacement or instantaneous velocitiy.

Files are saved according to the active_particles.naming.Velocity or
active_particles.naming.Trajectory or active_particles.naming.Displacement
standard or with custom name.

Environment modes
-----------------
MODE : string
    Plotting mode.
     _________________________________________________________________________
    | Mode           | Arrow direction        | Arrow length | Particle color |
    |________________|________________________|______________|________________|
    | 'd2min'        | None                   | None         | Amplitude of   |
    |                |                        |              | nonaffine      |
    |                |                        |              | squared        |
    |                |                        |              | displacement   |
    |________________|________________________|______________|________________|
    | 'velocity'     | Instantaneous velocity | Relative to  | Amplitude of   |
    |                | direction              | diameter     | velocity       |
    |________________|________________________|______________|________________|
    | 'trajectory'   | Displacement direction | Displacement | Only tracer    |
    |                |                        | amplitude    | particle       |
    |________________|________________________|______________|________________|
    | 'displacement' | Displacement direction | Relative to  | Amplitude of   |
    |                |                        | diameter     | displacement   |
    |________________|________________________|______________|________________|
    DEFAULT: displacement
PLOT : bool
    Plot single frame.
    DEFAULT: False
MOVIE : bool
    Make movie out of several plotted frames.
    DEFAULT: False
SHOW [PLOT mode] : bool
    Show figure.
    DEFAULT: False
SAVE [PLOT mode] : bool
    Save figure.
    DEFAULT: False
TRACER_PARTICLE ['trajecotry' mode] : bool
    Display tracer particle.
    DEFAULT: True
SUPTITLE : bool
    Display suptitle on figures.
    DEFAULT: True

Environment parameters
----------------------
DATA_DIRECTORY : string
    Data directory.
    DEFAULT: current working directory
PARAMETERS_FILE : string
    Simulation parameters file.
    DEFAULT: DATA_DIRECTORY/active_particles.naming.parameters_file
WRAPPED_FILE : string
    Wrapped trajectory file. (.gsd)
    DEFAULT: DATA_DIRECTORY/active_particles.naming.wrapped_trajectory_file
UNWRAPPED_FILE : string
    Unwrapped trajectory file. (.dat)
    NOTE: .dat files defined with active_particles.dat
    DEFAULT: DATA_DIRECTORY/active_particles.naming.unwrapped_trajectory_file
INITIAL_FRAME : int
    [PLOT mode] Frame to render.
    [MOVIE mode] Initial frame to render.
    NOTE: FRAME < 0 will be interpreted as the frame to render being the middle
          frame of the simulation.
    DEFAULT: -1
FINAL_FRAME [MOVIE mode] : int
    Final movie frame.
    DEFAULT: Final simulation frame.
FRAME_PERIOD [MOVIE mode] : int
    Frame rendering period.
    DEFAULT: active_particles.analysis.frame._frame_per
FRAME_MAXIMUM : int
    Maximum number of frames.
    DEFAULT: active_particles.analysis.frame._frame_max
DT : int
    Lag time for displacement.
    NOTE: [PLOT mode] TIME < 0 will be interpreted as a lag time corresponding
                      to the total number of simulation frames - FRAME + TIME.
          [MOVIE mode] TIME < 0 will be interpreted as a lag time
                       corresponding to the minimum distance between frames.
    DEFAULT: -1
BOX_SIZE : float
    Length of the square box to render.
    DEFAULT: simulation box size
X_ZERO : float
    1st coordinate of the centre of the square box to render.
    DEFAULT: 0
Y_ZERO : float
    2nd coordinate of the centre of the square box to render.
    DEFAULT: 0
V_MIN : float
    Minimum value of the colorbar.
    DEFAULT: ['d2min' mode] minimum nonaffine squared displacement value
             ['velocity' mode] 10^{E(log(||\\vec{v}||))-2*V(log(||\\vec{v}||))}
             [other modes] 10^{E(log(||\\vec{u}||))-2*V(log(||\\vec{u}||))}
    NOTE: Colorbar is represented in logarithmic scale so V_MIN > 0.
V_MAX : float
    Maximum value of the colorbar.
    DEFAULT: ['d2min' mode] maximum nonaffine squared displacement value
             ['velocity' mode] 10^{E(log(||\\vec{v}||))+2*V(log(||\\vec{v}||))}
             [other modes] 10^{E(log(||\\vec{u}||))+2*V(log(||\\vec{u}||))}
    NOTE: Colorbar is represented in logarithmic scale so V_MAX > 0.
ARROW_WIDTH : float
    Width of the arrows.
    DEFAULT: active_particles.analysis.frame._arrow_width
HEAD_WIDTH : float
    Width of the arrows' head.
    DEFAULT: active_particles.analysis.frame._arrow_head_width
HEAD_LENGTH : float
    Length of the arrows' head.
    DEFAULT: active_particles.analysis.frame._arrow_head_length
FRAME_VERTICAL_SIZE : float
    Vertical size of the frame (in inches).
    DEFAULT: active_particles.analysis.frame._frame_ver
FRAME_HORIZONTAL_SIZE : float
    Horizontal size of the frame (in inches).
    DEFAULT: active_particles.analysis.frame._frame_hor
FRAME_DEFINITION [SAVE mode] : float
    Definition of image (in dots per inches (dpi)).
    DEFAULT: active_particles.analysis.frame._frame_def
FONT_SIZE : int
    Font size.
    DEFAULT: active_particles.analysis.frame._font_size
PAD : float
    Separation between label and colormap.
    DEFAULT: active_particles.analysis.frame._colormap_label_pad
FIGURE_NAME [SAVE mode] : string
    Custom figure name.
    DEFAULT: according to naming standard

Output
------
> Prints execution time.
[PLOT mode]
> Plots system according to plotting mode.
[MOVIE mode]
> Creates movie from frame concatenation according to plotting mode.
[SHOW mode]
> Displays figure.
[SAVE mode]
> Saves figure in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output, mkdir
from active_particles.dat import Dat, Gsd
from active_particles.maths import normalise1D, amplogwidth

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

import sys

from math import ceil

import pickle

import numpy as np
np.seterr(divide='ignore')

import matplotlib as mpl
if not(get_env('SHOW', default=False, vartype=bool)):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as ColorsNormalise
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

from datetime import datetime

from collections import OrderedDict

import subprocess

# DEFAULT VARIABLES

_frame_per = 1      # default frame rendering period
_frame_max = 1000   # default maximum number of frames

_frame_ver = 12 # default vertical size of the frames (in inches)
_frame_hor = 16 # default horizontal size of the frames (in inches)
_frame_def = 80 # default definition of images (in dots per inches (dpi))

_arrow_width = 1e-3                         # default width of the arrows
_arrow_head_width = _arrow_width*3e2        # default width of the arrows' head
_arrow_head_length = _arrow_head_width*1.5  # default length of the arrows' head

_font_size = 15 # font size

_colormap_label_pad = 30    # default separation between label and colormap

# FUNCTIONS AND CLASSES

class _Frame:
    """
    This class is designed as the superclass of all other plotting classes
    specific to each mode. It initialises the figure and provides methods to
    plot circles representing particles and arrows at the particles' positions,
    and to add a colorbar.
    """

    def __init__(self, w_traj, frame, box_size, centre, arrow_width,
        arrow_head_width, arrow_head_length):
        """
        Initialises figure.

        Parameters
        ----------
        w_traj : active_particles.dat.Gsd
    		Wrapped trajectory object.
        frame : int
            Frame to render.
        box_size : float
            Length of the square box to render.
        centre : 2-uple like
            Centre of the box to render.
        arrow_width : float
            Width of the arrows.
        arrow_head_width : float
            Width of the arrows' head.
        arrow_head_length : float
            Length of the arrows' head.
        """

        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim([-1.1*box_size/2, 1.1*box_size/2])
        self.ax.set_xlabel(r'$x$')
        self.ax.set_ylim([-1.1*box_size/2, 1.1*box_size/2])
        self.ax.set_ylabel(r'$y$')
        self.ax.set_aspect('equal')

        self.positions = w_traj.position(frame, centre=centre)   # particles' positions at frame frame with centre as centre of frame
        self.diameters = w_traj.diameter(frame)                  # particles' diameters at frame frame

        self.particles = [particle for particle in range(len(self.positions))
            if (np.abs(self.positions[particle]) <= box_size/2).all()]  # particles inside box of centre centre and length box_size

        self.arrow_width = arrow_width
        self.arrow_head_width = arrow_head_width
        self.arrow_head_length = arrow_head_length

    def __del__(self):
        """
        Closes figure.
        """

        plt.close(self.fig)

    def draw_circle(self, particle, color='black', fill=False):
        """
        Draws circle at particle's position with particle's diameter.

        Parameters
        ----------
        particle : int
            Particle index.
        color : any matplotlib color
            Circle color. (default: 'black')
        fill : bool
            Filling the circle with same color. (default: False)
        """

        circle = plt.Circle(self.positions[particle],
            self.diameters[particle]/2, color=color, fill=fill, zorder=0)   # circle artist representing particle
        self.ax.add_artist(circle)

    def draw_arrow(self, particle, dx, dy, color='black'):
        """
        Draws arrow starting from particle's position.

        Parameters
        ----------
        particle : int
            Particle index.
        dx : float
            Arrow length in x-direction.
        dy : float
            Arrow length in y-direction.
        color : any matplotlib color
            Arrow color. (default: 'black')
        """

        length = np.sqrt(dx**2 + dy**2) # length of arrow
        if length == 0: return
        self.ax.arrow(*self.positions[particle], dx, dy, color=color,
            width=length*self.arrow_width,
            head_width=length*self.arrow_head_width,
            head_length=length*self.arrow_head_length, zorder=1)

    def colorbar(self, vmin, vmax, cmap=plt.cm.jet):
        """
        Adds colorbar to plot.

        Parameters
        ----------
        vmin : float
            Minimum value of the colorbar.
        vmax : float
            Maximum value of the colorbar.
        cmap : matplotlib colorbar
            Matplotlib colorbar to be used. (default: matplotlib.pyplot.cm.jet)
        """

        vNorm = ColorsNormalise(vmin=vmin, vmax=vmax)
        self.scalarMap = ScalarMappable(norm=vNorm, cmap=cmap)

        self.colormap_ax = make_axes_locatable(self.ax).append_axes('right',
            size='5%', pad=0.05)
        self.colormap = mpl.colorbar.ColorbarBase(self.colormap_ax, cmap=cmap,
            norm=vNorm, orientation='vertical')

class D2min(_Frame):
    """
    Plotting class specific to 'd2min' mode.
    """

    def __init__(self, u_traj, w_traj, frame, box_size, centre, arrow_width,
        arrow_head_width, arrow_head_length, pad, dt=0, **kwargs):
        """
        Initialises and plots figure.

        Parameters
        ----------
        u_traj : active_particles.dat.Dat
    		Unwrapped trajectory object.
        w_traj : active_particles.dat.Gsd
    		Wrapped trajectory object.
        frame : int
            Frame to render.
        box_size : float
            Length of the square box to render.
        centre : 2-uple like
            Centre of the box to render.
        arrow_width : float
            Width of the arrows.
            NOTE: not used.
        arrow_head_width : float
            Width of the arrows' head.
            NOTE: not used.
        arrow_head_length : float
            Length of the arrows' head.
            NOTE: not used.
        pad : float
            Separation between label and colormap.
        dt : int
            Lag time for displacement. (default=0)

        Optional keyword parameters
        ---------------------------
        vmin : float
            Minimum value of the colorbar.
        vmax : float
            Maximum value of the colorbar.
        """

        super().__init__(w_traj, frame, box_size, centre,
            arrow_width, arrow_head_width, arrow_head_length)   # initialise superclass

        self.d2min = w_traj.d2min(frame, frame + dt, *self.particles)   # particles' nonaffine squared displacement between frame and frame + dt

        self.vmin = np.log10(np.min(self.d2min[self.d2min != 0]))
        self.vmax = np.log10(np.max(self.d2min))
        try:
            self.vmin = np.log10(kwargs['vmin'])
        except (KeyError, AttributeError): pass # 'vmin' not in keyword arguments or None
        try:
            self.vmax = np.log10(kwargs['vmax'])
        except (KeyError, AttributeError): pass # 'vmax' not in keyword arguments or None

        self.colorbar(self.vmin, self.vmax, cmap=plt.cm.Greys)  # add colorbar to figure
        self.colormap.set_label(r'$\log D^2_{min}$',
            labelpad=pad, rotation=270)                         # colorbar legend

        self.draw()

    def draw(self):
        """
        Plots figure.
        """

        for particle, d2min in zip(self.particles, self.d2min): # for particle and particle's nonaffine squared displacement in rendered box
            self.draw_circle(particle, color=self.scalarMap.to_rgba(
                np.log10(np.linalg.norm(d2min))), fill=True)    # draw particle circle with color corresponding to nonaffine square displacement

class Velocity(_Frame):
    """
    Plotting class specific to 'velocity' mode.
    """

    def __init__(self, u_traj, w_traj, frame, box_size, centre, arrow_width,
        arrow_head_width, arrow_head_length, pad, **kwargs):
        """
        Initialises and plots figure.

        Parameters
        ----------
        u_traj : active_particles.dat.Dat
    		Unwrapped trajectory object.
        w_traj : active_particles.dat.Gsd
    		Wrapped trajectory object.
        frame : int
            Frame to render.
        box_size : float
            Length of the square box to render.
        centre : 2-uple like
            Centre of the box to render.
        arrow_width : float
            Width of the arrows.
        arrow_head_width : float
            Width of the arrows' head.
        arrow_head_length : float
            Length of the arrows' head.
        pad : float
            Separation between label and colormap.

        Optional keyword parameters
        ---------------------------
        vmin : float
            Minimum value of the colorbar.
        vmax : float
            Maximum value of the colorbar.
        """

        super().__init__(w_traj, frame, box_size, centre,
            arrow_width, arrow_head_width, arrow_head_length)   # initialise superclass

        self.velocities = u_traj.velocity(frame, *self.particles)   # particles' velocities at frame frame

        self.vmin, self.vmax = amplogwidth(self.velocities)
        try:
            self.vmin = np.log10(kwargs['vmin'])
        except (KeyError, AttributeError): pass # 'vmin' not in keyword arguments or None
        try:
            self.vmax = np.log10(kwargs['vmax'])
        except (KeyError, AttributeError): pass # 'vmax' not in keyword arguments or None


        self.colorbar(self.vmin, self.vmax) # add colorbar to figure
        self.colormap.set_label(r'$\log||\vec{v}(t)||$',
            labelpad=pad, rotation=270)     # colorbar legend

        self.draw()

    def draw(self):
        """
        Plots figure.
        """

        for particle, velocity in zip(self.particles, self.velocities): # for particle and particle's velocity in rendered box
            self.draw_circle(particle, color=self.scalarMap.to_rgba(
                np.log10(np.linalg.norm(velocity))), fill=True)         # draw particle circle with color corresponding to velocity
            self.draw_arrow(particle, *normalise1D(velocity)
                *0.75*self.diameters[particle])                         # draw velocity direction arrow

class Trajectory(_Frame):
    """
    Plotting class specific to 'trajectory' mode.
    """

    def __init__(self, u_traj, w_traj, frame, box_size, centre, arrow_width,
        arrow_head_width, arrow_head_length, dt=0, **kwargs):
        """
        Initialises and plots figure.

        Parameters
        ----------
        u_traj : active_particles.dat.Dat
    		Unwrapped trajectory object.
        w_traj : active_particles.dat.Gsd
    		Wrapped trajectory object.
        frame : int
            Frame to render.
        box_size : float
            Length of the square box to render.
        centre : 2-uple like
            Centre of the box to render.
        arrow_width : float
            Width of the arrows.
        arrow_head_width : float
            Width of the arrows' head.
        arrow_head_length : float
            Length of the arrows' head.
        dt : int
            Lag time for displacement. (default: 0)
        """

        super().__init__(w_traj, frame, box_size, centre,
            arrow_width, arrow_head_width, arrow_head_length)   # initialise superclass

        global trajectory_tracer_particle                               # index of tracer particle
        try:
            if not(trajectory_tracer_particle in self.particles):       # tracer particle not in frame
                raise NameError                                         # do as if tracer particle were not defined
        except NameError:                                               # tracer particle not defined
            if get_env('TRACER_PARTICLE', default=True, vartype=bool):  # TRACER_PARTICLE mode
                trajectory_tracer_particle = np.argmin(
                    np.sum(self.positions**2, axis=-1))                 # tracer particle at centre of frame
            else:
                trajectory_tracer_particle = -1                         # there is no particle with index -1

        self.displacements = u_traj.displacement(frame, frame + dt,
            *self.particles)   # particles' displacements between time and time + dt

        self.draw()

    def draw_circle(self, particle, color='black', fill=False):
        """
        Draws circle at particle's position with particle's diameter.

        Tracer particle is filled.

        Parameters
        ----------
        particle : int
            Particle index.
        color : any matplotlib color
            Circle color. (default: 'black')
        fill : bool
            Filling the circle with same color. (default: False)
        """

        if particle == trajectory_tracer_particle:
            super().draw_circle(particle, color=color, fill=True)   # fill circle if drawn particle is tracer particle
        else:
            super().draw_circle(particle, color=color, fill=fill)

    def draw(self):
        """
        Plots figure.
        """

        for particle, displacement in zip(self.particles,
            self.displacements):                        # for particle and particle's displacement in rendered box
            self.draw_circle(particle)                  # draw particle circle
            self.draw_arrow(particle, *displacement)    # draw particle dispalcement between frame and frame + dt

class Displacement(_Frame):
    """
    Plotting class specific to 'displacement' mode.
    """

    def __init__(self, u_traj, w_traj, frame, box_size, centre, arrow_width,
        arrow_head_width, arrow_head_length, pad, dt=0, **kwargs):
        """
        Initialises and plots figure.

        Parameters
        ----------
        u_traj : active_particles.dat.Dat
    		Unwrapped trajectory object.
        w_traj : active_particles.dat.Gsd
    		Wrapped trajectory object.
        frame : int
            Frame to render.
        box_size : float
            Length of the square box to render.
        centre : 2-uple like
            Centre of the box to render.
        arrow_width : float
            Width of the arrows.
        arrow_head_width : float
            Width of the arrows' head.
        arrow_head_length : float
            Length of the arrows' head.
        pad : float
            Separation between label and colormap.
        dt : int
            Lag time for displacement. (default=0)

        Optional keyword parameters
        ---------------------------
        vmin : float
            Minimum value of the colorbar.
        vmax : float
            Maximum value of the colorbar.
        """

        super().__init__(w_traj, frame, box_size, centre,
            arrow_width, arrow_head_width, arrow_head_length)   # initialise superclass

        self.displacements = u_traj.displacement(frame, frame + dt,
            *self.particles)    # particles' displacements between frame and frame + dt

        self.vmin, self.vmax = amplogwidth(self.displacements)
        try:
            self.vmin = np.log10(kwargs['vmin'])
        except (KeyError, AttributeError): pass # 'vmin' not in keyword arguments or None
        try:
            self.vmax = np.log10(kwargs['vmax'])
        except (KeyError, AttributeError): pass # 'vmax' not in keyword arguments or None

        self.colorbar(self.vmin, self.vmax) # add colorbar to figure
        self.colormap.set_label(r'$\log||\vec{u}(t, t+\Delta t)||$',
            labelpad=pad, rotation=270)     # colorbar legend

        self.draw()

    def draw(self):
        """
        Plots figure.
        """

        for particle, displacement in zip(self.particles,
            self.displacements):                                    # for particle and particle's displacement in rendered box
            self.draw_circle(particle, color=self.scalarMap.to_rgba(
                np.log10(np.linalg.norm(displacement))), fill=True) # draw particle circle with color corresponding to displacement amplitude
            self.draw_arrow(particle, *normalise1D(displacement)
                *0.75*self.diameters[particle])                     # draw displacement direction arrow

# SCRIPT

if __name__ == '__main__':  # executing as script

    startTime = datetime.now()

    # VARIABLE DEFINITIONS

    mode = get_env('MODE', default='displacement')          # plotting mode
    if mode == 'd2min':
        plotting_object = D2min
        naming_standard = naming.D2min()
    elif mode == 'velocity':
        plotting_object = Velocity
        naming_standard = naming.Velocity()
    elif mode == 'trajectory':
        plotting_object = Trajectory
        naming_standard = naming.Trajectory()
    elif mode == 'displacement':
        plotting_object = Displacement
        naming_standard = naming.Displacement()
    else: raise ValueError('Mode %s is not known.' % mode)  # mode is not known

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    wrap_file_name = get_env('WRAPPED_FILE',
        default=joinpath(data_dir, naming.wrapped_trajectory_file))	    # wrapped trajectory file (.gsd)
    unwrap_file_name = get_env('UNWRAPPED_FILE',
        default=joinpath(data_dir, naming.unwrapped_trajectory_file))	# unwrapped trajectory file (.dat)

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)  # initial frame to render
    dt = get_env('DT', default=-1, vartype=int)                     # displacement lag time

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    prep_frames = ceil(parameters['prep_steps']/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

    box_size = get_env('BOX_SIZE', default=parameters['box_size'],
		vartype=float)									# size of the square box to consider
    centre = (get_env('X_ZERO', default=0, vartype=float),
        get_env('Y_ZERO', default=0, vartype=float))	# centre of the box

    Nentries = parameters['N_steps']//parameters['period_dump'] # number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame

    # NAMING

    attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'init_frame': init_frame, 'dt': dt,
        'box_size': box_size, 'x_zero': centre[0], 'y_zero': centre[1]} # attributes displayed in filenames

    # STANDARD OUTPUT

    if 'SLURM_JOB_ID' in envvar:	# script executed from Slurm job scheduler
        slurm_output(joinpath(data_dir, 'out'), naming_standard, attributes)

    # FIGURE PARAMETERS

    vmin = get_env('V_MIN', vartype=float) # minimum value of the colorbar
    vmax = get_env('V_MAX', vartype=float) # maximum value of the colorbar

    frame_hor = get_env('FRAME_HORIZONTAL_SIZE', default=_frame_hor,
        vartype=float)  # horizontal size of the frame (in inches)
    frame_ver = get_env('FRAME_VERTICAL_SIZE', default=_frame_ver,
        vartype=float)  # vertical size of the frame (in inches)
    mpl.rcParams['figure.figsize'] = (frame_hor, frame_ver)

    frame_def = get_env('FRAME_DEFINITION', default=_frame_def,
        vartype=float)                                                  # definition of image (in dots per inches (dpi))
    font_size = get_env('FONT_SIZE', default=_font_size, vartype=float) # font size
    mpl.rcParams.update({'savefig.dpi': frame_def, 'font.size': font_size})

    arrow_width = get_env('ARROW_WIDTH', default=_arrow_width,
        vartype=float)  # width of the arrows
    arrow_head_width = get_env('HEAD_WIDTH', default=_arrow_head_width,
        vartype=float)  # width of the arrows' head
    arrow_head_length = get_env('HEAD_LENGTH', default=_arrow_head_length,
        vartype=float)  # length of the arrows' head

    pad = get_env('PAD', default=_colormap_label_pad, vartype=float)    # separation between label and colormap

    # LEGEND SUPTITLE

    display_suptitle = get_env('SUPTITLE', default=True, vartype=bool)  # display suptitle

    def suptitle(frame):
        """
        Returns figure suptitle.

        NOTE: Returns empty string if display_suptitle=False.

        Parameters
        ----------
        frame : int
            Index of rendered frame.
        """

        if not(display_suptitle): return ''

        suptitle = str(r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
    		% (parameters['N'], parameters['density'], parameters['vzero'],
    		parameters['dr']))
        suptitle += str(r'$, L=%.3e$' % parameters['box_size'])
        if 'BOX_SIZE' in envvar: suptitle += str(r'$, L_{new}=%.3e$' % box_size)
        suptitle += '\n'
        if 'X_ZERO' in envvar or 'Y_ZERO' in envvar:
            suptitle += str(r'$x_0 = %.3e, y_0 = %.3e$' % centre) + '\n'
        suptitle += str(r'$t = %.5e$'
            % (frame*parameters['period_dump']*parameters['time_step']))
        if mode == 'trajectory' or mode == 'displacement' or mode == 'd2min':
            suptitle += str(r'$, \Delta t = %.5e$'
                % (dt*parameters['period_dump']*parameters['time_step']))

        return suptitle

    # MODE SELECTION

    if get_env('PLOT', default=False, vartype=bool):    # PLOT mode

        Nframes = Nentries - init_frame  # number of frames available for the calculation
        dt = Nframes + dt if dt < 0 else dt

        with open(wrap_file_name, 'rb') as wrap_file,\
            open(unwrap_file_name, 'rb') as unwrap_file:    # opens wrapped and unwrapped trajectory files

            w_traj = Gsd(wrap_file, prep_frames=prep_frames)    # wrapped trajectory object
            u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object

            figure = plotting_object(u_traj, w_traj, init_frame, box_size,
                centre, arrow_width, arrow_head_width, arrow_head_length,
                pad=pad, dt=dt, vmin=vmin, vmax=vmax)
            figure.fig.suptitle(suptitle(init_frame))

            if get_env('SAVE', default=False, vartype=bool):    # SAVE mode
                figure_name, = naming_standard.filename(**attributes)
                figure.fig.savefig(joinpath(data_dir,
                    get_env('FIGURE_NAME', default=figure_name)))

    if get_env('MOVIE', default=False, vartype=bool):   # MOVIE mode

        frame_fin = get_env('FINAL_FRAME', default=Nentries, vartype=int)       # final movie frame
        frame_per = get_env('FRAME_PERIOD', default=_frame_per, vartype=int)    # frame rendering period
        frame_max = get_env('FRAME_MAXIMUM', default=_frame_max, vartype=int)   # maximum number of frames

        attributes = {**attributes,
            **{'frame_fin': frame_fin, 'frame_per': frame_per,
            'frame_max': frame_max}}
        movie_dir = joinpath(data_dir,
            naming_standard.movie(folder=True).filename(**attributes)[0])   # movie directory name
        mkdir(movie_dir)                                                    # create movie directory
        mkdir(joinpath(movie_dir, 'frames'), replace=True)                  # create frames directory (or replaces it if existing)

        Nframes = np.min([Nentries, frame_fin]) - init_frame                    # number of frames available for the movie
        Ntimes = Nframes//frame_per                                             # maximum number of rendered frames
        frames = np.array(list(OrderedDict.fromkeys(map(int, init_frame +
            np.linspace(0, Nframes - 1, Ntimes, endpoint=False, dtype=int)))))  # usable frames

        if dt < 0: dt = np.min(np.abs(frames - np.roll(frames, shift=1)))

        frames = frames[frames + dt < Nentries - 1][:frame_max] # rendered frames

        with open(wrap_file_name, 'rb') as wrap_file,\
            open(unwrap_file_name, 'rb') as unwrap_file:    # opens wrapped and unwrapped trajectory files

            w_traj = Gsd(wrap_file, prep_frames=prep_frames)    # wrapped trajectory object
            u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object

            for frame in frames:    # for rendered frames
                sys.stdout.write(
                    'Frame: %d' % (frames.tolist().index(frame) + 1)
                    + "/%d \r" % len(frames))

                figure = plotting_object(u_traj, w_traj, frame, box_size,
                    centre, arrow_width, arrow_head_width, arrow_head_length,
                    dt=dt, vmin=vmin, vmax=vmax)                        # plot frame
                figure.fig.suptitle(suptitle(frame))
                figure.fig.savefig(joinpath(movie_dir, 'frames',
                    '%010d' % frames.tolist().index(frame) + '.png'))   # save frame
                del figure                                              # delete (close) figure

        subprocess.call([
            'ffmpeg', '-r', '5', '-f', 'image2', '-s', '1280x960', '-i',
            joinpath(movie_dir , 'frames', '%10d.png'),
            '-pix_fmt', 'yuv420p', '-y',
            joinpath(movie_dir,
            naming_standard.movie().filename(**attributes)[0])
            ])  # generate movie from frames

    # EXECUTION TIME
    print("Execution time: %s" % (datetime.now() - startTime))

    if get_env('SHOW', default=False, vartype=bool):    # SHOW mode
        plt.show()
