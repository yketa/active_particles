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
    | 'velocity'     | Instantaneous velocity | Relative to  | Amplitude of   |
    |                | direction              | diameter     | velocity       |
    |________________|________________________|______________|________________|
    | 'trajectory'   | Displacement direction | Displacement | No             |
    |                |                        | amplitude    |                |
    |________________|________________________|______________|________________|
    | 'displacement' | Displacement direction | Relative to  | Amplitude of   |
    |                |                        | diameter     | displacement   |
    |________________|________________________|______________|________________|
    DEFAULT: velocity
SHOW : bool
	Show figure.
	DEFAULT: False
SAVE : bool
	Save figure.
	DEFAULT: False

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
FRAME : int
	Frame to render.
	NOTE: FRAME < 0 will be interpreted as the frame to render being the middle
	frame of the simulation.
	DEFAULT: -1
DT : int
	Lag time for displacement.
	NOTE: TIME < 0 will be interpreted as a lag time corresponding to the total
	number of simulation frames - FRAME + TIME.
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
    DEFAULT: ['velocity' mode] 1e-2 times the self-propulsion velocity
             [other modes] 1e-2 times the self-propulsion velocity times the
                           lag time
    NOTE: Colorbar is represented in logarithmic scale so V_MIN > 0.
V_MAX : float
    Maximum value of the colorbar.
    DEFAULT: ['velocity' mode] self-propulsion velocity
             [other modes] self-propulsion velocity times the lag time
    NOTE: Colorbar is represented in logarithmic scale so V_MAX > 0.
ARROW_WIDTH : float
    Width of the arrows.
    DEFAULT: active_articles.analysis.frame._arrow_width
HEAD_WIDTH : float
    Width of the arrows' head.
    DEFAULT: active_articles.analysis.frame._arrow_head_width
HEAD_LENGTH : float
    Length of the arrows' head.
    DEFAULT: active_articles.analysis.frame._arrow_head_length
FRAME_VERTICAL_SIZE : float
    Vertical size of the frame (in inches).
    DEFAULT: active_articles.analysis.frame._frame_ver
FRAME_HORIZONTAL_SIZE : float
    Horizontal size of the frame (in inches).
    DEFAULT: active_articles.analysis.frame._frame_hor
FRAME_DEFINITION [SAVE mode] : float
    Definition of image (in dots per inches (dpi)).
    DEFAULT: active_articles.analysis.frame._frame_def

Output
------
> Prints execution time.
> Plots system according to plotting mode.
[SHOW mode]
> Displays figure.
[SAVE mode]
> Saves figure in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env
from active_particles.dat import Dat, Gsd
from active_particles.maths import normalise1D

from os import getcwd
from os import environ as envvar

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
        self.ax.arrow(*self.positions[particle], dx, dy, color=color,
            width=length*self.arrow_width,
            head_width=length*self.arrow_head_width,
            head_length=length*self.arrow_head_length, zorder=1)

    def colorbar(self, vmin, vmax):
        """
        Adds colorbar to plot.

        Parameters
        ----------
        vmin : float
            Minimum value of the colorbar.
        vmax : float
            Maximum value of the colorbar.
        """

        cmap = plt.cm.jet
        vNorm = ColorsNormalise(vmin=vmin, vmax=vmax)
        self.scalarMap = ScalarMappable(norm=vNorm, cmap=cmap)

        self.colormap_ax = make_axes_locatable(self.ax).append_axes('right',
            size='5%', pad=0.05)
        self.colormap = mpl.colorbar.ColorbarBase(self.colormap_ax, cmap=cmap,
            norm=vNorm, orientation='vertical')

class Velocity(_Frame):
    """
    Plotting class specific to 'velocity' mode.
    """

    def __init__(self, u_traj, w_traj, frame, box_size, centre, vmin, vmax,
        arrow_width, arrow_head_width, arrow_head_length):
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
        vmin : float
            Minimum value of the colorbar.
        vmax : float
            Maximum value of the colorbar.
        arrow_width : float
            Width of the arrows.
        arrow_head_width : float
            Width of the arrows' head.
        arrow_head_length : float
            Length of the arrows' head.
        """

        super().__init__(w_traj, frame, box_size, centre,
            arrow_width, arrow_head_width, arrow_head_length)   # initialise superclass
        self.colorbar(vmin, vmax)                               # add colorbar to figure

        self.velocities = u_traj.velocity(frame)    # particles' velocities at frame frame

        self.draw()

    def draw(self):
        """
        Plots figure.
        """

        for particle in self.particles:         # for particle in rendered box
            self.draw_circle(particle, color=self.scalarMap.to_rgba(
            np.log10(np.linalg.norm(self.velocities[particle]))),
                fill=True)                      # draw particle circle with color corresponding to velocity
            self.draw_arrow(particle,
                *normalise1D(self.velocities[particle])
                *0.75*self.diameters[particle]) # draw velocity direction arrow

class Trajectory(_Frame):
    """
    Plotting class specific to 'trajectory' mode.
    """

    def __init__(self, u_traj, w_traj, frame, dt, box_size, centre,
        arrow_width, arrow_head_width, arrow_head_length):
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
        dt : int
            Lag time for displacement.
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

        super().__init__(w_traj, frame, box_size, centre,
            arrow_width, arrow_head_width, arrow_head_length)   # initialise superclass

        self.displacements = u_traj.displacement(frame, frame + dt)   # particles' displacements between time and time + dt

        self.draw()

    def draw(self):
        """
        Plots figure.
        """

        for particle in self.particles:                                 # for particle in rendered box
            self.draw_circle(particle)                                  # draw particle circle
            self.draw_arrow(particle, *self.displacements[particle])    # draw particle dispalcement between frame and frame + dt

class Displacement(_Frame):
    """
    Plotting class specific to 'displacement' mode.
    """

    def __init__(self, u_traj, w_traj, frame, dt, box_size, centre, vmin, vmax,
        arrow_width, arrow_head_width, arrow_head_length):
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
        dt : int
            Lag time for displacement.
        box_size : float
            Length of the square box to render.
        centre : 2-uple like
            Centre of the box to render.
        vmin : float
            Minimum value of the colorbar.
        vmax : float
            Maximum value of the colorbar.
        arrow_width : float
            Width of the arrows.
        arrow_head_width : float
            Width of the arrows' head.
        arrow_head_length : float
            Length of the arrows' head.
        """

        super().__init__(w_traj, frame, box_size, centre,
            arrow_width, arrow_head_width, arrow_head_length)   # initialise superclass
        self.colorbar(vmin, vmax)                               # add colorbar to figure

        self.displacements = u_traj.displacement(frame, frame + dt)   # particles' displacements between frame and frame + dt

        self.draw()

    def draw(self):
        """
        Plots figure.
        """

        for particle in self.particles:         # for particle in rendered box
            self.draw_circle(particle, color=self.scalarMap.to_rgba(
            np.log10(np.linalg.norm(self.displacements[particle]))),
            fill=True)                          # draw particle circle with color corresponding to displacement amplitude
            self.draw_arrow(particle,
                *normalise1D(self.displacements[particle])
                *0.75*self.diameters[particle]) # draw displacement direction arrow

# DEFAULT VARIABLES

_frame_ver = 12 # default vertical size of the frames (in inches)
_frame_hor = 16 # default horizontal size of the frames (in inches)
_frame_def = 80 # default definition of images (in dots per inches (dpi))

_arrow_width = 1e-3                         # default width of the arrows
_arrow_head_width = _arrow_width*3e2        # default width of the arrows' head
_arrow_head_length = _arrow_head_width*1.5  # default length of the arrows' head

if __name__ == '__main__':  # executing as script

    startTime = datetime.now()

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    wrap_file_name = get_env('WRAPPED_FILE',
        default=data_dir + '/' + naming.wrapped_trajectory_file)	# wrapped trajectory file (.gsd)
    unwrap_file_name = get_env('UNWRAPPED_FILE',
        default=data_dir + '/' + naming.unwrapped_trajectory_file)	# unwrapped trajectory file (.dat)

    frame = get_env('FRAME', default=-1, vartype=int)   # frame to render
    dt = get_env('DT', default=-1, vartype=int)         # displacement lag time

    parameters_file = get_env('PARAMETERS_FILE',
		default=data_dir + '/' + naming.parameters_file)	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
        parameters = pickle.load(param_file)				# parameters hash table

    prep_frames = ceil(parameters['prep_steps']/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

    box_size = get_env('BOX_SIZE', default=parameters['box_size'],
		vartype=float)									# size of the square box to consider
    centre = (get_env('X_ZERO', default=0, vartype=float),
        get_env('Y_ZERO', default=0, vartype=float))	# centre of the box

    Nentries = parameters['N_steps']//parameters['period_dump'] # number of time snapshots in unwrapped trajectory file
    frame = int(Nentries/2) if frame < 0 else frame
    Nframes = Nentries - frame                                  # number of frames available for the calculation
    dt = Nframes + dt if dt < 0 else dt

    # FIGURE PARAMETERS

    frame_ver = get_env('FRAME_VERTICAL_SIZE', default=_frame_ver,
        vartype=float)  # vertical size of the frame (in inches)
    frame_hor = get_env('FRAME_HORIZONTAL_SIZE', default=_frame_hor,
        vartype=float)  # horizontal size of the frame (in inches)
    frame_def = get_env('FRAME_DEFINITION', default=_frame_def,
        vartype=float)  # definition of image (in dots per inches (dpi))

    arrow_width = get_env('ARROW_WIDTH', default=_arrow_width,
        vartype=float)  # width of the arrows
    arrow_head_width = get_env('HEAD_WIDTH', default=_arrow_head_width,
        vartype=float)  # width of the arrows' head
    arrow_head_length = get_env('HEAD_LENGTH', default=_arrow_head_length,
        vartype=float)  # length of the arrows' head

    # NAMING

    attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'frame': frame, 'dt': dt, 'box_size': box_size,
        'x_zero': centre[0], 'y_zero': centre[1]}   # attributes displayed in filenames

    # MODE SELECTION

    mode = get_env('MODE', default='velocity')  # plotting mode

    with open(wrap_file_name, 'rb') as wrap_file,\
        open(unwrap_file_name, 'rb') as unwrap_file:    # opens wrapped and unwrapped trajectory files

        w_traj = Gsd(wrap_file, prep_frames=prep_frames)    # wrapped trajectory object
        u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object

        if mode == 'velocity':

            vmin = np.log10(get_env('V_MIN', default=1e-2*parameters['vzero'],
                vartype=float)) # minimum velocity for the colorbar
            vmax = np.log10(get_env('V_MAX', default=parameters['vzero'],
                vartype=float)) # maximum velocity for the colorbar

            figure = Velocity(u_traj, w_traj, frame, box_size, centre, vmin,
                vmax, arrow_width, arrow_head_width, arrow_head_length)
            figure.colormap.set_label(r'$\log||\vec{v}(t)||$',
                labelpad=30, rotation=270)

            naming_standard = naming.Velocity()

        elif mode == 'trajectory':

            figure = Trajectory(u_traj, w_traj, frame, dt, box_size, centre,
                arrow_width, arrow_head_width, arrow_head_length)

            naming_standard = naming.Trajectory()

        elif mode == 'displacement':

            vmin = np.log10(get_env('V_MIN',
                default=1e-2*parameters['vzero']*
                dt*parameters['period_dump']*parameters['time_step'],
                vartype=float)) # minimum velocity for the colorbar
            vmax = np.log10(get_env('V_MAX',
                default=parameters['vzero']*
                dt*parameters['period_dump']*parameters['time_step'],
                vartype=float)) # maximum velocity for the colorbar

            figure = Displacement(u_traj, w_traj, frame, dt, box_size, centre,
                vmin, vmax, arrow_width, arrow_head_width, arrow_head_length)
            figure.colormap.set_label(r'$\log||\vec{u}(t, t+\Delta t)||$',
                labelpad=30, rotation=270)

            naming_standard = naming.Displacement()

        else: raise ValueError('Mode %s is not known.' % mode)  # mode is not known

    # FIGURE DIMENSIONS
    figure.fig.set_size_inches(frame_hor, frame_ver)

    # LEGEND SUPTITLE
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
    if mode == 'trajectory' or mode == 'displacement':
        suptitle += str(r'$, \Delta t = %.5e$'
            % (dt*parameters['period_dump']*parameters['time_step']))
    figure.fig.suptitle(suptitle)

    # EXECUTION TIME
    print("Execution time: %s" % (datetime.now() - startTime))

    if get_env('SAVE', default=False, vartype=bool):    # SAVE mode
        figure_name, = naming_standard.filename(**attributes)
        figure.fig.savefig(data_dir + '/' +
            get_env('FIGURE_NAME', default=figure_name), dpi=frame_def)

    if get_env('SHOW', default=False, vartype=bool):    # SHOW mode
        plt.show()
