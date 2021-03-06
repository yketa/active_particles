"""
Module cuu calculates or plots displacements and their correlations.

Files are saved according to active_particles.naming.Cuu (displacement
correlation), active_particles.naming.Cww (relative displacement correlation),
active_particles.naming.Cdd (displacement norm correlation),
active_particles.naming.Cee (displacement direction correlation) and
active_particles.naming.Cnn (density correlation).

A brief description of the algorithm can be found at:
https://yketa.github.io/UBC_2018_Wiki/#Displacement%20correlations

Environment modes
-----------------
COMPUTE : bool
	Compute shear strain and displacement vorticity.
	DEFAULT: False
PLOT : bool
	Plot saved shear strain and displacement vorticity as well as their
	correlations.
	DEFAULT: False
SHOW [COMPUTE or PLOT mode] : bool
	Show graphs.
	DEFAULT: False
SAVE [COMPUTE or PLOT mode] : bool
	Save graphs.
	DEFAULT: False
GRID_CIRCLE [SHOW mode] : bool
	Analyse graphically values of corrected correlations at fixed radius.
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
INITIAL_FRAME : int
	Frame to consider as initial.
	NOTE: INITIAL_FRAME < 0 will be interpreted as the initial frame being
	      the middle frame of the simulation.
	DEFAULT: -1
TIME : int
	Lag time for displacement.
	NOTE: TIME < 0 will be interpreted as a lag time corresponding to the total
	      number of simulation frames - INITIAL_FRAME + TIME.
	DEFAULT: -1
INTERVAL_MAXIMUM : int
	Maximum number of intervals of length dt considered in correlations
	calculations.
	DEFAULT: 1
N_CASES : int
	Number of boxes in each direction to compute the shear strain and
	displacement vorticity grid.
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
R_MIN [PLOT or SHOW mode] : float
	Minimum radius for correlations plots.
	DEFAULT: active_particles.analysis.cuu._r_min
R_MAX [PLOT or SHOW mode] : float
	Maximum radius for correlations plots.
	DEFAULT: active_particles.analysis.cuu._r_max
CUU_MIN [PLOT or SHOW mode] : float
	Minimum displacement correlation for correlation plots.
	DEFAULT: active_particles.analysis.cuu._Cuu_min
CUU_MAX [PLOT or SHOW mode] : float
	Maximum displacement correlation for correlation plots.
	DEFAULT: active_particles.analysis.cuu._Cuu_max
CWW_MIN [PLOT or SHOW mode] : float
	Minimum relative displacement correlation for correlation plots.
	DEFAULT: active_particles.analysis.cuu._Cww_min
CWW_MAX [PLOT or SHOW mode] : float
	Maximum relative displacement correlation for correlation plots.
	DEFAULT: active_particles.analysis.cuu._Cww_max
CDD_MIN [PLOT or SHOW mode] : float
	Minimum displacement norm correlation for correlation plots.
	DEFAULT: active_particles.analysis.cuu._Cdd_min
CDD_MAX [PLOT or SHOW mode] : float
	Maximum displacement norm correlation for correlation plots.
	DEFAULT: active_particles.analysis.cuu._Cdd_max
CEE_MIN [PLOT or SHOW mode] : float
	Minimum displacement direction correlation for correlation plots.
	DEFAULT: active_particles.analysis.cuu._Cee_min
CEE_MAX [PLOT or SHOW mode] : float
	Maximum displacement direction correlation for correlation plots.
	DEFAULT: active_particles.analysis.cuu._Cee_max
AXIS [PLOT or SHOW mode] : string
	Axis scale for correlation plots.
	NOTE: 'LINLIN', 'LOGLIN', 'LINLOG' or 'LOGLOG'.
	DEFAULT: 'LOGLOG'

Output
------
[COMPUTE MODE]
> Prints execution time.
> Saves 2D and 1D density correlations according to active_particles.naming.Cnn
standards in DATA_DIRECTORY.
> Saves 2D, 1D, longitudinal and transversal displacement correlations and
1D correlations corrected with density correlations according to
active_particles.naming.Cuu standards in DATA_DIRECTORY.
> Saves 2D, 1D, longitudinal and transversal relative displacement correlations
and 1D correlations corrected with density correlations according to
active_particles.naming.Cww standards in DATA_DIRECTORY.
> Saves 2D and 1D displacement norm correlations and 1D correlations corrected
with density correlations according to active_particles.naming.Cdd standards in
DATA_DIRECTORY.
> Saves 2D, 1D, longitudinal and transversal displacement norm correlations and
1D correlations corrected with density correlations according to
active_particles.naming.Cee standards in DATA_DIRECTORY.
[SHOW or PLOT mode]
> Plots correlations for all variables.
[SAVE mode]
> Saves correlation figures in DATA_DIRECTORY.
"""

import active_particles.naming as naming

from active_particles.init import get_env, slurm_output
from active_particles.dat import Dat, Gsd
from active_particles.maths import relative_positions, wo_mean, g2Dto1Dsquare

from active_particles.analysis.correlations import corField2D_scalar_average,\
    corField2D_vector_average_Cnn, CorGrid

from os import getcwd
from os import environ as envvar
from os.path import join as joinpath

from math import ceil

import numpy as np

import pickle

from operator import itemgetter

from collections import OrderedDict

from datetime import datetime

import matplotlib as mpl
if not(get_env('SHOW', default=False, vartype=bool)):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable

from active_particles.plot.mpl_tools import GridCircle

# DEFAULT VARIABLES

_r_min = 1  # default minimum radius for correlations plots
_r_max = 20	# default maximum radius for correlations plots

_Cuu_min = 1e-3 # default minimum displacement correlation for correlation plots
_Cuu_max = 1    # default maximum displacement correlation for correlation plots

_Cww_min = 1e-3 # default minimum relative displacement correlation for correlation plots
_Cww_max = 1    # default maximum relative displacement correlation for correlation plots

_Cdd_min = 1e-1 # default minimum displacement norm correlation for correlation plots
_Cdd_max = 2    # default maximum displacement norm correlation for correlation plots

_Cee_min = 1e-3 # default minimum displacement direction correlation for correlation plots
_Cee_max = 1    # default maximum displacement direction correlation for correlation plots

# FUNCTIONS AND CLASSES

def displacement_grid(box_size, centre, Ncases, time, dt, w_traj, u_traj):
    """
    Calculates displcament grid from square uniform coarse-graining.

    Parameters
	----------
	box_size : float
		Length of the considered system's square box.
    centre : float array
        Centre of the box.
	Ncases : int
		Number of boxes in each direction to compute the displacements.
	time : int
		Frame at which displacements will be calculated.
	dt : int
		Length of the interval of time for which the displacements are
		calculated.
	w_traj : active_particles.dat.Gsd
		Wrapped trajectory object.
	u_traj : active_particles.dat.Dat
		Unwrapped trajectory object.

    Returns
    -------
    ugrid : 2D array like
        Displacement grid.
    """

    return w_traj.to_grid(
        time + dt*get_env('ENDPOINT', default=False, vartype=bool),
        u_traj.displacement(time, time + dt),
        Ncases=Ncases, box_size=box_size, centre=centre)

def displacement_related_grids(box_size, centre, Ncases, time,
	dt, w_traj, u_traj):
	"""
	Calculates grids of displacement (from
	active_particles.analysis.cuu.displacement_grid), density, relative
	displacement, displacement norm and displacement direction grids.

    Parameters
	----------
	box_size : float
		Length of the considered system's square box.
    centre : float array
        Centre of the box.
	Ncases : int
		Number of boxes in each direction to compute the displacements.
	time : int
		Frame at which displacements will be calculated.
	dt : int
		Length of the interval of time for which the displacements are
		calculated.
	w_traj : active_particles.dat.Gsd
		Wrapped trajectory object.
	u_traj : active_particles.dat.Dat
		Unwrapped trajectory object.

    Returns
    -------
    ddgrid : 2D array like
        Displacement norm grid concatenated with itself.
        NOTE: These are concatenated for dimension reasons.
    ugrid : 2D array like
        Displacement grid.
    wgrid : 2D array like
        Relative displacement grid.
    egrid : 2D array like
        Displacement direction grid.

    Output
	------
	Prints neighbours grid computation time.
	"""

	ugrid = displacement_grid(box_size, centre, Ncases, time, dt,
		w_traj, u_traj)	# displacement grid

	wgrid = ugrid - np.mean(ugrid, axis=(0, 1)) # relative displacement grid

	dgrid = np.sqrt(np.sum(ugrid**2, axis=-1))  # displacement norm grid
	dgridr = np.reshape(dgrid, dgrid.shape + (1,))

	egrid = np.divide(ugrid, dgridr, out=np.zeros(ugrid.shape),
        where=dgridr!=0) # displacement direction

	return np.concatenate((dgridr, dgridr), axis=-1), ugrid, wgrid, egrid  # displacement variable grid

class Cnn:
	"""
	Manipulates density self-correlations computed from displacement grids.
	"""

	def __init__(self, ugrid, box_size):
		"""
		Calculates grids of density and their averaged correlations.

		Parameters
		----------
		ugrid : array-like
			Array of displacements or list of array of displacements.
		box_size : float
			Length of the system's square box.
		"""

		self.ugrid = np.array(ugrid, ndmin=4)	# list of array of displacements
		self.box_size = box_size

		self.ngrid = (self.ugrid != 0).any(axis=-1)*1			# density grid
		self.cnn2D = corField2D_scalar_average(self.ngrid)		# 2D density correlation
		self.cnn1D = g2Dto1Dsquare(self.cnn2D, self.box_size)	# 1D averaged density correlation

	def save(self, attributes, dir=getcwd()):
		"""
		Saves 2D and 1D density correlation grids to directory dir.

		Parameters
		----------
		attributes : hash-table
			Attributes displayed in filename.
		dir : string
			Saving directory. (default: current working directory)
		"""

		self.filename, =  naming.Cnn().filename(**attributes)	# density correlation file name

		with open(joinpath(dir, self.filename), 'wb') as dump_file:
			pickle.dump([self.cnn2D, self.cnn1D], dump_file)

def c2Dtochi(c2D, box_size, r_min=None, r_max=None):
	"""
	For the 2D correlation grid c2D, this function returns the susceptibility
	of the square box system of length L, defined as
		chi = 1/L^2 \int dx dy c2D(x, y).

	Parameters
	----------
	c2D : 2D array
		2D correlation grid.
	box_size : float
		System square box size.
	r_min : float
		Minimum radius to consider in correlation integration.
		NOTE: if r_min == None, no minimum is considered.
		DEFAULT: None
	r_max : float
		Maximum radius to consider in correlation integration.
		NOTE: if r_max == None, r_max is considered to be box_size/2.
		DEFAULT: None

	Returns
	-------
	chi : float
		Susceptibility.
	"""

	c2Dgrid = CorGrid(c2D, box_size)
	if r_max == None: r_max = box_size/2

	c2Dr = np.sqrt(np.sum(
		c2Dgrid.display_grid.get_grid_coordinates()**2,
		axis=-1))

	c = c2Dgrid.display_grid[
		(r_min == None or c2Dr >= r_min) & (c2Dr <= r_max)]
	return np.sum(c)/np.prod(c2D.shape)

def c1Dtochi(c1D, box_size, r_min=None, r_max=None):
	"""
	For the cylindrically-averaged 1D correlation function c1D(r), this
	function returns the 2D susceptibility of the square box system of length
	L, defined as chi = 2\\pi/L^2 \\int_{r_min}^{r_max} dr c1D(r).
	NOTE: for displacement-related correlations, this susceptibility
	corresponds to the cooperativity.

	Parameters
	----------
	c1D : 1D array
        Correlation 1D average.
        NOTE: This has to be of the form (r, c1D(r)) with c1D(r) the averaged
        2D grid at radius r.
	box_size : float
		Simulation box size.
	r_min : float
		Lower bound of cooperativity integral. (default: None)
		NOTE: r_min=None corresponds to r_min=min(c1D[:, 0])
	r_max : float
		Higher bound of cooperativity integral. (default: None)
		NOTE: r_max=None corresponds to r_max=max(c1D[:, 0])

	Returns
	-------
	chi : float
		Susceptibility.
	"""

	c1D = np.array(sorted(c1D, key=lambda el: el[0]))

	if r_min == None: r_min = np.min(c1D[:, 0])
	if r_max == None: r_max = np.max(c1D[:, 0])

	c1Drmin = np.interp(r_min, c1D[:, 0], c1D[:, 1])	# value of c1D at r_min by linear interpolation
	c1Drmax = np.interp(r_max, c1D[:, 0], c1D[:, 1])	# value of c1D at r_max by linear interpolation

	r, c = np.transpose([[r, c] for r, c in c1D if r >= r_min and r <= r_max])	# values of radii and c1D at these radii in integration interval
	r = np.array([r_min, *r, r_max])											# radii with r_min and r_max
	c = np.array([c1Drmin, *c, c1Drmax])										# c1D values with ones at r_min and r_max

	return np.trapz(2*np.pi*r*c, r)/(box_size**2)

def plot_correlation(C, C2D, C1D, C1Dcor, C_min, C_max, naming_standard,
    **directional_correlations):
    """
    Plot correlations.

    Parameters
    ----------
    C : string
        Correlation name.
    C2D : 2D array
        Correlation 2D grid.
    C1D : 1D array
        Correlation 1D average.
        NOTE: This has to be of the form (r, C1D(r)) with C1D(r) the averaged
        2D grid at radius r.
    C1Dcor : 1D array
        Correlation 1D average, correction with density correlation.
        NOTE: This has to be of the form (r, C1D(r)) with C1D(r) the averaged
        2D grid at radius r.
    C_min : float
        Correlation minimum for plot.
    C_max : float
        Correlation maximum for plot.
    naming_standard : active_particles.naming standard
		Standard naming object.

    Optional keyword arguments
    --------------------------
    CL : float
        Longitudinal correlation.
    CT : float
        Transversal correlation.
    NOTE: These two variables have to be provided together.

	Returns
	-------
	fig : matplotlib figure
		Main figure.
	axs : array of matplotlib axis
		Main figure's axis.
	gc [GRID_CIRCLE mode] : active_particles.plot.mpl_tools.GridCircle object
		Grid circle object.
    """
    cmap = plt.cm.jet

    fig, axs = plt.subplots(2, 2)

    fig.set_size_inches(16, 16)
    fig.subplots_adjust(wspace=0.3)
    fig.subplots_adjust(hspace=0.3)

    suptitle = str(r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$'
        % (parameters['N'], parameters['density'], parameters['vzero'],
		parameters['dr']) + '\n' +
        r'$S_{init}=%.2e, \Delta t=%.2e$' % (init_frame,
		dt*parameters['period_dump']*parameters['time_step']) +
        r'$, S_{max}=%.2e, N_{cases}=%.2e$' % (int_max, Ncases))
    fig.suptitle(suptitle)

    # C2D

    cgrid = CorGrid(C2D, box_size, display_size=2*r_max)

    Cmin = np.min(C2D)
    Cmax = np.max(C2D)

    CvNorm = colors.Normalize(vmin=Cmin, vmax=Cmax)
    CscalarMap = cmx.ScalarMappable(norm=CvNorm, cmap=cmap)

    axs[0, 0].imshow(cgrid.display_grid.grid, cmap=cmap, norm=CvNorm,
        extent=[-r_max, r_max, -r_max, r_max])

    axs[0, 0].set_xlabel(r'$x$')
    axs[0, 0].set_ylabel(r'$y$')
    axs[0, 0].set_title('2D ' + r'$%s$' % C + ' ' +
        (r'$(%s^T/%s^L(\frac{r}{a} = %.3e) = %.3e)$'
        % (C, C, (box_size/Ncases)/parameters['a'],
        directional_correlations['CT']/directional_correlations['CL'])
        if 'CL' in directional_correlations
        and 'CT' in directional_correlations else ''))

    divider = make_axes_locatable(axs[0, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=CvNorm, orientation='vertical')
    cb.set_label(r'$%s$' % C, labelpad=20, rotation=270)

    # C1D shifted

    fplot(axs[1, 0])(C1D[1:, 0], C1D[1:, 1]/Cnn1D[-1, 1])

    axs[1, 0].set_xlabel(r'$r$')
    axs[1, 0].set_ylabel(r'$%s$' % C + r'$/C_{\rho\rho}(r=r_{max})$')
    axs[1, 0].set_title('radial ' + r'$%s$' % C + r'$/C_{\rho\rho}(r=r_{max})$'
        + ' ' + r'$(C_{\rho\rho}(r=r_{max}) = %.3e)$' % Cnn1D[-1, 1])

    axs[1, 0].set_xlim(r_min, r_max)
    axs[1, 0].set_ylim(C_min, C_max)

    # Cnn1D and C1D

    axs[0, 1].set_title('radial ' + r'$C_{\rho\rho}$' + ' and ' + r'$%s$' % C)
    axs[0, 1].set_xlabel(r'$r$')
    axs[0, 1].set_xlim(r_min, r_max)

    axs[0, 1].plot(Cnn1D[1:, 0], Cnn1D[1:, 1], color='#1f77b4')
    axs[0, 1].set_ylabel(r'$C_{\rho\rho}$', color='#1f77b4')
    axs[0, 1].tick_params('y', colors='#1f77b4')

    ax_right = axs[0, 1].twinx()
    ax_right.semilogy(C1D[1:, 0], C1D[1:, 1], color='#ff7f0e')
    ax_right.set_ylabel(r'$%s$' % C, color='#ff7f0e', rotation=270,
        labelpad=10)
    ax_right.tick_params('y', colors='#ff7f0e')
    ax_right.set_ylim(C_min*Cnn1D[-1, 1], C_max*Cnn1D[-1, 1])

    # C1D/Cnn

    fplot(axs[1, 1])(C1Dcor[1:, 0], C1Dcor[1:, 1])

    axs[1, 1].set_xlabel(r'$r$')
    axs[1, 1].set_ylabel(r'$%s$' % C + r'$/C_{\rho\rho}$')
    axs[1, 1].set_title('radial ' + r'$%s$' % C + r'$/C_{\rho\rho}$')

    axs[1, 1].set_xlim(r_min, r_max)
    axs[1, 1].set_ylim(C_min, C_max)

    # SAVING

    if get_env('SAVE', default=False, vartype=bool):	# SAVE mode
        image_name, = naming_standard.image().filename(**attributes)
        fig.savefig(joinpath(data_dir, image_name))

	# GRID CIRCLE

    if get_env('GRID_CIRCLE', default=False, vartype=bool):	# GRID_CIRCLE mode

        ccorgrid = CorGrid(C2D/Cnn2D, box_size, display_size=2*r_max)	# correlation corrected with density correlation

        gc = GridCircle(ccorgrid.display_grid.grid,
            extent=(-r_max, r_max, -r_max, r_max))
        gc.fig.set_size_inches(fig.get_size_inches())

        gc.fig.suptitle(suptitle)
        gc.fig.subplots_adjust(wspace=0.4)	# width space
        gc.fig.subplots_adjust(hspace=0.3)	# height space

        gc.ax_grid.set_xlabel(r'$x$')
        gc.ax_grid.set_ylabel(r'$y$')
        gc.ax_grid.set_title('2D ' + r'$%s$' % C + ' ' +
            (r'$(%s^T/%s^L(\frac{r}{a} = %.3e) = %.3e)$'
            % (C, C, (box_size/Ncases)/parameters['a'],
            directional_correlations['CT']/directional_correlations['CL'])
            if 'CL' in directional_correlations
            and 'CT' in directional_correlations else ''))

        gc.colormap.set_label(r'$%s$' % C, labelpad=20, rotation=270)

        gc.ax_plot.set_xlabel(r'$\theta$')
        gc.ax_plot.set_ylabel(r'$%s(\theta)$' % C)

    try:
        return fig, axs, gc
    except NameError: return fig, axs

# SCRIPT

if __name__ == '__main__':  # executing as script

    # VARIABLE DEFINITIONS

    data_dir = get_env('DATA_DIRECTORY', default=getcwd())	# data directory

    dt = get_env('TIME', default=-1, vartype=int)	# lag time for displacement

    init_frame = get_env('INITIAL_FRAME', default=-1, vartype=int)	# frame to consider as initial
    int_max = get_env('INTERVAL_MAXIMUM', default=1, vartype=int)	# maximum number of intervals of length dt considered in correlations calculations

    parameters_file = get_env('PARAMETERS_FILE',
		default=joinpath(data_dir, naming.parameters_file))	# simulation parameters file
    with open(parameters_file, 'rb') as param_file:
	       parameters = pickle.load(param_file)				# parameters hash table

    box_size = get_env('BOX_SIZE', default=parameters['box_size'],
		vartype=float)									# size of the square box to consider
    centre = (get_env('X_ZERO', default=0, vartype=float),
		get_env('Y_ZERO', default=0, vartype=float))	# centre of the box

    prep_frames = ceil(parameters['prep_steps']/parameters['period_dump'])	# number of preparation frames (FIRE energy minimisation)

    Ncases = get_env('N_CASES', default=ceil(np.sqrt(parameters['N'])),
		vartype=int)        # number of boxes in each direction with which to compute the displacement grid
    dL = box_size/Ncases    # boxes separation

    Nentries = parameters['N_steps']//parameters['period_dump']		# number of time snapshots in unwrapped trajectory file
    init_frame = int(Nentries/2) if init_frame < 0 else init_frame	# initial frame
    Nframes = Nentries - init_frame									# number of frames available for the calculation

    dt = Nframes + dt if dt <= 0 else dt	# length of the interval of time for which displacements are calculated

    attributes = {'density': parameters['density'],
		'vzero': parameters['vzero'], 'dr': parameters['dr'],
		'N': parameters['N'], 'init_frame': init_frame, 'dt': dt,
		'int_max': int_max, 'Ncases': Ncases, 'box_size': box_size,
        'x_zero': centre[0], 'y_zero': centre[1]}		# attributes displayed in filenames
    naming_Cnn = naming.Cnn()                           # Cnn naming object
    Cnn_filename, = naming_Cnn.filename(**attributes)   # Cnn filename
    naming_Cuu = naming.Cuu()                           # Cuu naming object
    Cuu_filename, = naming_Cuu.filename(**attributes)   # Cuu filename
    naming_Cww = naming.Cww()                           # Cww naming object
    Cww_filename, = naming_Cww.filename(**attributes)   # Cww filename
    naming_Cdd = naming.Cdd()                           # Cdd naming object
    Cdd_filename, = naming_Cdd.filename(**attributes)   # Cdd filename
    naming_Cee = naming.Cee()                           # Cee naming object
    Cee_filename, = naming_Cee.filename(**attributes)   # Cee filename

	# STANDARD OUTPUT

    if 'SLURM_JOB_ID' in envvar:	# script executed from Slurm job scheduler
        slurm_output(joinpath(data_dir, 'out'), naming_Cuu, attributes)

    # MODE SELECTION

    if get_env('COMPUTE', default=False, vartype=bool):	# COMPUTE mode

        startTime = datetime.now()

		# VARIABLE DEFINITIONS

        wrap_file_name = get_env('WRAPPED_FILE',
			default=joinpath(data_dir, naming.wrapped_trajectory_file))		# wrapped trajectory file (.gsd)
        unwrap_file_name = get_env('UNWRAPPED_FILE',
			default=joinpath(data_dir, naming.unwrapped_trajectory_file))	# unwrapped trajectory file (.dat)

        times = np.array(list(OrderedDict.fromkeys(map(
			lambda x: int(x),
			np.linspace(init_frame, Nentries - dt - 1, int_max)
			))))	# frames at which shear strain will be calculated

        # DISPLACEMENT CORRELATIONS

        with open(wrap_file_name, 'rb') as wrap_file,\
			open(unwrap_file_name, 'rb') as unwrap_file:	# opens wrapped and unwrapped trajectory files

            w_traj = Gsd(wrap_file, prep_frames=prep_frames)	# wrapped trajectory object
            u_traj = Dat(unwrap_file, parameters['N'])			# unwrapped trajectory object
            DDgrid, Ugrid, Wgrid, Egrid = tuple(np.transpose(list(map(
                lambda time: displacement_related_grids(
                	box_size, centre, Ncases, time, dt, w_traj, u_traj),
                times)), (1, 0, 2, 3, 4)))						# lists of displacement variables

        Dgrid = DDgrid[:, :, :, 0]  				# list of displacement norm grids
        Cdd2D = corField2D_scalar_average(Dgrid)	# displacement norm correlation grids

        Cnn_object = Cnn(Ugrid, box_size)	# density correlation object
        Cnn2D = Cnn_object.cnn2D			# 2D density correlation grid
        Cnn1D = Cnn_object.cnn1D			# 1D averaged density correlation grid

        (Cuu2D, CuuL, CuuT), (Cww2D, CwwL, CwwT), (Cee2D, CeeL, CeeT) = tuple(
            map(lambda Grid: corField2D_vector_average_Cnn(Grid, Cnn2D),
            [Ugrid, Wgrid, Egrid]))                                             # displacement, relative displacement and displacement direction correlation grids

        (Cuu1D, Cuu1Dcor), (Cww1D, Cww1Dcor), (Cdd1D, Cdd1Dcor),\
            (Cee1D, Cee1Dcor) = tuple(map(
            lambda C2D:
            tuple(map(lambda C: g2Dto1Dsquare(C, box_size),
            [C2D, np.divide(C2D, Cnn2D, out=np.zeros(C2D.shape),
            where=Cnn2D!=0)]
            )), [Cuu2D, Cww2D, Cdd2D, Cee2D]))  # 1D displacement variables correlations

        # SAVING

		# density correlations
        Cnn_object.save(attributes, dir=data_dir)
		# everything else
        with open(joinpath(data_dir, Cuu_filename), 'wb') as Cuu_dump_file,\
            open(joinpath(data_dir, Cww_filename), 'wb') as Cww_dump_file,\
            open(joinpath(data_dir, Cdd_filename), 'wb') as Cdd_dump_file,\
            open(joinpath(data_dir, Cee_filename), 'wb') as Cee_dump_file:
            pickle.dump([Cuu2D, Cuu1D, Cuu1Dcor, CuuL, CuuT], Cuu_dump_file)
            pickle.dump([Cww2D, Cww1D, Cww1Dcor, CwwL, CwwT], Cww_dump_file)
            pickle.dump([Cdd2D, Cdd1D, Cdd1Dcor], Cdd_dump_file)
            pickle.dump([Cee2D, Cee1D, Cee1Dcor, CeeL, CeeT], Cee_dump_file)

        # EXECUTION TIME

        print("Execution time: %s" % (datetime.now() - startTime))

    if get_env('PLOT', default=False, vartype=bool):	# PLOT mode

		# DATA

        with open(joinpath(data_dir, Cnn_filename), 'rb') as Cnn_dump_file,\
            open(joinpath(data_dir, Cuu_filename), 'rb') as Cuu_dump_file,\
            open(joinpath(data_dir, Cww_filename), 'rb') as Cww_dump_file,\
            open(joinpath(data_dir, Cdd_filename), 'rb') as Cdd_dump_file,\
            open(joinpath(data_dir, Cee_filename), 'rb') as Cee_dump_file:
            Cnn2D, Cnn1D = pickle.load(Cnn_dump_file)
            Cuu2D, Cuu1D, Cuu1Dcor, CuuL, CuuT = pickle.load(Cuu_dump_file)
            Cww2D, Cww1D, Cww1Dcor, CwwL, CwwT = pickle.load(Cww_dump_file)
            Cdd2D, Cdd1D, Cdd1Dcor = pickle.load(Cdd_dump_file)
            Cee2D, Cee1D, Cee1Dcor, CeeL, CeeT = pickle.load(Cee_dump_file)

    if get_env('PLOT', default=False, vartype=bool) or\
		get_env('SHOW', default=False, vartype=bool):	# PLOT or SHOW mode

		# PLOT

        plot_axis = get_env('AXIS', default='LOGLOG')   # plot scales
        fplot = lambda ax: ax.loglog if plot_axis == 'LOGLOG'\
            else ax.semilogy if plot_axis == 'LINLOG'\
            else ax.semilogx if plot_axis == 'LOGLIN'\
            else ax.plot

        r_min = get_env('R_MIN', default=_r_min, vartype=float) # minimum radius for correlations plots
        r_max = get_env('R_MAX', default=_r_max, vartype=float)	# maximum radius for correlations plots
        r_max = box_size/2 if r_max < 0 else r_max	# half size of the box showed for 2D correlation

        Cuu_min = get_env('CUU_MIN', default=_Cuu_min, vartype=float)   # minimum displacement correlation for correlation plots
        Cuu_max = get_env('CUU_MAX', default=_Cuu_max, vartype=float)   # maximum displacement correlation for correlation plots

        Cww_min = get_env('CWW_MIN', default=_Cww_min, vartype=float)   # minimum relative displacement correlation for correlation plots
        Cww_max = get_env('CWW_MAX', default=_Cww_max, vartype=float)   # maximum relative displacement correlation for correlation plots

        Cdd_min = get_env('CDD_MIN', default=_Cdd_min, vartype=float)   # minimum displacement norm correlation for correlation plots
        Cdd_max = get_env('CDD_MAX', default=_Cdd_max, vartype=float)   # maximum displacement norm correlation for correlation plots

        Cee_min = get_env('CEE_MIN', default=_Cee_min, vartype=float)   # minimum displacement direction correlation for correlation plots
        Cee_max = get_env('CEE_MAX', default=_Cee_max, vartype=float)   # maximum displacement direction correlation for correlation plots

        plot_Cuu = plot_correlation('C_{uu}', Cuu2D, Cuu1D, Cuu1Dcor, Cuu_min,
			Cuu_max, naming_Cuu,
			CL=CuuL, CT=CuuT)
        plot_Cww = plot_correlation('C_{\delta u \delta u}', Cww2D, Cww1D,
			Cww1Dcor, Cww_min, Cww_max, naming_Cww,
			CL=CwwL, CT=CwwT)
        plot_Cdd = plot_correlation('C_{|u||u|}', Cdd2D, Cdd1D, Cdd1Dcor,
			Cdd_min, Cdd_max, naming_Cdd)
        plot_Cee = plot_correlation('C_{\hat{u}\hat{u}}', Cee2D, Cee1D,
			Cee1Dcor, Cee_min, Cee_max, naming_Cee,
			CL=CeeL, CT=CeeT)

        if get_env('SHOW', default=False, vartype=bool):	# SHOW mode
            plt.show()
