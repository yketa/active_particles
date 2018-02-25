#! /home/yketa/miniconda3/bin/python3.6

# ALGORITHM BASED ON FILY, HENKES AND MARCHETTI, SOFT MATTER, 2014, 10, 2132

import os

import subprocess

import hoomd
import hoomd.md

import numpy as np

import pickle

hoomd.context.initialize(''); # initialise hoomd

# PARAMETERS

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else '/home/yketa/hoomd/data' # name of the simulation data directory
call = subprocess.call(['mkdir', '-p', data_dir]) # creation of the simulation data directory if it does not exist

name_par = data_dir + '/' + os.environ['NAME_PAR'] if 'NAME_PAR' in os.environ else data_dir + '/param.pickle' # name of the simulation parameters saving file

name_log = data_dir + '/' + os.environ['NAME_LOG'] if 'NAME_LOG' in os.environ else data_dir + '/log-output.log' # name of the log output file
name_trajectory = data_dir + '/' + os.environ['NAME_TRAJECTORY'] if 'NAME_TRAJECTORY' in os.environ else data_dir + '/trajectory' # name of the trajectory output gsd file
name_vel = data_dir + '/' + os.environ['NAME_VELOCITY'] if 'NAME_VELOCITY' in os.environ else data_dir + '/velocity.csv' # name of the velocity output file
name_xml = data_dir + '/' + os.environ['NAME_XML'] if 'NAME_XML' in os.environ else '' # name of the trajectory output xml file

N = int(eval(os.environ['N'])) if 'N' in os.environ else 2000 # number of particles

a = float(eval(os.environ['MEAN_RADIUS'])) if 'MEAN_RADIUS' in os.environ else 1 # mean radius of the particles
pdi = float(eval(os.environ['POLYDISPERSITY'])) if 'POLYDISPERSITY' in os.environ else 0.2 # polydispersity index (between 0 and 1) for the uniformly distributed radii
N_sizes = int(eval(os.environ['N_SIZES'])) if 'N_SIZES' in os.environ else min(10, N) # number of different sizes of particles

density = float(eval(os.environ['DENSITY'])) if 'DENSITY' in os.environ else 0.7 # (area) density of the particles
box_size = np.sqrt((np.pi/density)*((N - 1)*(a**2)*((1 - pdi)**2) + 2*pdi*(1 - pdi)*N*(a**2) + (2/(3*(N - 1)))*(pdi**2)*N*(2*N - 1)*(a**2))) # size of the square box

kT = float(eval(os.environ['TEMPERATURE'])) if 'TEMPERATURE' in os.environ else 0 # temperature for actual dynamics

mu = float(eval(os.environ['MOBILITY'])) if 'MOBILITY' in os.environ else 1 # mobility of the particles
k = float(eval(os.environ['SPRING_CONSTANT'])) if 'SPRING_CONSTANT' in os.environ else 1 # spring constant of the harmonic contact repulsion of the particles

vzero = float(eval(os.environ['SELF_PROPULSION_SPEED']))*(a*mu*k) if 'SELF_PROPULSION_SPEED' in os.environ else 1e-2*(a*mu*k) # self-propulsion speed (to be entred in dimensionless units)
dr = float(eval(os.environ['ROTATION_DIFFUSION']))*(mu*k) if 'ROTATION_DIFFUSION' in os.environ else 5e-4*(mu*k) # rotational diffusion constant (to be entred in dimensionless units)

damp_bro = float(eval(os.environ['DAMPING_BROWNIAN'])) if 'DAMPING_BROWNIAN' in os.environ else 1 # damping for the Brownian dynamics

shear_rate = float(eval(os.environ['SHEAR_RATE'])) if 'SHEAR_RATE' in os.environ else 0 # shear rate along the x axis

time_step = float(eval(os.environ['TIME_STEP'])) if 'TIME_STEP' in os.environ else 1e-2 # integration time step
N_steps = int(eval(os.environ['N_STEPS'])) if 'N_STEPS' in os.environ else int(1e4) # number of integration steps

period_dump = int(eval(os.environ['PERIOD_DUMP'])) if 'PERIOD_DUMP' in os.environ else 100 # period of dumping to gsd file

init_gsd = os.environ['INITALISATION_GSD'] if 'INITALISATION_GSD' in os.environ else '' # initialisation gsd file
init_frame = int(eval(os.environ['INITIALISATION_FRAME'])) if 'INITIALISATION_FRAME' in os.environ else 0 # initialisation frame in the gsd file

if 'INITIALISATION_GSD' in os.environ and init_frame == -1: # initialise from the last frame
        init_gsd_file = gsd.pygsd.GSDFile(open(init_gsd, 'rb'));
        init_gsd_trajectory = gsd.hoomd.HOOMDTrajectory(init_gsd_file);
        init_frame = len(init_gsd_trajectory) - 1 # last frame

# BOX PARAMETRISATION

if 'INITIALISATION_GSD' in os.environ:
	system = hoomd.init.read_gsd(filename=init_gsd, frame=init_frame) # initialisation from gsd file
	snapshot = system.take_snapshot(all=True)
	N_cases = len(snapshot.particles.types) # number of different sizes
else:
	snapshot = hoomd.data.make_snapshot(N=N, particle_types=[str(particle) for particle in range(N_sizes)], box=hoomd.data.boxdim(L=box_size, dimensions=2)) # creating box with N particles
	snapshot.particles.position[:] = np.array([np.ndarray.tolist(np.array([np.random.rand(1)[0], np.random.rand(1)[0]])*box_size - box_size/2) + [0] for particle in range(N)]) # random position of the particles
	diameters = 2*np.linspace(a*(1 - pdi), a*(1 + pdi), N_sizes)
	snapshot.particles.diameter[:] = np.array(list(map(lambda particle: diameters[particle%N_sizes], range(N)))) # diameters of the particles (uniformly distributed polidispersity with polydispersity index pdi and average radius a)
	snapshot.particles.typeid[:] = np.array(list(map(lambda particle: particle%N_sizes, range(N)))) # index of the types of the particles
	system = hoomd.init.read_snapshot(snapshot) # initialisation from the created snapshot

# NEIGHBOR LIST

all = hoomd.group.all() # all particles
nl = hoomd.md.nlist.cell(); # neighbour list

# OUTPUT

hoomd.analyze.log(filename=name_log, quantities=['potential_energy', 'temperature', 'pressure_xx', 'pressure_xy', 'pressure_xz', 'pressure_yy', 'pressure_yz', 'pressure_zz','pressure','xy','lx','ly'], period=period_dump, overwrite=False if 'INITIALISATION_GSD' in os.environ else True)
hoomd.dump.gsd(filename=name_trajectory + '.gsd', period=period_dump, group=all, overwrite=False if 'INITIALISATION_GSD' in os.environ else True, dynamic=['momentum', 'attribute'])
hoomd.dump.dcd(filename=name_trajectory + '.dcd', period=period_dump, group=all, overwrite=False if 'INITIALISATION_GSD' in os.environ else True, dynamic=['momentum', 'attribute'])
if 'NAME_XML' in os.environ: hoomd.deprecated.dump.xml(filename=name_xml, period=period_dump, group=all)

# REPULSIVE FORCE

dpdc = hoomd.md.pair.dpd_conservative(nlist=nl, r_cut=0) # DPD conservative linear repulsive force
for type1 in range(N_sizes):
	for type2 in range(type1, N_sizes):
		dpdc.pair_coeff.set(str(type1), str(type2), A=mu*k*(snapshot.particles.diameter[type1] + snapshot.particles.diameter[type2])/2, r_cut=(snapshot.particles.diameter[type1] + snapshot.particles.diameter[type2])/2) # repulsive force parameters between particle1 and particle2

# MINIMISATION OF THE ENERGY [PREPARATION]

fire = hoomd.md.integrate.mode_minimize_fire(dt=time_step) # FIRE energy minimisation
nve = hoomd.md.integrate.nve(group=all)

prep_steps = 0
while not(fire.has_converged()):
	hoomd.run(100) # run FIRE until it has converged
	prep_steps += 100

# PARAMETERS FILE

par_file = open(name_par, 'w+b') # parameters saving file
pickle.dump([N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps], par_file)
par_file.close()

# INTEGRATION MODE

hoomd.md.integrate.mode_standard(dt=time_step); # enables integration with time step dt

# ACTUAL DYNAMICS

nve.disable() # disabling NVE integration
hoomd.md.integrate.brownian(group=all, kT=kT, dscale=damp_bro, seed=123); # integration with [overdamped] Brownian dynamics a temperature kT with gamma = damping*d_i where d_i is the diameter of particle i
hoomd.update.box_resize(xy=hoomd.variant.linear_interp([(0, 0), (N_steps, shear_rate*time_step*N_steps)])) # Lees-Edwards boundary conditions with constant shear rate

# SELF-PROPELLING FORCE INITIALISATION

activity = [] # active force applied on each particle
for particle in range(N):
	force = np.array([np.random.rand(1)[0] - 0.5 for axis in range(2)] + [0]) # non-normalised direction of the force
	force *= vzero/np.sqrt(np.sum(force**2)) # active force 
	activity += [tuple(force)]

propulsion = hoomd.md.force.active(group=all, seed=123, f_lst=activity, rotation_diff=dr, orientation_link=False)

# RUN

vel_dump = open(name_vel, 'w')
for runs in range(N_steps//period_dump):
	snap1 = system.take_snapshot(all=True)
	hoomd.run(1)
	snap2 = system.take_snapshot(all=True)
	# calculation and saving of the velocities
	velocities = (snap2.particles.position - snap1.particles.position)/time_step # manual calculation of the velocities
	dump_str = ""
	for value in list(snap1.particles.position.flatten()) + list(velocities.flatten()):
		dump_str += str("%e," % value)
	dump_str += "\n"
	vel_dump.write(dump_str)

	hoomd.run(period_dump - 1)
vel_dump.close()
