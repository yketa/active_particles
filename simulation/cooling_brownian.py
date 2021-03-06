#! /home/yketa/miniconda3/bin/python3.6

import os

import subprocess

import hoomd
import hoomd.md

import numpy as np

import pickle
import struct

hoomd.context.initialize(''); # initialise hoomd

# PARAMETERS

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else '/home/yketa/hoomd/colmig_DPD_P_A/data/cooling_brownian' # name of the simulation data directory
call = subprocess.call(['mkdir', '-p', data_dir]) # creation of the simulation data directory if it does not exist

name_par = data_dir + '/' + os.environ['NAME_PAR'] if 'NAME_PAR' in os.environ else data_dir + '/param.pickle' # name of the simulation parameters saving file
name_par_cool = data_dir + '/' + os.environ['NAME_PAR_COOL'] if 'NAME_PAR_COOL' in os.environ else data_dir + '/param_cool.pickle' # name of the simulation parameters specific to cooling saving file

name_log = data_dir + '/' + os.environ['NAME_LOG'] if 'NAME_LOG' in os.environ else data_dir + '/log-output.log' # name of the log output file
name_trajectory = data_dir + '/' + os.environ['NAME_TRAJECTORY'] if 'NAME_TRAJECTORY' in os.environ else data_dir + '/trajectory' # name of the trajectory output files (.gsd and .data)
name_xml = data_dir + '/' + os.environ['NAME_XML'] if 'NAME_XML' in os.environ else '' # name of the trajectory output xml file

N = int(eval(os.environ['N'])) if 'N' in os.environ else 2000 # number of particles

a = float(eval(os.environ['MEAN_RADIUS'])) if 'MEAN_RADIUS' in os.environ else 1 # mean radius of the particles
pdi = float(eval(os.environ['POLYDISPERSITY'])) if 'POLYDISPERSITY' in os.environ else 0.2 # polydispersity index (between 0 and 1) for the uniformly distributed radii
N_sizes = int(eval(os.environ['N_SIZES'])) if 'N_SIZES' in os.environ else min(10, N) # number of different sizes of particles

density = float(eval(os.environ['DENSITY'])) if 'DENSITY' in os.environ else 0.7 # (area) density of the particles
box_size = np.sqrt((np.pi/density)*((N - 1)*(a**2)*((1 - pdi)**2) + 2*pdi*(1 - pdi)*N*(a**2) + (2/(3*(N - 1)))*(pdi**2)*N*(2*N - 1)*(a**2))) # size of the square box

kT0 = float(eval(os.environ['TEMPERATURE_INITIAL'])) if 'TEMPERATURE_INITIAL' in os.environ else 1 # initial temperature
kT1 = float(eval(os.environ['TEMPERATURE_FINAL'])) if 'TEMPERATURE_FINAL' in os.environ else 0 # final temperature
prep_time = float(eval(os.environ['PREPARATION_TIME'])) if 'PREPARATION_TIME' in os.environ else 10 # preparation time of the sample at initial temperature
cool_rate = float(eval(os.environ['COOLING_RATE'])) if 'COOLING_RATE' in os.environ else 10 # cooling rate towards final temperature

mu = float(eval(os.environ['MOBILITY'])) if 'MOBILITY' in os.environ else 1 # mobility of the particles
k = float(eval(os.environ['SPRING_CONSTANT'])) if 'SPRING_CONSTANT' in os.environ else 1 # spring constant of the harmonic contact repulsion of the particles

vzero = float(eval(os.environ['SELF_PROPULSION_SPEED']))*(a*mu*k) if 'SELF_PROPULSION_SPEED' in os.environ else 0*(a*mu*k) # self-propulsion speed (to be entred in dimensionless units)
dr = float(eval(os.environ['ROTATION_DIFFUSION']))*(mu*k) if 'ROTATION_DIFFUSION' in os.environ else 0*(mu*k) # rotational diffusion constant (to be entred in dimensionless units)

damp_bro = float(eval(os.environ['DAMPING_BROWNIAN'])) if 'DAMPING_BROWNIAN' in os.environ else 1 # damping for the Brownian dynamics

time_step = float(eval(os.environ['TIME_STEP'])) if 'TIME_STEP' in os.environ else 1e-2 # integration time step

period_dump = int(eval(os.environ['PERIOD_DUMP'])) if 'PERIOD_DUMP' in os.environ else 100 # period of dumping to gsd file

init_gsd = os.environ['INITALISATION_GSD'] if 'INITALISATION_GSD' in os.environ else '' # initialisation gsd file
init_frame = int(eval(os.environ['INITIALISATION_FRAME'])) if 'INITIALISATION_FRAME' in os.environ else 0 # initialisation frame in the gsd file

if 'INITIALISATION_GSD' in os.environ and init_frame == -1: # initialise from the last frame
	with gsd.pygsd.GSDFile(open(init_gsd, 'rb')) as init_gsd_file:
		init_gsd_trajectory = gsd.hoomd.HOOMDTrajectory(init_gsd_file)
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

hoomd.analyze.log(filename=name_log, quantities=['potential_energy', 'temperature', 'pressure_xx', 'pressure_xy', 'pressure_xz', 'pressure_yy', 'pressure_yz', 'pressure_zz','pressure','xy','lx','ly'], period=period_dump, overwrite=False if 'INITIALISATION_GSD' in os.environ else True) # log output file
hoomd.dump.gsd(filename=name_trajectory + '.gsd', period=period_dump, group=all, overwrite=False if 'INITIALISATION_GSD' in os.environ else True, dynamic=['momentum', 'attribute']) # trajectory gsd file
if 'NAME_XML' in os.environ: hoomd.deprecated.dump.xml(filename=name_xml, period=period_dump, group=all)

def dump_trajectory(dump_file, N, positions, velocities):
	# writes data to binary trajectory file
	for data in [positions, velocities]:
		for particle in range(N):
			for coord in range(2):
					dump_file.write(struct.pack('d', data[particle, coord]))

snaps = [[None, hoomd.data.make_snapshot(N=N, box=hoomd.data.boxdim(L=box_size))], None] # list of system snapshots ([[time, time + period_dump], time + period_dump + 1])

increments = np.zeros((N, 3)) # increments of space to cancel wrapping due to periodic boundary conditions
inc = lambda increments, snaps, L: increments + np.sign(snaps[0][0].particles.position[:])*(snaps[0][0].particles.position[:]*snaps[0][1].particles.position[:] < 0)*(abs(snaps[0][0].particles.position[:]) > L/4)*L # function to update increments array

positions = lambda snaps, increments: snaps[0][1].particles.position[:] + increments # returns unwrapped positions array
velocities = lambda snaps, dt: (snaps[1].particles.position[:] - snaps[0][1].particles.position[:])/dt # returns velocities array

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

# PARAMETERS FILES

with open(name_par, 'wb') as par_file: # parameters saving file
	pickle.dump([N, a, pdi, N_sizes, density, box_size, kT0, mu, k, vzero, dr, damp_bro, 0, time_step, int((abs(kT1 - kT0)/(cool_rate * time_step))), period_dump, prep_steps], par_file)

with open(name_par_cool, 'wb') as par_cool_file: # parameters specific to cooling file
	pickle.dump([kT1, prep_time, cool_rate], par_cool_file)

# INTEGRATION MODE

hoomd.md.integrate.mode_standard(dt=time_step); # enables integration with time step dt

# SELF-PROPELLING FORCE INITIALISATION

activity = [] # active force applied on each particle
for particle in range(N):
	force = np.array([np.random.rand(1)[0] - 0.5 for axis in range(2)] + [0]) # non-normalised direction of the force
	force *= vzero/np.sqrt(np.sum(force**2)) # active force
	activity += [tuple(force)]

propulsion = hoomd.md.force.active(group=all, seed=123, f_lst=activity, rotation_diff=dr, orientation_link=False)

# PREPARATION

nve.disable() # disabling NVE integration
preparation = hoomd.md.integrate.brownian(group=all, kT=kT0, dscale=damp_bro, seed=123); # integration with [overdamped] Brownian dynamics a temperature kT0 with gamma = damping*d_i where d_i is the diameter of particle i
print('\nINITIALISATION AT T=%e\n' % kT0)
hoomd.run(prep_time/time_step) # running preparation steps
preparation.disable() # end preparation

# ACTUAL COOLING (or HEATING)

cooling = hoomd.md.integrate.brownian(group=all, kT=hoomd.variant.linear_interp([(0, kT0), (abs(kT0 - kT1)/(cool_rate * time_step), kT1)]), dscale=damp_bro, seed=123); # integration with [overdamped] Brownian dynamics, coolong from temperature kT0 to kT1, with gamma = damping*d_i where d_i is the diameter of particle i
print('\nCOOLING TO T=%e' % kT1 + ' AT COOLING RATE K=%e\n' % cool_rate)

# RUN

def run(dump_file, snaps, increments, N, L, dt, period_dump):
	# runs hoomd for period_dump iterations and saves corresponding data

	snaps[0][0] = snaps[0][1] # reference time in snaps increases by period_dump
	# for positions
	snaps[0][1] = system.take_snapshot(all=True)
	# for velocities
	hoomd.run(1)
	snaps[1] = system.take_snapshot(all=True)
	# updating increments
	increments = inc(increments, snaps, L)
	# dump
	dump_trajectory(dump_file, N, positions(snaps, increments), velocities(snaps, dt))
	# additional run
	hoomd.run(period_dump - 1)

	return snaps, increments

with open(name_trajectory + '.dat', 'wb') as output_trajectory: # trajectory data file
	for runs in range(int((abs(kT1 - kT0)/(cool_rate * time_step))//period_dump)):
		snaps, increments = run(output_trajectory, snaps, increments, N, box_size, time_step, period_dump)
