import numpy as np
import pickle
import os

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else '/home/yketa/hoomd' # data directory

def mean_square_displacement(param_file='param.pickle', vel_file='velocity.csv'):
	# Returns mean square displacements at associated times.

	file = open(param_file, 'rb')
	N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(file)
	file.close()

	data_vel = np.genfromtxt(fname=vel_file, delimiter=',')
	pos = lambda time, particle: data_vel[time, particle*3:particle*3 + 2] # function giving the positions

	inc = lambda time, particle: np.sign(pos(time - 1, particle))*(pos(time, particle)*pos(time - 1, particle) < 0)*(abs(pos(time - 1, particle)) > box_size/4)*box_size # increment in coordinates that has to be applied to particle particle between times time - 1 and time because of periodic boundary conditions

	increments = np.zeros((N, 2)) # increments to be applied
	positions = [[pos(0, particle) for particle in range(N)]] # array of the positions with the correct increments
	for time in range(1, len(data_vel)):
		increments += np.array([inc(time, particle) for particle in range(N)])
		positions += [np.array([pos(time, particle) for particle in range(N)]) + increments]
	positions = np.array(positions)

	time = [0] # array of increment of time
	msdis = [0] # array of mean square displacements

	for dt in range(1, len(data_vel)): # loop over all different increment of time
		time += [time_step*period_dump*dt]
		msdis += [np.mean([np.sum(((positions[time + dt] - np.mean(positions[time + dt], axis=0)) - (positions[time] - np.mean(positions[time], axis=0)))**2, axis=1) for time in range(len(data_vel) - dt)])]

	return time, msdis

file = open(data_dir + '/msd.pickle', 'wb')
time, msdis = mean_square_displacement(data_dir + '/param.pickle', data_dir + '/velocity.csv') # calculation of the mean square displacement
pickle.dump([time, msdis], file)
file.close()

