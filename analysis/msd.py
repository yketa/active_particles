import numpy as np
import pickle
import os

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else '/home/yketa/hoomd' # data directory

def mean_square_displacement(param_file='param.pickle', pos_file='position.csv'):
	# Returns mean square displacements at associated times.

	file = open(param_file, 'rb')
	N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(file)
	file.close()

	data_pos = np.genfromtxt(fname=pos_file, delimiter=',')
	positions = lambda time: np.array([data_pos[time, particle*2:particle*2 + 2] for particle in range(N)]) # function giving all the positions at a given time

	time = [0] # array of increment of time
	msdis = [0] # array of mean square displacements

	for dt in range(1, len(data_vel)): # loop over all different increment of time
		time += [time_step*period_dump*dt]
		msdis += [np.mean([np.sum(((positions(time + dt) - np.mean(positions(time + dt), axis=0)) - (positions(time) - np.mean(positions(time), axis=0)))**2, axis=1) for time in range(len(data_vel) - dt)])]

	return time, msdis

file = open(data_dir + '/msd.pickle', 'wb')
time, msdis = mean_square_displacement(data_dir + '/param.pickle', data_dir + '/position.csv') # calculation of the mean square displacement
pickle.dump([time, msdis], file)
file.close()
