import numpy as np
import pickle
import os

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else '/home/yketa/hoomd' # data directory

def mean_square_displacement(param_file='param.pickle', pos_file='position.csv'):
	# Returns mean square displacements at associated times.

	file = open(param_file, 'rb')
	N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(file)
	file.close()

	Nentries = N_steps//period_dump # number of time snapshots in velocity and position files
	Ntimes = Nentries - int(Nentries/2) # number of time snapshots considered in the calculation

	time = [0] # array of increment of time
	msdis = [0] # array of mean square displacements

	pos = lambda time: np.reshape(np.genfromtxt(fname=pos_file, delimiter=',', skip_header=int(Nentries/2) + time, max_rows=1)[:-1], (N, 2))

	for dt in range(1, Ntimes): # loop over all different increment of time
		time += [time_step*period_dump*dt]
		msdis += [0]
		for time_ in range(Ntimes - dt):
			r_t_dt = (lambda positions: positions - np.mean(positions, axis=0))(pos(time + dt))
			r_t = (lambda positions: positions - np.mean(positions, axis=0))(pos(time))
			dr2 = (r_t_dt - r_t)**2
			msdis[-1] += np.mean(np.sum(dr2, axis=1))
		msdis[-1] /= Ntimes - dt
		# msdis += [np.mean([np.sum(((lambda positions: positions - np.mean(positions, axis=0))(pos(time + dt)) - (lambda positions: positions - np.mean(positions, axis=0))(pos(time)))**2, axis=1) for time in range(Ntimes - dt)])]

	return time, msdis

file = open(data_dir + '/msd.pickle', 'wb')
time, msdis = mean_square_displacement(data_dir + '/param.pickle', data_dir + '/position.csv') # calculation of the mean square displacement
pickle.dump([time, msdis], file)
file.close()
