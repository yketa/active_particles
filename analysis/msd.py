#! /home/yketa/miniconda3/bin/python3.6

import numpy as np

import pickle

import os
import sys

sys.path.append('/home/yketa')
from exponents import float_to_letters

def mean_square_displacement(data_dir='/home/yketa/hoomd', param_file='param.pickle', pos_file='position.csv', init_frame=-1, snap_max=10, snap_period=1):
	# Returns mean square displacements at associated times.

	file = open(param_file, 'rb')
	N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(file)
	file.close()

	Nentries = N_steps//period_dump # number of time snapshots in velocity and position files
	init_frame = int(Nentries/2) if init_frame < 0 else init_frame # initial frame
	Ntimes = (Nentries - init_frame)//snap_period # number of time snapshots considered in the calculation

	pos = lambda time: np.reshape(np.genfromtxt(fname=pos_file, delimiter=',', skip_header=init_frame + time*snap_period, max_rows=1)[:-1], (N, 2)) # array of positions at time time
	times = lambda N, n: np.linspace(0, N - n - 1, min(snap_max, N - n), dtype=int) # list of time snapshots at which to evaluate the mean square displacement at n*dumps for N time snapshots

	out_file = data_dir + '/msd_std_D' + float_to_letters(density) + '_V' + float_to_letters(vzero) + '_R' + float_to_letters(dr) + '_N' + float_to_letters(N) + '_I' + float_to_letters(init_frame) + '_M' + float_to_letters(snap_max) + '_P' + float_to_letters(snap_period) + '.csv' # MSD file
	msd_dump = open(out_file, 'w') # MSD dump

	for dt in range(1, Ntimes): # loop over all different increment of time
		# CALCULATION
		time = time_step*period_dump*dt*snap_period # time at which is evaluated the mean square displacement
		msdis, stdev = (lambda val: (np.mean(val), np.std(val)))([np.sum(((lambda positions: positions - np.mean(positions, axis=0))(pos(time + dt)) - (lambda positions: positions - np.mean(positions, axis=0))(pos(time)))**2, axis=1) for time in times(Ntimes, dt)])
		# OUTPUT
		msd_dump.write(str("%e,%e,%e\n" % (time, msdis, stdev)))

	msd_dump.close()

	return time, msdis

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else os.getcwd() # data directory
param_file = os.environ['PARAMETERS_FILE'] if 'PARAMETERS_FILE' in os.environ else data_dir + '/param.pickle' # parameters file
pos_file = os.environ['POSITION_FILE'] if 'POSITION_FILE' in os.environ else data_dir + '/position.csv' # positions file

init_frame = int(eval(os.environ['INITIAL_FRAME'])) if 'INITIAL_FRAME' in os.environ else -1 # initial time for the calculation of the mean square displacement
snap_max = int(eval(os.environ['SNAP_MAXIMUM'])) if 'SNAP_MAXIMUM' in os.environ else 10 # maximum number of time snapshots taken for the calculation of the mean square displacement at each time
snap_period = int(eval(os.environ['SNAP_PERIOD'])) if 'SNAP_PERIOD' in os.environ else 1 # mean square displacement will be calculated for each snap_period dumps period of time

mean_square_displacement(data_dir=data_dir, param_file=param_file, pos_file=pos_file, init_frame=init_frame, snap_max=snap_max, snap_period=snap_period) # calculation of the mean square displacement
