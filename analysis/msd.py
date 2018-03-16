#! /home/yketa/miniconda3/bin/python3.6

import numpy as np

import pickle

import os
import sys

sys.path.append('/home/yketa')
from exponents import float_to_letters
sys.path.append('/home/yketa/hoomd/colmig_DPD_P_A/data')
from readdat import *

from collections import OrderedDict

def mean_square_displacement(data_dir='/home/yketa/hoomd', parameters_file='param.pickle', unwrapped_file='trajectory.dat', init_frame=-1, snap_max=1, snap_period=1):
	# Returns mean square displacements at associated times.

	with open(parameters_file, 'rb') as param_file:
		N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(param_file)

	Nentries = N_steps//period_dump # number of time snapshots in velocity and position files
	init_frame = int(Nentries/2) if init_frame < 0 else init_frame # initial frame
	Nframes = Nentries - init_frame # number of frames available for the calculation
	Ntimes = Nframes//snap_period # number of time intervals considered in the calculation

	pos = lambda time: np.reshape(np.genfromtxt(fname=pos_file, delimiter=',', skip_header=init_frame + time, max_rows=1)[:-1], (N, 2)) # array of positions at time time

	times = lambda Nframes, dt: np.linspace(0, Nframes - dt - 1, min(snap_max, Nframes - dt), dtype=int) # list of time snapshots at which to evaluate the mean square displacement at n*dumps for N time snapshots
	intervals = lambda Nframes, Ntimes: np.array(list(OrderedDict.fromkeys(map(lambda x: int(x), np.exp(np.linspace(np.log(1), np.log(Nframes - 1), Ntimes)))))) # intervals of times logarithmically spaced for the calculation

	out_file = data_dir + '/msd_std_D' + float_to_letters(density) + '_V' + float_to_letters(vzero) + '_R' + float_to_letters(dr) + '_N' + float_to_letters(N) + '_I' + float_to_letters(init_frame) + '_M' + float_to_letters(snap_max) + '_P' + float_to_letters(snap_period) + '.csv' # MSD file

	time = list(map(lambda dt: time_step*period_dump*dt, intervals(Nframes, Ntimes)))
	with open(unwrapped_file, 'rb') as unwrap_file:
		pos = lambda time: getarray(unwrap_file, N, time) # positions of the particles at time time

		def displacements(time, dt): # square displacements of the particles between time and time + dt
			wo_drift = lambda values: values - np.mean(values, axis=0) # deviation from the mean of array values
			return np.sum((wo_drift(pos(time + dt)) - wo_drift(pos(time)))**2, axis=1)
		def msd_sterr(dt): # mean square displacement and standard error for time intervals of size dt
			sterr = lambda values: np.std(values)/np.sqrt(np.prod(values.shape)) # standard error of array values
			return (lambda disp: (np.mean(disp), sterr(disp)))(
				np.array(list(map(lambda time: displacements(time, dt), times(Nframes, dt))))
			)

		msd, sterr = np.transpose(list(map(lambda dt: msd_sterr(dt), intervals(Nframes, Ntimes))))

	with open(out_file, 'w') as msd_dump:
		msd_dump.write("time, MSD, sterr\n") # header
		for t, m, s in zip(time, msd, sterr):
				msd_dump.write("%e,%e,%e\n" % (t, m, s)) # writing output to csv file

	return time, msd, sterr

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else os.getcwd() # data directory
parameters_file = os.environ['PARAMETERS_FILE'] if 'PARAMETERS_FILE' in os.environ else data_dir + '/param.pickle' # parameters file
unwrapped_file = os.environ['UNWRAPPED_FILE'] if 'UNWRAPPED_FILE' in os.environ else data_dir + '/trajectory.dat' # trajectory binary file

init_frame = int(eval(os.environ['INITIAL_FRAME'])) if 'INITIAL_FRAME' in os.environ else -1 # initial time for the calculation of the mean square displacement
snap_max = int(eval(os.environ['SNAP_MAXIMUM'])) if 'SNAP_MAXIMUM' in os.environ else 1 # maximum number of time snapshots taken for the calculation of the mean square displacement at each time
snap_period = int(eval(os.environ['SNAP_PERIOD'])) if 'SNAP_PERIOD' in os.environ else 1 # mean square displacement will be calculated for each snap_period dumps period of time

mean_square_displacement(data_dir=data_dir, parameters_file=parameters_file, unwrapped_file=unwrapped_file, init_frame=init_frame, snap_max=snap_max, snap_period=snap_period) # calculation of the mean square displacement

