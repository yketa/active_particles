#! /home/yketa/miniconda3/bin/python3.6

import os

import numpy as np

import pickle

vel_file = os.environ['VELOCITY_FILE'] if 'VELOCITY_FILE' in os.environ else 'velocity.csv' # velocity file
param_file = os.environ['PARAMETERS_FILE'] if 'PARAMETERS_FILE' in os.environ else 'param.pickle' # parameters file

name_pos = os.environ['POSITION_FILE'] if 'POSITION_FILE' in os.environ else 'position.csv' # position file

file = open(param_file, 'rb')
N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(file)
file.close()

pos = lambda time: np.genfromtxt(fname=vel_file, delimiter=',', skip_header=time, max_rows=1)[:3*N] # function giving the positions
inc = lambda pos0, pos1: np.sign(pos0)*(pos1*pos0 < 0)*(abs(pos0) > box_size/4)*box_size # increment in coordinates that has to be applied to position array between times 0 and 1

pos_dump = open(name_pos, 'w')
def dump_positions(positions, pos_dump):
	# Writes positions to output file.
	for value in positions[:-1]:
		pos_dump.write(str("%e," % value))
	pos_dump.write("\n")

increments = np.zeros(3*N) # increments to be applied
pos0 = pos(0) # positions at time 0
dump_positions(pos0, pos_dump)
for time in range(1, N_steps//period_dump):
	pos1 = pos(time) # positions at time 1
	increments += inc(pos0, pos1)
	dump_positions(pos1 + increments, pos_dump)
	pos0 = np.copy(pos1)
pos_dump.close()
