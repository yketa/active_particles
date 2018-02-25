#! /home/yketa/miniconda3/bin/python3.6

import numpy as np
import pickle
import os

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else '/home/yketa/hoomd' # data directory

param_file = data_dir + '/' + os.environ['PARAMETERS_FILE'] if 'PARAMETERS_FILE' in os.environ else data_dir + '/param.pickle' # parameters file
vel_file = data_dir + '/' + os.environ['VELOCITY_FILE'] if 'VELOCITY_FILE' in os.environ else data_dir + '/velocity.csv' # velocity data file
varN_file = data_dir + '/' + os.environ['OUTPUT_FILE'] if 'OUTPUT_FILE' in os.environ else data_dir + '/varN.pickle' # number variance output file

file = open(param_file, 'rb')
N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(file)
file.close()

time = int(eval(os.environ['TIME'])) if 'TIME' in os.environ else -1 # time to evaluate the spatial variance of the number of particles
Ncases = int(eval(os.environ['N_CASES'])) if 'N_CASES' in os.environ else int(np.sqrt(N)) + (1 - int(np.sqrt(N))%2) # number of cases in each direction to compute the number of neighbouring particles
eval_points = int(eval(os.environ['EVALUATION_POINTS'])) if 'EVALUATION_POINTS' in os.environ else 100 # number of points of evaluation

data = np.genfromtxt(fname=vel_file, delimiter=',', skip_header=time if time >= 0 else N_steps//period_dump - 1, max_rows=1) # positions and velocity at time time (time = -1 if time < 0 in input)
positions = np.array([data[particle*3:particle*3 + 2] for particle in range(N)]) # array of the positions at time time

number = lambda point, size: np.sum((abs((positions - np.array(point) + box_size/2)%box_size - box_size/2) <= size/2).all(axis=1)) # returns the number of particles in a square of size size centred on the point point
points = lambda nb: np.array([(x, y) for y in np.linspace(- box_size/2, box_size/2, nb, endpoint=False) for x in np.linspace(- box_size/2, box_size/2, nb, endpoint=False)]) # returns nb x nb points linearly placed in the box
varNumber = lambda size, nb: (lambda list: np.mean(list**2) - np.mean(list))(np.array(list(map(lambda point: number(point, size), points(nb))))) # returns the number variance for points in squares of size size centred around points linearly spaced in the box

out_file = open(varN_file, 'wb')
size, varN = (lambda array: (array, list(map(lambda size: varNumber(size, Ncases), array))))(np.linspace(2*a, box_size, eval_points))
pickle.dump([size, varN], out_file)
out_file.close()

