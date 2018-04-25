#! /home/yketa/miniconda3/bin/python3.6

import os

import matplotlib as mpl
if not('SHOW' in os.environ and eval(os.environ['SHOW'])):
	mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

import gsd
import gsd.hoomd
import gsd.pygsd

import pickle

import sys

sys.path.append('/home/yketa')
from exponents import *

import math
from collections import OrderedDict

from datetime import datetime
startTime = datetime.now()

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else os.getcwd() # data directory

parameters_file = os.environ['PARAMETERS_FILE'] if 'PARAMETERS_FILE' in os.environ else data_dir + '/param.pickle' # parameters pickle file
wrap_file_name = os.environ['WRAPPED_FILE'] if 'WRAPPED_FILE' in os.environ else data_dir + '/trajectory.gsd' # trajectory gsd file

init_frame = int(eval(os.environ['INITIAL_FRAME'])) if 'INITIAL_FRAME' in os.environ else -1 # initial time for the calculation of the displacement correlation
int_max = int(eval(os.environ['INTERVAL_MAXIMUM'])) if 'INTERVAL_MAXIMUM' in os.environ else 1 # maximum number of intervals taken for the calculation of the displacement correlation

Ncases = int(eval(os.environ['N_CASES'])) if 'N_CASES' in os.environ else -1 # number of cases in each direction to compute the number of neighbouring particles

max_box_size = float(eval(os.environ['MAX_BOX_SIZE'])) if 'MAX_BOX_SIZE' in os.environ else 10 # size of the square box in which particles are counted

Nbins = int(eval(os.environ['N_BINS'])) if 'N_BINS' in os.environ else 10 # number of bins for the histogram
phimax = float(eval(os.environ['PHIMAX'])) if 'PHIMAX' in os.environ else 1

with open(parameters_file, 'rb') as param_file:
	N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(param_file)
Nentries = N_steps//period_dump # number of frames
prep_frames = prep_steps//period_dump # number of preparation frames

init_frame = int(Nentries/2) if init_frame < 0 else init_frame
frames = np.array(list(OrderedDict.fromkeys(map(lambda x: int(x), np.linspace(init_frame, Nentries - 1, int_max))))) # frames at which to calculate the histogram

Ncases = int(np.sqrt(N)) + (1 - int(np.sqrt(N))%2) if Ncases < 0 else Ncases
points = np.array([(x, y) for y in np.linspace(- box_size/2, box_size/2, Ncases, endpoint=False) for x in np.linspace(- box_size/2, box_size/2, Ncases, endpoint=False)])


with gsd.pygsd.GSDFile(open(wrap_file_name, 'rb')) as wrap_file:
	w_traj = gsd.hoomd.HOOMDTrajectory(wrap_file);

	surface = lambda frame, point: np.sum((np.pi/4)*(w_traj[int(prep_frames + frame)].particles.diameter**2)*(abs((w_traj[int(prep_frames + frame)].particles.position[:, :2] - np.array(point) + box_size/2)%box_size - box_size/2) <= max_box_size/2).all(axis=1)) # get the numbers of particle in the square of centre point and length max_box_size at time frame
	densities = np.array(list(map(lambda frame: list(map(lambda point: surface(frame, point), points)), frames))).flatten()/(max_box_size**2) # densities

filename = lambda ext: str(data_dir + '/varN_D%s_V%s_R%s_N%s_I%s_M%s_C%s_B%s.' % tuple(map(float_to_letters, [density, vzero, dr, N, init_frame, int_max, Ncases, max_box_size])) + ext) # filename
if 'SAVE' in os.environ and eval(os.environ['SAVE']):
	with open(filename('pickle'), 'wb') as dump_file:
		pickle.dump(densities, dump_file)

# PLOT

histogram = list(map(lambda bin: (np.sum(densities >= phimax*bin/Nbins) - np.sum(densities > phimax*(bin + 1)/Nbins))/len(densities), range(Nbins)))
plt.semilogy(np.linspace(0, phimax, len(histogram)), histogram, '.')
plt.title(r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, \tilde{\nu}_r=%.2e$' % (N, density, vzero, dr) + '\n' + r'$S_{init}=%.2e, S_{max}=%.2e, N_{cases}=%.2e, r_{max}=%.2e$' % (init_frame, int_max, Ncases, max_box_size))
plt.xlabel(r'$\phi_{loc}$')
plt.ylabel(r'$P(\phi_{loc})$')
if 'SAVEFIG' in os.environ and eval(os.environ['SAVEFIG']):
	plt.savefig(filename('eps'))
	print('Figure saved to %s.' % filename('eps'))
# EXECUTION TIME
print("Execution time: %s" % (datetime.now() - startTime))
if 'SHOW' in os.environ and eval(os.environ['SHOW']):
	plt.show()
