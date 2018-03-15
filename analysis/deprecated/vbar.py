#! /home/yketa/miniconda3/bin/python3.6

import os
import numpy as np

input = os.environ['INPUT_FILE'] if 'INPUT_FILE' in os.environ else '/home/yketa/hoomd/velocity.csv'
data = np.genfromtxt(fname=input, delimiter=',')

N = data.shape[1]//6

velocities = np.array([[data[time, N*3 + particle*3:N*3 + particle*3 + 2]/np.sqrt(np.sum(data[time, N*3 + particle*3:N*3 + particle*3 + 2]**2)) for particle in range(N)] for time in range(len(data))])

vbar = np.mean(np.sqrt(np.sum(np.sum(velocities, axis=1)**2, axis=1))/N)
print(vbar)
