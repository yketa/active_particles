"""
Prints parameter from parameters file.

Simulation directory name is considered to be
active_particles.naming.sim_directory and parameters file name is considered to
be active_particles.naming.parameters_file.

Script is called according to the following scheme:
$ python param.py [DIRECTORY NAME] [PARAMETER NAME]
"""

import sys
import pickle
from active_particles.naming import sim_directory, parameters_file
from os.path import join as joinpath


try:
    with open(joinpath(sim_directory, sys.argv[1], parameters_file),
        'rb') as param_file:
        param = pickle.load(param_file)
except IndexError:          # sys.argv[1] does not exist
    print('No directory name found.')
    sys.exit()
except FileNotFoundError:   # parameters file not found
    print('Parameters file %s not found in directory %s.' %
        (sys.argv[1], parameters_file))
    sys.exit()

try:
    print(param[sys.argv[2]])
except IndexError:  # sys.argv[2] does not exist
    print('No parameter name found.')
except KeyError:    # param[sys.argv[2]] does not exist
    print('Parameter %s does not exist.' % sys.argv[2])
