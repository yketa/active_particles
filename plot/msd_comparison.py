#! /home/yketa/miniconda3/bin/python3.6

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider, RadioButtons
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches

import numpy as np
import os
import pickle
import sys
sys.path.append('/home/yketa')
from exponents import *
sys.path.append('/home/yketa/mpl_fitting_line')
from mpl_fitting_line import *

font_size = int(eval(os.environ['FONT_SIZE'])) if 'FONT_SIZE' in os.environ else 10
marker_size = int(eval(os.environ['MARKER_SIZE'])) if 'MARKER_SIZE' in os.environ else 20
mpl.rcParams.update({'font.size': font_size, 'lines.markersize': marker_size})

ncol_legend = int(eval(os.environ['NCOL_LEGEND'])) if 'NCOL_LEGEND' in os.environ else 1 # number of columns for the legend

wspace = float(eval(os.environ['WSPACE'])) if 'WSPACE' in os.environ else 0.2
hspace = float(eval(os.environ['HSPACE'])) if 'HSPACE' in os.environ else 0.05

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else os.getcwd() # data directory
os.chdir(data_dir) # change working directory to data directory

snap_max = int(eval(os.environ['SNAP_MAXIMUM'])) if 'SNAP_MAXIMUM' in os.environ else 1 # maximum number of time snapshots taken for the calculation of the mean square displacement at each time
snap_period = int(eval(os.environ['SNAP_PERIOD'])) if 'SNAP_PERIOD' in os.environ else 1 # mean square displacement will be calculated for each snap_period dumps period of time

with open(data_dir + '/param.pickle', 'rb') as param_file:
	N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(param_file)
av_p_sep = box_size/np.sqrt(N) # average particle separation

snip0 = 'msd_sterr_'
snip1 = '_M%s_P%s' % tuple(map(float_to_letters, [snap_max, snap_period]))
files = [file for file in os.listdir() if (snip0 in file) and (snip1 in file)] # msd csv
init = lambda file: eval(letters_to_float(file.split("_I")[1][:5])) # lag time from file name
init_list = sorted(list(map(init, files))) # list of lag times

slope0 = float(eval(os.environ['SLOPE'])) if 'SLOPE' in os.environ else 1 # default slope for slider
slope_min = float(eval(os.environ['SLOPE_MIN'])) if 'SLOPE_MIN' in os.environ else 0 # minimum slope for slider
slope_max = float(eval(os.environ['SLOPE_MAX'])) if 'SLOPE_MAX' in os.environ else 2 # maximum slope for slider

colormap = os.environ['COLORMAP'] if 'COLORMAP' in os.environ else 'jet' # colormap for curves
cm = plt.get_cmap(colormap)
cNorm  = colors.Normalize(vmin=0, vmax=len(init_list) - 1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
colors = {init_list[i]: scalarMap.to_rgba(i) for i in range(len(init_list))} # curves colors

# msd

msd = {} # Cuu1D with renormalisation by Cnn
for file in files:
	msd[file] = np.genfromtxt(fname=file, delimiter=',', skip_header=True)

# PLOT

fig, ax = plt.subplots()
fig.set_size_inches(30, 30)
fig.subplots_adjust(wspace=wspace)
fig.subplots_adjust(hspace=hspace)

ax.set_xlabel(r'$\Delta t$')
ax.set_ylabel(r'$<|\Delta r(\Delta t)|^2>$')

for file in sorted(files, key=init):
	ax.loglog(msd[file][:, 0]/av_p_sep, msd[file][:, 1], color=colors[init(file)], label=r'$S_{init} = %.0e$' % (init(file)))
ax.legend(loc=7, bbox_to_anchor=(0.5, -0.1))

fitting_line = FittingLine(ax, slope0, slope_min, slope_max) # interactive fitting line

title = r'$N=%.2e, \phi=%1.2f, $' % (N, density)
title += r'$\tilde{v}=%.2e, \tilde{\nu}_r=%.2e$' % (vzero, dr) if not('TEMPERATURE' in os.environ and eval(os.environ['TEMPERATURE'])) else r'$kT=%.2e, k=%.2e$' % (kT, k)
title += '\n'
title += r'$S_{max}=%.2e, S_{period}=%.2e$' % (snap_max, snap_period)
fig.suptitle(title)

plt.show()
