#! /home/yketa/miniconda3/bin/python3.6

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
import operator
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable
sys.path.append('/home/yketa')
from exponents import *

N = int(eval(os.environ['N'])) if 'N' in os.environ else int(1e5)
density = int(eval(os.environ['DENSITY'])) if 'DENSITY' in os.environ else 0.8
vzero = int(eval(os.environ['VZERO'])) if 'VZERO' in os.environ else 1e-2

dr_min = float(eval(os.environ['DR_MIN'])) if 'DR_MIN' in os.environ else 1e-5
dr_max = float(eval(os.environ['DR_MAX'])) if 'DR_MAX' in os.environ else 1e-2

phimax = float(eval(os.environ['PHIMAX'])) if 'PHIMAX' in os.environ else 1

Pphilocmin = float(eval(os.environ['PPHILOC_MIN'])) if 'PPHILOC_MIN' in os.environ else 1e-4
Pphilocmax = float(eval(os.environ['PPHILOC_MAX'])) if 'PPHILOC_MAX' in os.environ else 1e-1
contours = int(eval(os.environ['CONTOURS'])) if 'CONTOURS' in os.environ else 20

init_frame = int(eval(os.environ['INITIAL_FRAME'])) if 'INITIAL_FRAME' in os.environ else 5000
int_max = int(eval(os.environ['INTERVAL_MAXIMUM'])) if 'INTERVAL_MAXIMUM' in os.environ else 1
Ncases = int(eval(os.environ['N_CASES'])) if 'N_CASES' in os.environ else 500
max_box_size = float(eval(os.environ['MAX_BOX_SIZE'])) if 'MAX_BOX_SIZE' in os.environ else 10

Nbins = int(eval(os.environ['N_BINS'])) if 'N_BINS' in os.environ else 200
bins = np.linspace(0, phimax, Nbins, endpoint=False)

dirs = [dir for dir in os.listdir() if str('D%s_V%s' % tuple(map(float_to_letters, [density, vzero]))) in dir and str('N%s_Ll0000' % float_to_letters(N)) in dir]
filename = lambda dir: str(dir + '/varN_' + dir[:-7] + '_I%s_M%s_C%s_B%s.pickle' % tuple(map(float_to_letters, [init_frame, int_max, Ncases, max_box_size])))
dr, densities = {}, {}
for dir in dirs:
	with open(dir + '/param.pickle', 'rb') as p_file, open(filename(dir), 'rb') as d_file:
		dr[dir] = pickle.load(p_file)[10]
		densities[dir] = pickle.load(d_file)
dirs = sorted(dirs, key=lambda dir: 1/dr[dir])

histogram = []
for bin in bins:
	for dir in dirs:
		Pphiloc = (np.sum(densities[dir] >= bin) - np.sum(densities[dir] > (bin + phimax/Nbins)))/len(densities[dir])
		histogram += [[np.log10(1/dr[dir]), bin, np.log10(Pphiloc) if Pphiloc > Pphilocmin else np.log10(Pphilocmin)]]
histogram = np.array(histogram)

fig, ax = plt.subplots()
fig.set_size_inches(10, 10)

colormap = os.environ['COLORMAP'] if 'COLORMAP' in os.environ else 'inferno'
cm = plt.get_cmap(colormap)
Norm  = colors.Normalize(vmin=np.log10(Pphilocmin), vmax=np.log10(Pphilocmax))
scalarMap = cmx.ScalarMappable(norm=Norm, cmap=cm)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = mpl.colorbar.ColorbarBase(cax, cmap=cm, norm=Norm, orientation='vertical')
cb.set_label(r'$\log P(\phi_{loc})$', labelpad=20, rotation=270)

ax.tricontourf(histogram[:, 0], histogram[:, 1], histogram[:, 2], contours, cmap=cm, norm=Norm)
ax.plot(histogram[:, 0], [density]*len(histogram), linestyle='--', color='red', linewidth=4)

ax.set_title(r'$N=%.1e, \phi=%1.2f, \tilde{v} = %.2e$' % (N, density, vzero) + '\n' + r'$S_{init} = %.1e, S_{max} = %.1e, N_{cases} = %.1e, r_{max} = %.1e$' % (init_frame, int_max, Ncases, max_box_size))
ax.set_xlabel(r'$\log( \tau_r \equiv \tilde{\nu}_r^{-1})$')
ax.set_ylabel(r'$\phi_{loc}$')
plt.show()

