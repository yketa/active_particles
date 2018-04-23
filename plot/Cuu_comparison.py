#! /home/yketa/miniconda3/bin/python3.6

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import pickle
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mp

import os
import sys
sys.path.append('/home/yketa')
from exponents import *

from scipy import stats as st

os.chdir('/home/yketa/hoomd/colmig_DPD_P_A/data')

font_size = int(eval(os.environ['FONT_SIZE'])) if 'FONT_SIZE' in os.environ else 15
mp.rcParams.update({'font.size': font_size})

ratio_legend = int(eval(os.environ['RATIO_LEGEND'])) if 'RATIO_LEGEND' in os.environ else 10

plot_axis = os.environ['AXIS'] if 'AXIS' in os.environ else 'LOGLOG'
fplot = lambda ax: ax.loglog if plot_axis == 'LOGLOG' else ax.semilogy if plot_axis == 'LINLOG' else ax.semilogx if plot_axis == 'LOGLIN' else ax.plot

density = float(eval(os.environ['DENSITY'])) if 'DENSITY' in os.environ else 0.80
vzero = float(eval(os.environ['VZERO'])) if 'VZERO' in os.environ else 1e-2
number = float(eval(os.environ['NUMBER'])) if 'NUMBER' in os.environ else 1e5
init = float(eval(os.environ['INITIAL_FRAME'])) if 'INITIAL_FRAME' in os.environ else 5000
smax = float(eval(os.environ['SNAP_MAXIMUM'])) if 'SNAP_MAXIMUM' in os.environ else 100
Ncases = int(eval(os.environ['N_CASES'])) if 'N_CASES' in os.environ else 500

Cuu_min = float(eval(os.environ['CUU_MIN'])) if 'CUU_MIN' in os.environ else 5e-3
r_min = float(eval(os.environ['R_MIN'])) if 'R_MIN' in os.environ else 1
r_max = float(eval(os.environ['R_MAX'])) if 'R_MAX' in os.environ else 20

thresh = float(eval(os.environ['THRESHOLD'])) if 'THRESHOLD' in os.environ else np.exp(-1)

drdt = (np.array([[i*(10**j) for i in [1, 2]] for j in range(-2, 3)])/2).flatten() if not 'SINGLE' in os.environ else [float(eval(os.environ['SINGLE']))]
# dr = ['h1000', 'h5000', 'i1000', 'i5000', 'j1000', 'j5000'] + (['k1000'] if vzero == 1e-2 and density == 0.8 else [])
dr = ['h1000', 'h2000', 'h5000', 'i1000', 'i2000', 'i5000', 'j1000', 'j2000', 'j5000'] if density == 1 else list(map(float_to_letters, [i*(10**j) for j in range(-5, -2) for i in range(1, 10)] + [1e-2]))
dt = ['l1000', 'l2000', 'l4000', 'l5000', 'm1000', 'm2000', 'm4000', 'm5000', 'n1000', 'n2000', 'n4000', 'n5000', 'o1000', 'o2000', 'o4000']

variable = os.environ['CORRELATION'] if 'CORRELATION' in os.environ else 'Cuu'
C = {'Cuu':'C_{uu}', 'Cww':'C_{\delta u \delta u}', 'Cdd':'C_{|u||u|}', 'Cee':'C_{\hat{u}\hat{u}}'}[variable]

colormap = os.environ['COLORMAP'] if 'COLORMAP' in os.environ else 'jet'

fname = lambda var, d, t: str(('/home/yketa/hoomd/colmig_DPD_P_A/data/D%s_V%s_R%s_N%s_Ll0000/' + var + 'b_D%s_V%s_R%s_N%s_I%s_T%s_M%s_C%s.pickle') % tuple(map(float_to_letters, [density, vzero, float(eval(letters_to_float(d))), number, density, vzero, float(eval(letters_to_float(d))), number, init, float(eval(letters_to_float(t))), smax, Ncases])))
Cuu, g1D, rho, Cuu_red, Cuu_red_linreg, Cuu_red_type = {}, {}, {}, {}, {}, {}
dr_list = []
for d in dr:
	for t in dt:
		if float(eval(letters_to_float(d))) * 10 * float(eval(letters_to_float(t))) in drdt:
			if not(d in dr_list):
				dr_list += [d]
			with open(fname(variable, d, t), 'rb') as Cuu_file, open(fname('Cnn', d, t), 'rb') as g_file:
				Cuu[(d, t)] = pickle.load(Cuu_file)[1]
				g1D[(d, t)], rho[(d, t)] = pickle.load(g_file)[1], 1
				Cuu_red[(d, t)] = np.array([[Cuu[(d, t)][i, 0], Cuu[(d, t)][i, 1]/(rho[(d, t)]*g1D[(d, t)][i, 1])] for i in range(len(Cuu[(d, t)])) if Cuu[(d, t)][i, 0] >= r_min and Cuu[(d, t)][i, 0] <= r_max and Cuu[(d, t)][i, 1]/(rho[(d, t)]*g1D[(d, t)][i, 1]) >= Cuu_min])
				Cuu_red_linreg[(d, t)] = np.array([st.linregress(np.log(Cuu_red[(d, t)][:, 0]), np.log(Cuu_red[(d, t)][:, 1])), st.linregress(Cuu_red[(d, t)][:, 0], np.log(Cuu_red[(d, t)][:, 1])), st.linregress(np.log(Cuu_red[(d, t)][:, 0]), np.log(-np.log(Cuu_red[(d, t)][:, 1])))])
				Cuu_red_type[(d, t)] = np.argmax(Cuu_red_linreg[(d, t)][:, 2]**2)

types = {0: 'powerlaw', 1: 'exponential', 2: 'delayed exponential'}
fit = lambda x, d, t: {0: np.exp(Cuu_red_linreg[(d, t)][0, 1])*(x**Cuu_red_linreg[(d, t)][0, 0]), 1: np.exp(Cuu_red_linreg[(d, t)][1, 0]*x + Cuu_red_linreg[(d, t)][1, 1]), 2: np.exp(- (x**Cuu_red_linreg[(d, t)][2, 0])*np.exp(Cuu_red_linreg[(d, t)][2, 1]))}[Cuu_red_type[(d, t)]]
cor_length = lambda d, t: {0: (thresh*np.exp(- Cuu_red_linreg[(d, t)][0, 1]))**(1/Cuu_red_linreg[(d, t)][0, 0]), 1: np.log(thresh*np.exp(- Cuu_red_linreg[(d, t)][1, 1]))/Cuu_red_linreg[(d, t)][1, 0], 2: (- np.exp(- Cuu_red_linreg[(d, t)][2, 1])*np.log(thresh))**(1/Cuu_red_linreg[(d, t)][2, 0])}[Cuu_red_type[(d, t)]]

default_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
colors = {dr[i]:default_colors[i] for i in range(len(dr))}
"""jet = cm = plt.get_cmap(colormap)
cNorm  = colors.Normalize(vmin=0, vmax=len(dr_list) - 1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
colors = {dr_list[i]:scalarMap.to_rgba(i) for i in range(len(dr_list))}"""
linestyle = {0: '.-', 1: 's-', 2: '*-'}

plot_lines, plot_columns = (2, 5) if not 'SINGLE' in os.environ else (1, 1)
# fig, ax = plt.subplots(plot_lines, plot_columns)
fig = plt.figure()
gs = GridSpec(plot_lines, plot_columns + 1, width_ratios=[1]*plot_columns + [plot_columns/ratio_legend])
ax = np.array([[plt.subplot(gs[i, j]) for  j in range(plot_columns)] for i in range(plot_lines)])
leg = plt.subplot(gs[:, -1])
for plot in range(plot_lines*plot_columns):
	axis = ax[plot//plot_columns, plot%plot_columns] if not 'SINGLE' in os.environ else ax

	drdt_list = [(d, t) for d in dr for t in dt if float(eval(letters_to_float(d))) * 10 * float(eval(letters_to_float(t))) == drdt[plot]]

	for d, t in drdt_list:
		# fplot(axis)(Cuu[(d, t)][:, 0], Cuu[(d, t)][:, 1], color=colors[d], label=str(r'$\tilde{\nu}_r = %s$' % letters_to_float(d)))
		# fplot(axis)(Cuu_red[(d, t)][:, 0], Cuu_red[(d, t)][:, 1], linestyle[Cuu_red_type[(d, t)]], color=colors[d])
		fplot(axis)(Cuu_red[(d, t)][:, 0], Cuu_red[(d, t)][:, 1], '.-', color=colors[d])
		if 'FIT' in os.environ and eval(os.environ['FIT']):
			fplot(axis)(Cuu_red[(d, t)][:, 0], list(map(lambda x: fit(x, d, t), Cuu_red[(d, t)][:, 0])), '--', color=colors[d], label=r'$\xi = %s$' % cor_length(d, t))

	axis.set_xlim([r_min, r_max])
	axis.set_ylim([Cuu_min, 1])

	axis.set_title(str(r'$\tilde{\nu}_r \Delta t = %.3e$' % drdt[plot]))
	axis.set_xlabel(r'$r$')
	axis.set_ylabel(r'$%s$' % C)
	# if 'LEGEND' in os.environ and eval(os.environ['LEGEND']):
	# 	axis.legend()

for axis in ax.flat if not 'SINGLE' in os.environ else [ax]:
	axis.label_outer()

lines = list(map(lambda d: Line2D([0], [0], color=colors[d], lw=2, label=r'$\tilde{\nu}_r = %.0e$' % float(eval(letters_to_float(d)))), dr_list))
leg.legend(handles=lines, loc='center')
leg.axis('off')

fig.suptitle(r'$N=%.2e, \phi=%1.2f, \tilde{v}=%.2e, S_{init}=%.2e, S_{max}=%.2e, N_{cases}=%.2e$' % (number, density, vzero, init, smax, Ncases) + '\n' + r'$%s,min} = %.2e, r_{min} = %.2e, r_{max} = %.2e$' % (C[:-1], Cuu_min, r_min, r_max))

fig.set_size_inches(12*plot_columns,6*plot_lines)
fig.subplots_adjust(wspace=0.1)
fig.subplots_adjust(hspace=0.15)

plt.show()
