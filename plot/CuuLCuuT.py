#! /home/yketa/miniconda3/bin/python3.6

import numpy as np
import math
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/yketa')
from exponents import *
import pickle
import operator
import os
from scipy import stats as st
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import matplotlib.cm as cmx
from collections import OrderedDict

var = os.environ['CORRELATION'] if 'CORRELATION' in os.environ else 'Cww'
C = {'Cuu':'C_{uu}', 'Cww':'C_{\delta u \delta u}', 'Cdd':'C_{|u||u|}', 'Cee':'C_{\hat{u}\hat{u}}'}[var]
var += 'b' if not('ENDPOINT' in os.environ and not(eval(os.environ['ENDPOINT']))) else ''

dr_min = float(eval(os.environ['DR_MIN'])) if 'DR_MIN' in os.environ else 1e-10
dr_max = float(eval(os.environ['DR_MAX'])) if 'DR_MAX' in os.environ else 1e10
dr_c = float(eval(os.environ['DR_C'])) if 'DR_C' in os.environ else 3e-4

colormap = os.environ['COLORMAP'] if 'COLORMAP' in os.environ else 'jet'
default_markers = ['.', 's', '*', 'o']

density = float(eval(os.environ['DENSITY'])) if 'DENSITY' in os.environ else 0.8
vzero = float(eval(os.environ['VZERO'])) if 'VZERO' in os.environ else 1e-2
number = float(eval(os.environ['NUMBER'])) if 'NUMBER' in os.environ else 1e5

dt_list = np.array([1, 2, 4, 5, 10, 20, 40, 50, 100, 200, 400, 500, 1000, 2000, 4000])
# dr_list = np.array([2e-5, 7e-5, 2e-4, 7e-4, 2e-3, 7e-3])

dirs = [dir for dir in os.listdir() if 'D%s_V%s' % tuple(map(float_to_letters, [density, vzero])) in dir and 'N%s_Ll' % float_to_letters(number) in dir]

r_min = float(eval(os.environ['R_MIN'])) if 'R_MIN' in os.environ else 1
r_max = float(eval(os.environ['R_MAX'])) if 'R_MAX' in os.environ else 20

init_frame = float(eval(os.environ['INITIAL_FRAME'])) if 'INITIAL_FRAME' in os.environ else 5000
N_cases = int(eval(os.environ['N_CASES'])) if 'N_CASES' in os.environ else 500
int_max = int(eval(os.environ['INTERVAL_MAXIMUM']))if 'INTERVAL_MAXIMUM' in os.environ else 100

pname = lambda dir: str('%s/param.pickle' % dir)
fname = lambda dir, t: str('%s/%s_%s_I%s_T%s_M%s_C%s.pickle' % ((dir, var, dir[:-7]) + tuple(list(map(float_to_letters, [init_frame, t, int_max, N_cases])))))
iname = lambda dir: str('%s/int%s_%s_I%s_M%s_C%s_RMIN%s_RMAX%s.pickle' % ((dir, var, dir[:-7]) + tuple(list(map(float_to_letters, [init_frame, int_max, N_cases, r_min, r_max])))))

box_size, drdtmax, dr, time_step, period_dump, ra, CL, CT = {}, {}, {}, {}, {}, {}, {}, {}
for dir in dirs:
	with open(iname(dir), 'rb') as i_file:
		drdtmax[dir] = pickle.load(i_file)[1]
	with open(pname(dir), 'rb') as p_file:
		box_size[dir], dr[dir], time_step[dir], period_dump[dir] = operator.itemgetter(5, 10, 13, 15)(pickle.load(p_file))
	ra[dir] = box_size[dir]/N_cases
	for t in dt_list:
		with open(fname(dir, t), 'rb') as C_file:
			CL[(dir, t)], CT[(dir, t)] = pickle.load(C_file)[-2:]

dr_list = sorted(list(OrderedDict.fromkeys([dr[dir] for dir in dr if dr[dir] >= dr_min and dr[dir] <= dr_max])))
dr0_list = sorted(list(OrderedDict.fromkeys([dr[dir] for dir in dr if dr[dir] < dr_min or dr[dir] > dr_max])))
time_step_list = sorted(list(OrderedDict.fromkeys([time_step[dir] for dir in time_step])))

jet = cm = plt.get_cmap(colormap)
cNorm  = colors.Normalize(vmin=0, vmax=len(dr_list) - 1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
jet0 = cm0 = plt.get_cmap('Greys')
cNorm0 = colors.Normalize(vmin=0, vmax=len(dr0_list))
scalarMap0 = cmx.ScalarMappable(norm=cNorm0, cmap=jet0)
colors = {**{dr_list[i]:scalarMap.to_rgba(i) for i in range(len(dr_list))}, **{dr0_list[i]:scalarMap0.to_rgba(i + 1) for i in range(len(dr0_list))}}

markers = {time_step_list[i]:default_markers[i] for i in range(len(time_step_list))}

def plot(C, CL, CT):

	# RATIO

	fig = plt.figure()
	fig.suptitle(r'$\phi=%1.2f, \tilde{v}=%.2e, N=%.2e$' % (density, vzero, number)+ '\n' + r'$S_{init} = %.2e, S_{max}=%.2e, N_{cases}=%.2e$' % (init_frame, N_cases, int_max))

	if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']):
		gs = GridSpec(2, 2, width_ratios=[10, 1])
		ax = plt.subplot(gs[0, 0])

		ax1 = plt.subplot(gs[1, 0])
		ax1.set_xlabel(r'$\tau_r \equiv \tilde{\nu}_r^{-1}$')
		ax1.set_ylabel(r'$%s^T/%s^L(\frac{r}{a}=%.3e, \Delta t^*)$' % (C, C, ra[dirs[0]]))
		ax1.plot(np.linspace(1/np.max(np.append(dr_list, dr0_list)), 1/np.min(np.append(dr_list, dr0_list)), 100), [1]*100, color='black', linestyle='--')
		ax1.set_xscale('log')

		lines1 = list(map(lambda dt: Line2D([0], [0], lw=0, marker=markers[dt], color='black', label=r'$dt=%.0e$' % dt), time_step_list))
		lines1 += [Line2D([0], [0], lw=0, label=''), Line2D([0], [0], color='black', linestyle='--', label=r'$1$')]
		lines = [Line2D([0], [0], lw=0, marker='s', color='black', label=r'$\tilde{\nu}_r \Delta t^*$')]
		# ax1.legend(handles=lines1)
		ax.legend(handles=lines)

		leg = plt.subplot(gs[:, 1])

	else:
		gs = GridSpec(1, 2, width_ratios=[10, 1])
		ax = plt.subplot(gs[0])
		leg = plt.subplot(gs[:, 1])

	ax.set_xlabel(r'$\tilde{\nu}_r\Delta t$')
	ax.set_ylabel(r'$%s^T/%s^L(\frac{r}{a}=%.3e, \Delta t)$' % (C, C, ra[dirs[0]]))

	leg.axis('off')
	lines0 = list(map(lambda d: Line2D([0], [0], color=colors[d], lw=2, label=r'$\tilde{\nu}_r = %.0e$' % d), sorted(np.append(dr_list, dr0_list))))
	if not('INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX'])):
		lines0 += [Line2D([0], [0], lw=0, label='')]
		lines0 += list(map(lambda dt: Line2D([0], [0], marker=markers[dt], color='black', label=r'$dt=%.0e$' % dt), time_step_list))
	leg.legend(handles=lines0, loc='center')

	for dir in dirs:
		ax.semilogx(dr[dir]*time_step[dir]*period_dump[dir]*dt_list, list(map(lambda t: CT[(dir, t)]/CL[(dir, t)], dt_list)), color=colors[dr[dir]], marker='' if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']) else markers[time_step[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
		if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']):
			ax.scatter(drdtmax[dir]*dr[dir]*time_step[dir]*period_dump[dir], CT[(dir, drdtmax[dir])]/CL[(dir, drdtmax[dir])], color=colors[dr[dir]], marker='s', label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
			ax1.scatter(1/dr[dir], CT[(dir, drdtmax[dir])]/CL[(dir, drdtmax[dir])], color=colors[dr[dir]], marker=markers[time_step[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])

	if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']) and 'FIT' in os.environ and eval(os.environ['FIT']):
		dirs_to_fit = sorted([dir for dir in dirs if dr[dir] in dr_list and dr[dir] > dr_c], key=lambda dir: 1/dr[dir])
		slope, intercept, stderr = operator.itemgetter(0, 1, 4)(st.linregress(list(map(lambda dir: np.log(1/dr[dir]), dirs_to_fit)), list(map(lambda dir: CT[(dir, drdtmax[dir])]/CL[(dir, drdtmax[dir])], dirs_to_fit))))
		ax1.plot(list(map(lambda dir: 1/dr[dir], dirs_to_fit)), list(map(lambda dir: slope*np.log(1/dr[dir]) + intercept, dirs_to_fit)), linestyle='-.', color='black')
		lines1 += [Line2D([0], [0], color='black', linestyle='-.', label='slope = ' + r'$%.1e \pm %.0e$' % (slope, stderr))]
	if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']):
		ax1.legend(handles=lines1)

	# RAW CL AND CT

	fig0 = plt.figure()
	fig0.suptitle(r'$\phi=%1.2f, \tilde{v}=%.2e, N=%.2e$' % (density, vzero, number)+ '\n' + r'$S_{init} = %.2e, S_{max}=%.2e, N_{cases}=%.2e$' % (init_frame, N_cases, int_max))

	if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']):
		gs0 = GridSpec(2, 3, width_ratios=[5, 5, 1])
		axL = plt.subplot(gs0[0, 0])
		axT = plt.subplot(gs0[1, 0])

		axL1 = plt.subplot(gs0[0, 1])
		axL1.set_xlabel(r'$\tau_r \equiv \tilde{\nu}_r^{-1}$')
		axL1.set_ylabel(r'$%s^L(\frac{r}{a}=%.3e, \Delta t^*)$' % (C, ra[dirs[0]]))
		# axL1.plot(np.linspace(1/np.max(np.append(dr_list, dr0_list)), 1/np.min(np.append(dr_list, dr0_list)), 100), [1]*100, color='black', linestyle='--')
		axL1.set_xscale('log')
		axT1 = plt.subplot(gs0[1, 1])
		axT1.set_xlabel(r'$\tau_r \equiv \tilde{\nu}_r^{-1}$')
		axT1.set_ylabel(r'$%s^T(\frac{r}{a}=%.3e, \Delta t^*)$' % (C, ra[dirs[0]]))
		# axT1.plot(np.linspace(1/np.max(np.append(dr_list, dr0_list)), 1/np.min(np.append(dr_list, dr0_list)), 100), [1]*100, color='black', linestyle='--')
		axT1.set_xscale('log')

		lines1 = list(map(lambda dt: Line2D([0], [0], lw=0, marker=markers[dt], color='black', label=r'$dt=%.0e$' % dt), time_step_list))
		# lines1 += [Line2D([0], [0], lw=0, label=''), Line2D([0], [0], color='black', linestyle='--', label=r'$1$')]
		lines = [Line2D([0], [0], lw=0, marker='s', color='black', label=r'$\tilde{\nu}_r \Delta t^*$')]
		# ax1.legend(handles=lines1)
		axL.legend(handles=lines)
		axT.legend(handles=lines)

		leg = plt.subplot(gs0[:, 2])

	else:
		gs0 = GridSpec(1, 3, width_ratios=[5, 5, 1])
		axL = plt.subplot(gs0[0])
		axT = plt.subplot(gs0[1])
		leg = plt.subplot(gs0[2])

	axL.set_xlabel(r'$\tilde{\nu}_r\Delta t$')
	axL.set_ylabel(r'$%s^L(\frac{r}{a}=%.3e, \Delta t)$' % (C, ra[dirs[0]]))
	axT.set_xlabel(r'$\tilde{\nu}_r\Delta t$')
	axT.set_ylabel(r'$%s^T(\frac{r}{a}=%.3e, \Delta t)$' % (C, ra[dirs[0]]))

	leg.axis('off')
	lines0 = list(map(lambda d: Line2D([0], [0], color=colors[d], lw=2, label=r'$\tilde{\nu}_r = %.0e$' % d), sorted(np.append(dr_list, dr0_list))))
	if not('INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX'])):
		lines0 += [Line2D([0], [0], lw=0, label='')]
		lines0 += list(map(lambda dt: Line2D([0], [0], marker=markers[dt], color='black', label=r'$dt=%.0e$' % dt), time_step_list))
	leg.legend(handles=lines0, loc='center')

	for dir in dirs:
		axL.semilogx(dr[dir]*time_step[dir]*period_dump[dir]*dt_list, list(map(lambda t: CL[(dir, t)], dt_list)), color=colors[dr[dir]], marker='' if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']) else markers[time_step[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
		axT.semilogx(dr[dir]*time_step[dir]*period_dump[dir]*dt_list, list(map(lambda t: CT[(dir, t)], dt_list)), color=colors[dr[dir]], marker='' if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']) else markers[time_step[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
		if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']):
			axL.scatter(drdtmax[dir]*dr[dir]*time_step[dir]*period_dump[dir], CL[(dir, drdtmax[dir])], color=colors[dr[dir]], marker='s', label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
			axL1.scatter(1/dr[dir], CL[(dir, drdtmax[dir])], color=colors[dr[dir]], marker=markers[time_step[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
			axT.scatter(drdtmax[dir]*dr[dir]*time_step[dir]*period_dump[dir], CT[(dir, drdtmax[dir])], color=colors[dr[dir]], marker='s', label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
			axT1.scatter(1/dr[dir], CT[(dir, drdtmax[dir])], color=colors[dr[dir]], marker=markers[time_step[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])

	linesL1, linesT1 = [], []
	# if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']) and 'FIT' in os.environ and eval(os.environ['FIT']):
	# 	dirs_to_fit = sorted([dir for dir in dirs if dr[dir] in dr_list and dr[dir] > dr_c], key=lambda dir: 1/dr[dir])
	# 	slopeL, interceptL, stderrL = operator.itemgetter(0, 1, 4)(st.linregress(list(map(lambda dir: np.log(1/dr[dir]), dirs_to_fit)), list(map(lambda dir: CL[(dir, drdtmax[dir])], dirs_to_fit))))
	# 	axL1.plot(list(map(lambda dir: 1/dr[dir], dirs_to_fit)), list(map(lambda dir: slopeL*np.log(1/dr[dir]) + interceptL, dirs_to_fit)), linestyle='-.', color='black')
	# 	linesL1 = [Line2D([0], [0], color='black', linestyle='-.', label='slope = ' + r'$%.1e \pm %.0e$' % (slopeL, stderrL))]
	# 	slopeT, interceptT, stderrT = operator.itemgetter(0, 1, 4)(st.linregress(list(map(lambda dir: np.log(1/dr[dir]), dirs_to_fit)), list(map(lambda dir: CT[(dir, drdtmax[dir])], dirs_to_fit))))
	# 	axT1.plot(list(map(lambda dir: 1/dr[dir], dirs_to_fit)), list(map(lambda dir: slopeT*np.log(1/dr[dir]) + interceptT, dirs_to_fit)), linestyle='-.', color='black')
	# 	linesT1 = [Line2D([0], [0], color='black', linestyle='-.', label='slope = ' + r'$%.1e \pm %.0e$' % (slopeT, stderrT))]
	if 'INTCUUMAX' in os.environ and eval(os.environ['INTCUUMAX']):
		axL1.legend(handles=lines1 + linesL1)
		axT1.legend(handles=lines1 + linesT1)

plot(C, CL, CT)

plt.show()
