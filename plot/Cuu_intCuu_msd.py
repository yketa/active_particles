#! /home/yketa/miniconda3/bin/python3.6

import pickle
import matplotlib.pyplot as plt
import numpy as np
import operator
import os
import sys
sys.path.append('/home/yketa')
from exponents import *
from matplotlib.lines import Line2D
from collections import OrderedDict
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy import stats as st

density = float(eval(os.environ['DENSITY'])) if 'DENSITY' in os.environ else 0.8
number = int(eval(os.environ['NUMBER'])) if 'NUMBER' in os.environ else int(1e5)
vzero = float(eval(os.environ['VZERO'])) if 'VZERO' in os.environ else 1e-2

r_min = float(eval(os.environ['R_MIN'])) if 'R_MIN' in os.environ else 1
r_max = float(eval(os.environ['R_MAX'])) if 'R_MAX' in os.environ else 20

dr_min = float(eval(os.environ['DR_MIN'])) if 'DR_MIN' in os.environ else 1e-10
dr_max = float(eval(os.environ['DR_MAX'])) if 'DR_MAX' in os.environ else 1e10
dr_c = float(eval(os.environ['DR_C'])) if 'DR_C' in os.environ else 3e-4

axis = os.environ['AXIS'] if 'AXIS' in os.environ else 'LOGLOG'
fplot = lambda ax: ax.plot if axis == 'LINLIN' else ax.semilogx if axis == 'LOGLIN' else ax.semilogy if axis == 'LINLOG' else ax.loglog

drdtmax_xs = os.environ['DRDTMAX_XSCALE'] if 'DRDTMAX_XSCALE' in os.environ else 'log'
drdtmax_ys = os.environ['DRDTMAX_YSCALE'] if 'DRDTMAX_YSCALE' in os.environ else 'log'
chimax_xs = os.environ['CHIMAX_XSCALE'] if 'CHIMAX_XSCALE' in os.environ else 'log'
chimax_ys = os.environ['CHIMAX_YSCALE'] if 'CHIMAX_YSCALE' in os.environ else 'log'
msd_xs = os.environ['MSD_XSCALE'] if 'MSD_XSCALE' in os.environ else 'log'
msd_ys = os.environ['MSD_YSCALE'] if 'MSD_YSCALE' in os.environ else 'log'

colormap = os.environ['COLORMAP'] if 'COLORMAP' in os.environ else 'jet'

dirs = [dir for dir in os.listdir() if str('D%s_V%s' % tuple(map(float_to_letters, [density, vzero]))) in dir and str('N%s_Ll' % float_to_letters(number)) in dir]
# dirs = ['Dk8000_Vj1000_Rg3000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rg2000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rj1000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rk5000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rf5000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rg1000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rj5000_Nq1000_Ll1000', 'Dk8000_Vj1000_Rk1000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rk1000_Nq1000_Ll1000', 'Dk8000_Vj1000_Rj5000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rk1000_Nq1000_Ll3000', 'Dk8000_Vj1000_Rg5000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rg4000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rk1000_Nq1000_Ll2000', 'Dk8000_Vj1000_Ri5000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rh1000_Nq1000_Ll1000', 'Dk8000_Vj1000_Rh1000_Nq1000_Ll0000', 'Dk8000_Vj1000_Ri1000_Nq1000_Ll1000', 'Dk8000_Vj1000_Rh5000_Nq1000_Ll0000', 'Dk8000_Vj1000_Rh5000_Nq1000_Ll1000', 'Dk8000_Vj1000_Ri1000_Nq1000_Ll0000']
dt_list = ['l1000', 'l2000', 'l4000', 'l5000', 'm1000', 'm2000', 'm4000', 'm5000', 'n1000', 'n2000', 'n4000', 'n5000', 'o1000', 'o2000', 'o4000']
# dt_list = list(map(float_to_letters, [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 16, 19, 20, 22, 27, 32, 38, 40, 45, 50, 54, 64, 77, 91, 100, 109, 129, 154, 183, 200, 218, 260, 309, 368, 400, 438, 500, 521, 620, 738, 879, 1000, 1045, 1244, 1480, 1761, 2000, 2096, 2494, 2967, 3531, 4000, 4201, 4999]))
init_list = [int(eval(os.environ['SINGLE_MSD']))] if 'SINGLE_MSD' in os.environ else [500, 1000, 2000, 5000] # initial frames in MSD calculation
# init_list = [5000]

# default_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
# colors = {dr_list[i]:default_colors[i] for i in range(len(dr_list))}
default_markers = ['.', 's', '*', 'o']
default_linestyles = ['-', '--', '-.', ':']
linestyles = {init_list[i]:default_linestyles[i] for i in range(len(init_list))}

init_frame = float_to_letters(int(eval(os.environ['INITIAL_FRAME'])) if 'INITIAL_FRAME' in os.environ else 5000)
N_cases = int(eval(os.environ['N_CASES'])) if 'N_CASES' in os.environ else 500
int_max = float_to_letters(int(eval(os.environ['INTERVAL_MAXIMUM'])) if 'INTERVAL_MAXIMUM' in os.environ else 100) # maximum number of intervals for Cuu calculation
fname = lambda var, dir, t: str('%s/%s_%s_I%s_T%s_M%s_C%s.pickle' % (dir, var, dir[:-7], init_frame, t, int_max, float_to_letters(N_cases)))

snap_max = float_to_letters(int(eval(os.environ['SNAP_MAXIMUM'])) if 'SNAP_MAXIMUM' in os.environ else 1) # maximum number of intervals for MSD calculation
snap_per = float_to_letters(int(eval(os.environ['SNAP_PERIOD'])) if 'SNAP_PERIOD' in os.environ else 50) # period of intervals to take to calculate MSD

def cut(arr, r_min, r_max=-1):
	r_max = arr[-1, 0] if r_max < 0 else r_max
	ind_min = next(ind for ind in range(len(arr)) if arr[ind, 0] >= r_min)
	ind_max = next(ind for ind in range(len(arr)) if arr[ind, 0] >= r_max)
	return np.reshape(np.append(np.append([[r_min, arr[ind_min - 1, 1] + (r_min - arr[ind_min - 1, 0])*(arr[ind_min, 1] - arr[ind_min - 1, 1])/(arr[ind_min, 0] - arr[ind_min - 1, 0])]], arr[ind_min:ind_max]), [[r_max, arr[ind_max - 1, 1] + (r_max - arr[ind_max - 1, 0])*(arr[ind_max, 1] - arr[ind_max - 1, 1])/(arr[ind_max, 0] - arr[ind_max - 1, 0])]]), (ind_max - ind_min + 2, 2))

N, box_size, time_step, period_dump, dr, Cuu, Cuu2D, a = {}, {}, {}, {}, {}, {}, {}, {}
for dir in dirs:
	with open(dir + '/param.pickle', 'rb') as param_file:
		N[dir], box_size[dir], time_step[dir], period_dump[dir], dr[dir], a[dir] = operator.itemgetter(0, 5, 13, 15, 10, 1)(pickle.load(param_file))
	for t in dt_list:
		with open(fname('Cwwb', dir, t), 'rb') as Cuu_file, open(fname('Cnnb', dir, t), 'rb') as Cnn_file:
			Cuu1D, Cuu2D[(dir, t)] = operator.itemgetter(2, 0)(pickle.load(Cuu_file))

			if 'CUT' in os.environ and eval(os.environ['CUT']):
				Cuu[(dir, t)] = cut(Cuu1D, r_min, r_max)
			else:
				Cuu[(dir, t)] = cut(Cuu1D, a[dir])

			Cuu2D[(dir, t)] /= pickle.load(Cnn_file)[0]

dr_list = sorted(list(OrderedDict.fromkeys([dr[dir] for dir in dr if dr[dir] >= dr_min and dr[dir] <= dr_max])))
dr0_list = sorted(list(OrderedDict.fromkeys([dr[dir] for dir in dr if dr[dir] < dr_min or dr[dir] > dr_max])))
time_step_list = sorted(list(OrderedDict.fromkeys([time_step[dir] for dir in time_step])))

# colors = {dr_list[i]:plt.cm.RdYlBu(np.linspace(0, 1, len(dr_list))[i]) for i in range(len(dr_list))}
jet = cm = plt.get_cmap(colormap) 
cNorm  = colors.Normalize(vmin=0, vmax=len(dr_list) - 1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
jet0 = cm0 = plt.get_cmap('Greys')
cNorm0 = colors.Normalize(vmin=0, vmax=len(dr0_list))
scalarMap0 = cmx.ScalarMappable(norm=cNorm0, cmap=jet0)
colors = {**{dr_list[i]:scalarMap.to_rgba(i) for i in range(len(dr_list))}, **{dr0_list[i]:scalarMap0.to_rgba(i + 1) for i in range(len(dr0_list))}}

markers = {time_step_list[i]:default_markers[i] for i in range(len(time_step_list))}

# CHI, CHIMAX, DRDTMAX

fig = plt.figure()
gs = GridSpec(2, 3, width_ratios=[5, 5, 1])
ax0 = plt.subplot(gs[:, 0])
ax1 = plt.subplot(gs[0, 1])
ax2 = plt.subplot(gs[1, 1], sharex=ax1)
plt.setp(ax1.get_xticklabels(), visible=False)
ax3 = plt.subplot(gs[:, 2])
ax3.axis('off')

ax0.set_xlabel((r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t$')
ax0.set_ylabel(r'$\chi(\Delta t) = \frac{N}{L^2}$' + (r'$ \int_{r=r_{min}}^{r=r_{max}} dr$' if 'CUT' in os.environ and eval(os.environ['CUT']) else r'$\int_{r=a}^{r=L/2} dr$') + ' '+ r'$2 \pi r C_{uu}(r, \Delta t)$')

ax1.set_ylabel((r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t^*$')
ax1.set_xscale(drdtmax_xs)
ax1.set_yscale(drdtmax_ys)

ax2.set_xlabel(r'$\tau_r \equiv \tilde{\nu}_r^{-1}$')
ax2.set_ylabel(r'$\chi(\Delta t^*) = \frac{N}{L^2}$' + (r'$ \int_{r=r_{min}}^{r=r_{max}} dr$' if 'CUT' in os.environ and eval(os.environ['CUT']) else r'$\int_{r=a}^{r=L/2} dr$') + ' '+ r'$2 \pi r C_{uu}(r, \Delta t^*)$')
ax2.set_xscale(chimax_xs)
ax2.set_yscale(chimax_ys)

exploitables = [[], []]
intCuu, drdtmax, intCuumax = {}, {}, {}
for dir in dirs:
	# intCuu[dir] = np.array(list(map(lambda t: np.trapz(Cuu[(dir, t)][:, 1], Cuu[(dir, t)][:, 0]), dt_list)))
	intCuu[dir] = np.array(list(map(lambda t: 2*np.pi*(N[dir]/(box_size[dir]**2))*np.trapz(Cuu[(dir, t)][:, 1]*Cuu[(dir, t)][:, 0], Cuu[(dir, t)][:, 0]), dt_list)))
	if (intCuu[dir][np.argmax(intCuu[dir]):] < np.max(intCuu[dir])).any() and dr[dir] in dr_list: # if there exists a point after the max of intCuu which is lower than the max
		exploitables[dr[dir] < dr_c] = np.append(exploitables[dr[dir] < dr_c], dir)
	drdtmax[dir] = (dr[dir] if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else 1)*period_dump[dir]*time_step[dir]*float(eval(letters_to_float(dt_list[np.argmax(intCuu[dir])])))
	intCuumax[dir] = np.max(intCuu[dir])
	with open(str('%s/intCuu_%s.pickle' % (dir, dir[:-7])), 'wb') as dump_file:
		pickle.dump([np.transpose([list(map(lambda t: dr[dir] * period_dump[dir] * time_step[dir] * float(eval(letters_to_float(t))), dt_list)), intCuu[dir]]), float(eval(letters_to_float(dt_list[np.argmax(intCuu[dir])]))), dr[dir]*period_dump[dir]*time_step[dir]*float(eval(letters_to_float(dt_list[np.argmax(intCuu[dir])]))), np.max(intCuu[dir])], dump_file)
	
	if dr[dir] <= dr_max and dr[dir] >= dr_min:
		fplot(ax0)(list(map(lambda t: (dr[dir] if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else 1) * period_dump[dir] * time_step[dir] * float(eval(letters_to_float(t))), dt_list)), intCuu[dir], marker=markers[time_step[dir]], linestyle='-', color=colors[dr[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
		ax1.scatter(1/dr[dir], drdtmax[dir], marker=markers[time_step[dir]], color=colors[dr[dir]])
		ax2.scatter(1/dr[dir], intCuumax[dir], marker=markers[time_step[dir]], color=colors[dr[dir]])

if 'FIT' in os.environ and eval(os.environ['FIT']):
	slopes, intercepts, stderr = {}, {}, {}
	for part in range(2):
		exploitables[part] = sorted(exploitables[part], key=lambda dir: 1/dr[dir])
	
		slopes[(part, ax1)], intercepts[(part, ax1)], stderr[(part, ax1)] = operator.itemgetter(0, 1, 4)(st.linregress(list(map(lambda dir: np.log(1/dr[dir]), exploitables[part])), list(map(lambda dir: np.log(drdtmax[dir]), exploitables[part]))))
		ax1.plot(list(map(lambda dir: 1/dr[dir], exploitables[part])), list(map(lambda x: x**slopes[(part, ax1)]*np.exp(intercepts[(part, ax1)]), list(map(lambda dir: 1/dr[dir], exploitables[part])))), ['--', '-.'][part], color='black', label=(r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t^* \propto (\tilde{\nu}_r^{-1})^{%.1e \pm %.0e}$' % (slopes[(part, ax1)], stderr[(part, ax1)]))

		slopes[(part, ax2)], intercepts[(part, ax2)], stderr[(part, ax2)] = operator.itemgetter(0, 1, 4)(st.linregress(list(map(lambda dir: np.log(1/dr[dir]), exploitables[part])), list(map(lambda dir: np.log(intCuumax[dir]), exploitables[part]))))
		ax2.plot(list(map(lambda dir: 1/dr[dir], exploitables[part])), list(map(lambda x: x**slopes[(part, ax2)]*np.exp(intercepts[(part, ax2)]), list(map(lambda dir: 1/dr[dir], exploitables[part])))), ['--', '-.'][part], color='black', label=r'$\chi(\Delta t^*) \propto (\tilde{\nu}_r^{-1})^{%.1e \pm %.0e}$' % (slopes[(part, ax2)], stderr[(part, ax2)]))

	ax1.legend()
	ax2.legend()

# for dir in list(map(lambda d: 'Dk8000_Vj1000_R%s_Nq1000_Ll0000' % d, ['g2000', 'g7000', 'h2000', 'h7000', 'i2000', 'i7000'])):
# 	print(dir, drdtmax[dir]/(dr[dir]*period_dump[dir]*time_step[dir]))

lines = list(map(lambda d: Line2D([0], [0], color=colors[d], lw=2, label=r'$\tilde{\nu}_r = %.0e$' % d), dr_list))
lines0 = list(map(lambda d: Line2D([0], [0], color=colors[d], lw=2, label=r'$\tilde{\nu}_r = %.0e$' % d), sorted(np.append(dr_list, dr0_list))))
lines += [Line2D([0], [0], lw=0, label='')]
# lines0 += [Line2D([0], [0], lw=0, label='')]
lines += list(map(lambda dt: Line2D([0], [0], lw=0, marker=markers[dt], color='black', label=r'$dt=%.0e$' % dt), time_step_list))
linesdt = list(map(lambda dt: Line2D([0], [0], marker=markers[dt], color='black', label=r'$dt=%.0e$' % dt), time_step_list))
# lines0 += list(map(lambda dt: Line2D([0], [0], marker=markers[dt], color='black', label=r'$dt=%.0e$' % dt), time_step_list))
ax3.legend(handles=lines, loc='center')

fig.subplots_adjust(wspace=0.4)
fig.subplots_adjust(hspace=0.05)
fig.suptitle(r'$N=%.1e, \phi=%1.2f, \tilde{v}=%.1e$' % (number, density, vzero) + (r'$, r_{min}=%.2e, r_{max}=%.2e$' % (r_min, r_max) if 'CUT' in os.environ and eval(os.environ['CUT']) else '') + '\n' + r'$N_{cases}=5\cdot10^2, S_{init}=5\cdot10^3, S_{max}=1\cdot10^2$')

# CHI, CHIMAX, DRDTMAX (integrated from 2D)

if 'TWOD' in os.environ and eval(os.environ['TWOD']):
	fig1 = plt.figure()
	gs1 = GridSpec(2, 3, width_ratios=[5, 5, 1])
	ax10 = plt.subplot(gs[:, 0])
	ax11 = plt.subplot(gs[0, 1])
	ax12 = plt.subplot(gs[1, 1], sharex=ax1)
	plt.setp(ax11.get_xticklabels(), visible=False)
	ax13 = plt.subplot(gs[:, 2])
	ax13.axis('off')

	ax10.set_title('integrated from 2D ' + r'$C_{uu}$')
	ax10.set_xlabel((r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t$')
	ax10.set_ylabel(r'$\chi(\Delta t) = \frac{N}{L^2}$' + (r'$ \int_{r=0}^{r=r_{max}} dr$' if 'CUT' in os.environ and eval(os.environ['CUT']) else r'$\int_{r=0}^{r=L/2} dr$') + ' '+ r'$2 \pi r C_{uu}(r, \Delta t)$')

	ax11.set_ylabel((r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t^*$')
	ax11.set_xscale(drdtmax_xs)
	ax11.set_yscale(drdtmax_ys)

	ax12.set_xlabel(r'$\tau_r \equiv \tilde{\nu}_r^{-1}$')
	ax12.set_ylabel(r'$\chi(\Delta t^*) = \frac{N}{L^2}$' + (r'$ \int_{r=0}^{r=r_{max}} dr$' if 'CUT' in os.environ and eval(os.environ['CUT']) else r'$\int_{r=0}^{r=L/2} dr$') + ' '+ r'$2 \pi r C_{uu}(r, \Delta t^*)$')
	ax12.set_xscale(chimax_xs)
	ax12.set_yscale(chimax_ys)

	exploitables1 = [[], []]
	intCuu1, drdtmax1, intCuumax1 = {}, {}, {}
	cases_for_int = [(i, j) for i in range(-int(N_cases/2), int(N_cases/2)) for j in range(-int(N_cases/2), int(N_cases/2)) if i**2 + j**2 <= (((r_max/box_size[dir])*N_cases)**2 if 'CUT' in os.environ and eval(os.environ['CUT']) else (N_cases/2)**2)]
	for dir in dirs:
		intCuu1[dir] = np.array(list(map(lambda t: (N[dir]/(N_cases**2))*np.sum(list(map(lambda case: Cuu2D[(dir, t)][case], cases_for_int))), dt_list)))
		if (intCuu1[dir][np.argmax(intCuu1[dir]):] < np.max(intCuu1[dir])).any() and dr[dir] in dr_list: # if there exists a point after the max of intCuu which is lower than the max
			exploitables1[dr[dir] < dr_c] = np.append(exploitables1[dr[dir] < dr_c], dir)
		drdtmax1[dir] = (dr[dir] if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else 1)*period_dump[dir]*time_step[dir]*float(eval(letters_to_float(dt_list[np.argmax(intCuu1[dir])])))
		intCuumax1[dir] = np.max(intCuu1[dir])
		
		if dr[dir] <= dr_max and dr[dir] >= dr_min:
			fplot(ax10)(list(map(lambda t: (dr[dir] if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else 1) * period_dump[dir] * time_step[dir] * float(eval(letters_to_float(t))), dt_list)), intCuu1[dir], marker=markers[time_step[dir]], linestyle='-', color=colors[dr[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
			ax11.scatter(1/dr[dir], drdtmax1[dir], marker=markers[time_step[dir]], color=colors[dr[dir]])
			ax12.scatter(1/dr[dir], intCuumax1[dir], marker=markers[time_step[dir]], color=colors[dr[dir]])

	if 'FIT' in os.environ and eval(os.environ['FIT']):
		slopes1, intercepts1, stderr1 = {}, {}, {}
		for part in range(2):
			exploitables1[part] = sorted(exploitables1[part], key=lambda dir: 1/dr[dir])
		
			slopes1[(part, ax11)], intercepts1[(part, ax11)], stderr1[(part, ax11)] = operator.itemgetter(0, 1, 4)(st.linregress(list(map(lambda dir: np.log(1/dr[dir]), exploitables1[part])), list(map(lambda dir: np.log(drdtmax1[dir]), exploitables1[part]))))
			ax11.plot(list(map(lambda dir: 1/dr[dir], exploitables1[part])), list(map(lambda x: x**slopes1[(part, ax11)]*np.exp(intercepts1[(part, ax11)]), list(map(lambda dir: 1/dr[dir], exploitables1[part])))), ['--', '-.'][part], color='black', label=(r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t^* \propto (\tilde{\nu}_r^{-1})^{%.1e \pm %.0e}$' % (slopes1[(part, ax11)], stderr1[(part, ax11)]))

			slopes1[(part, ax12)], intercepts1[(part, ax12)], stderr1[(part, ax12)] = operator.itemgetter(0, 1, 4)(st.linregress(list(map(lambda dir: np.log(1/dr[dir]), exploitables1[part])), list(map(lambda dir: np.log(intCuumax1[dir]), exploitables1[part]))))
			ax12.plot(list(map(lambda dir: 1/dr[dir], exploitables1[part])), list(map(lambda x: x**slopes1[(part, ax12)]*np.exp(intercepts1[(part, ax12)]), list(map(lambda dir: 1/dr[dir], exploitables1[part])))), ['--', '-.'][part], color='black', label=r'$\chi(\Delta t^*) \propto (\tilde{\nu}_r^{-1})^{%.1e \pm %.0e}$' % (slopes1[(part, ax12)], stderr1[(part, ax12)]))

		ax11.legend()
		ax12.legend()

	ax13.legend(handles=lines, loc='center')

	fig1.subplots_adjust(wspace=0.4)
	fig1.subplots_adjust(hspace=0.05)
	fig1.suptitle(r'$N=%1.e, \phi=%1.2f, \tilde{v}=%.1e$' % (number, density, vzero) + (r'$, r_{min}=%.2e, r_{max}=%.2e$' % (r_min, r_max) if 'CUT' in os.environ and eval(os.environ['CUT']) else '') + '\n' + r'$N_{cases}=5\cdot10^2, S_{init}=5\cdot10^3, S_{max}=1\cdot10^2$')

# CHI, MSD

fig0 = plt.figure()
gs0 = GridSpec(1, 3, width_ratios=[5, 5, 1])
ax00 = plt.subplot(gs0[0])
ax01 = plt.subplot(gs0[1])
ax02 = plt.subplot(gs0[2])
ax02.axis('off')

ax00.set_title(r'$N_{cases}=5\cdot10^2, S_{init}=5\cdot10^3, S_{max}=1\cdot10^2$')
ax00.set_xlabel((r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t$')
ax00.set_ylabel(r'$\chi(\Delta t) = \frac{N}{L^2}$' + (r'$ \int_{r=r_{min}}^{r=r_{max}} dr$' if 'CUT' in os.environ and eval(os.environ['CUT']) else r'$\int_{r=a}^{r=L/2} dr$') + ' '+ r'$2 \pi r C_{uu}(r, \Delta t)$')

ax01.set_title((r'$S_{init}=%.1e, $' % int(eval(os.environ['SINGLE_MSD'])) if 'SINGLE_MSD' in os.environ else '') + r'$S_{max}=%1.e$' % float(eval(letters_to_float(snap_max))))
ax01.set_xscale(msd_xs)
ax01.set_yscale(msd_ys)
ax01.set_xlabel((r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t$')
ax01.set_ylabel(r'$<|\Delta r(\Delta t)|^2>/\Delta t$')

msd_fname = lambda dir, init: str('%s/msd_sterr_%s_I%s_M%s_P%s.csv' % (dir, dir[:-7], float_to_letters(init), snap_max, snap_per))

msd = {}
for dir in dirs:
	fplot(ax00)(list(map(lambda t: (dr[dir] if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else 1) * period_dump[dir] * time_step[dir] * float(eval(letters_to_float(t))), dt_list)), intCuu[dir], marker=markers[time_step[dir]], linestyle='-', color=colors[dr[dir]], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])
	for init in init_list:
		msd[(dir, init)] = np.genfromtxt(fname=msd_fname(dir, init), delimiter=',', skip_header=True)
		ax01.errorbar(msd[(dir, init)][:, 0]*(dr[dir] if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else 1), msd[(dir, init)][:, 1]/msd[(dir, init)][:, 0], yerr=msd[(dir, init)][:, 2]/msd[(dir, init)][:, 0], color=colors[dr[dir]], linestyle=linestyles[init], label=r'$\tilde{\nu}_r = %.0e$' % dr[dir])

if not('SINGLE_MSD' in os.environ):
	lines1 = list(map(lambda init: Line2D([0], [0], color='black', lw=2, linestyle=linestyles[init], label=r'$S_{init} = %.0e$' % init), init_list))
	ax01.legend(handles=lines1)
ax02.legend(handles=lines0, loc='center')
ax00.legend(handles=linesdt)

fig0.subplots_adjust(wspace=0.4)
fig0.subplots_adjust(hspace=0.05)
fig0.suptitle(r'$N=%.1e, \phi=%1.2f, \tilde{v}=%.1e$' % (number, density, vzero) + (r'$, r_{min}=%.2e, r_{max}=%.2e$' % (r_min, r_max) if 'CUT' in os.environ and eval(os.environ['CUT']) else ''))

plt.show()

