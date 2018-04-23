#! /home/yketa/miniconda3/bin/python3.6

import os
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/yketa')
from exponents import *
import pickle
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mpcolors
import matplotlib.cm as cmx
import matplotlib as mp
from collections import OrderedDict
import operator

font_size = int(eval(os.environ['FONT_SIZE'])) if 'FONT_SIZE' in os.environ else 20
marker_size = int(eval(os.environ['MARKER_SIZE'])) if 'MARKER_SIZE' in os.environ else 20
mp.rcParams.update({'font.size': font_size, 'lines.markersize': marker_size})

ratio_legend = int(eval(os.environ['RATIO_LEGEND'])) if 'RATIO_LEGEND' in os.environ else 8

wspace = float(eval(os.environ['WSPACE'])) if 'WSPACE' in os.environ else 0.4
hspace = float(eval(os.environ['HSPACE'])) if 'HSPACE' in os.environ else 0.05

systems = [(0.8, 1e-2), (1.0, 1e-1)] # (\phi, \tilde{v}) for each system
tau_rc = {(0.8, 1e-2): 3e3, (1.0, 1e-1): 7e2}
linestyles = {(0.8, 1e-2): '--', (1.0, 1e-1): '-.'}
dr_min = {(0.8, 1e-2): 1e-5, (1.0, 1e-1): 2e-4}
dr_max = {(0.8, 1e-2): 1e-2, (1.0, 1e-1): 2e-2}

colormaps = {systems[system]: [os.environ['COLORMAP0'] if 'COLORMAP0' in os.environ else 'cool', os.environ['COLORMAP0'] if 'COLORMAP0' in os.environ else 'hot'][system] for system in range(len(systems))}

ncol_legend = {systems[system]: [int(eval(os.environ['NCOL_LEGEND'])) if 'NCOL_LEGEND' in os.environ else 1, int(eval(os.environ['NCOL_LEGEND'])) if 'NCOL_LEGEND' in os.environ else 1][system] for system in range(len(systems))}

default_markers = ['.', 's', '*', 'o']

vars = ['Cuu', 'Cee']
C = {'Cuu': 'C_{uu}', 'Cee': 'C_{\hat{u}\hat{u}}'}

dirs, dr_list = [{}]*2
for system in systems:
    dirs[system] = [dir for dir in os.listdir() if 'D%s_V%s' % tuple(map(float_to_letters, system)) in dir and 'Nq1000_Ll' in dir]

r_min_C = float(eval(os.environ['R_MIN'])) if 'R_MIN' in os.environ else 1
r_max_C = float(eval(os.environ['R_MAX'])) if 'R_MAX' in os.environ else 20
init_frame_C = float(eval(os.environ['INITIAL_FRAME_C'])) if 'INITIAL_FRAME_C' in os.environ else 5000
N_cases_C = int(eval(os.environ['N_CASES_C'])) if 'N_CASES_C' in os.environ else 500
int_max_C = int(eval(os.environ['INTERVAL_MAXIMUM_C']))if 'INTERVAL_MAXIMUM_C' in os.environ else 100

init_frame_N = int(eval(os.environ['INITIAL_FRAME_N'])) if 'INITIAL_FRAME_N' in os.environ else 5000
int_max_N = int(eval(os.environ['INTERVAL_MAXIMUM_N'])) if 'INTERVAL_MAXIMUM_N' in os.environ else 1
Ncases_N = int(eval(os.environ['N_CASES_N'])) if 'N_CASES_N' in os.environ else 500
max_box_size_N = float(eval(os.environ['MAX_BOX_SIZE_N'])) if 'MAX_BOX_SIZE_N' in os.environ else 10
Nbins = int(eval(os.environ['N_BINS'])) if 'N_BINS' in os.environ else 200
phimax = lambda system: 1 if system[0] < 1 else 1.5
bins = lambda phimax: np.linspace(0, phimax, Nbins, endpoint=False)

pname = lambda dir: str('%s/param.pickle' % dir)
fname = lambda dir, var, t: str('%s/%sb_%s_I%s_T%s_M%s_C%s.pickle' % ((dir, var, dir[:-7]) + tuple(list(map(float_to_letters, [init_frame_C, t, int_max_C, N_cases_C])))))
iname = lambda dir, var: str('%s/int%sb_%s_I%s_M%s_C%s_RMIN%s_RMAX%s.pickle' % ((dir, var, dir[:-7]) + tuple(list(map(float_to_letters, [init_frame_C, int_max_C, N_cases_C, r_min_C, r_max_C])))))
dname = lambda dir: str('%s/varN_%s' % (dir, dir[:-7]) + '_I%s_M%s_C%s_B%s.pickle' % tuple(map(float_to_letters, [init_frame_N, int_max_N, Ncases_N, max_box_size_N])))

dr_list, dr0_list, time_step_list = {}, {}, {}
box_size, dr, time_step, period_dump, densities, philocmax = {system: {} for system in systems}, {system: {} for system in systems}, {system: {} for system in systems}, {system: {} for system in systems}, {system: {} for system in systems}, {system: {} for system in systems}
dt_max, int_max, CL, CT = {system: {var: {} for var in vars} for system in systems}, {system: {var: {} for var in vars} for system in systems}, {system: {var: {} for var in vars} for system in systems}, {system: {var: {} for var in vars} for system in systems}
all_dr = []

for system in systems:
    for dir in dirs[system]:

        # box_size, dr, time_step, period_dump
        with open(pname(dir), 'rb') as p_file:
            box_size[system][dir], dr[system][dir], time_step[system][dir], period_dump[system][dir] = operator.itemgetter(5, 10, 13, 15)(pickle.load(p_file))

        # densities, philocmax
        with open(dname(dir), 'rb') as d_file:
            densities[system][dir] = pickle.load(d_file)
        histogram = []
        for bin in bins(phimax(system)):
            histogram += [(np.sum(densities[system][dir] >= bin) - np.sum(densities[system][dir] > (bin + phimax(system)/Nbins)))/len(densities[system][dir])]
        philocmax[system][dir] = bins(phimax(system))[np.argmax(histogram)]

        for var in vars:

            # dt_max, int_max
            with open(iname(dir, var), 'rb') as i_file:
                dt_max[system][var][dir], int_max[system][var][dir] = operator.itemgetter(1, 3)(pickle.load(i_file))

            # CL, CT
            with open(fname(dir, var, dt_max[system][var][dir]), 'rb') as f_file:
                CL[system][var][dir], CT[system][var][dir] = operator.itemgetter(3, 4)(pickle.load(f_file))

    dr_list[system] = sorted(list(OrderedDict.fromkeys([dr[system][dir] for dir in dr[system] if dr[system][dir] >= dr_min[system] and dr[system][dir] <= dr_max[system]])))
    dr0_list[system] = sorted(list(OrderedDict.fromkeys([dr[system][dir] for dir in dr[system] if dr[system][dir] < dr_min[system] or dr[system][dir] > dr_max[system]])))
    all_dr += dr_list[system]
    all_dr += dr0_list[system]
    time_step_list[system] = sorted(list(OrderedDict.fromkeys([time_step[system][dir] for dir in dirs[system]])))
time_steps = sorted(list(OrderedDict.fromkeys(np.append(*[time_step_list[system] for system in systems]))))
all_dr = sorted(list(OrderedDict.fromkeys(all_dr)))

# PLOT

fig = plt.figure()
fig.set_size_inches(30, 30)
fig.subplots_adjust(wspace=wspace)
fig.subplots_adjust(hspace=hspace)
fig.suptitle(r'$N=1\cdot10^5$' + '\n' + r'$[C]: S_{init}=%.2e, S_{max}=%.2e, N_{cases}=%.2e, r_{min}=%.2e, r_{max}=%.2e$' % (init_frame_C, int_max_C, N_cases_C, r_min_C, r_max_C) + '\n' + r'$[\phi]: S_{init}=%.2e, S_{max}=%.2e, N_{cases}=%.2e, r_{max}=%.2e$' % (init_frame_N, int_max_N, Ncases_N, max_box_size_N))

gs = GridSpec(4, 4, width_ratios=[1, 1, 2/ratio_legend, 2/ratio_legend])
ax = np.array([[plt.subplot(gs[i, j]) for j in range(2)] for i in range(4)])
plt.setp([a.get_xticklabels() for a in ax[:-1, :].flatten()], visible=False)
ax[0, 0].set_title('Displacement')
ax[0, 1].set_title('Displacement direction')

leg = {systems[system]: plt.subplot(gs[:, 2 + system]) for system in range(len(systems))}
cmap, norm, scalarMap = {system: [] for system in systems}, {system: [] for system in systems}, {system: [] for system in systems}
colors, lines = {}, {}
# markers
markers = {time_steps[i]:default_markers[i] for i in range(len(time_steps))}
for system in systems:

    # colormaps
    cmap[system] += [plt.get_cmap(colormaps[system]), plt.get_cmap('Greys')]
    norm[system] += [mpcolors.Normalize(vmin=0, vmax=len(dr_list[system]) + 1), mpcolors.Normalize(vmin=0, vmax=len(dr0_list[system]))]
    scalarMap[system] = list(map(lambda norm, cmap: cmx.ScalarMappable(norm=norm, cmap=cmap), *[norm[system], cmap[system]]))
    colors[system] = {**{dr_list[system][i]:scalarMap[system][0].to_rgba(i + 1) for i in range(len(dr_list[system]))}, **{dr0_list[system][i]:scalarMap[system][1].to_rgba(i + 1) for i in range(len(dr0_list[system]))}}

    # legends
    lines[system] = [Line2D([0], [0], lw=0, label=r'$\phi=%1.2f, \tilde{v}=%.2f$' % system), Line2D([0], [0], lw=0), Line2D([0], [0], color='black', linestyle=linestyles[system], label=r'$\tau_{r,c} = %.1e$' % tau_rc[system]), Line2D([0], [0], lw=0)]
    lines[system] += list(map(lambda d: Line2D([0], [0], color=colors[system][d], lw=2, label=r'$\tilde{\nu}_r = %.0e$' % d), sorted(np.append(dr_list[system], dr0_list[system]))))
    lines[system] += [Line2D([0], [0], lw=0)]
    lines[system] += list(map(lambda time_step: Line2D([0], [0], color='black', marker=markers[time_step], label=r'$dt = %.0e$' % time_step), time_step_list[system]))
    leg[system].axis('off')
    leg[system].legend(handles=lines[system], loc='center', ncol=ncol_legend[system])

for var in range(len(vars)):

    # xlabel
    ax[-1, var].set_xlabel(r'$\tau_r \equiv \tilde{\nu}_r^{-1}$')
    # ylabel
    ax[0, var].set_ylabel(r'$\phi^*_{loc}$')
    ax[1, var].set_ylabel(r'$\chi(\Delta t^*) = \frac{N}{L^2} \int_{r=r_{min}}^{r=r_{max}} dr %s(r, \Delta t^*)$' % C[vars[var]])
    ax[2, var].set_ylabel((r'$\tilde{\nu}_r$' if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else '') + r'$\Delta t^*$')
    ax[3, var].set_ylabel(r'$%s^T/%s^L(\frac{r}{a} = %.2e, \Delta t^*)$' % (C[vars[var]], C[vars[var]], box_size[systems[var]][dirs[systems[var]][0]]/N_cases_C))
    # xscale
    for line in range(4):
        ax[line, var].set_xscale('log')
        ax[line, var].set_xlim(1/np.max(all_dr), 1/np.min(all_dr))
    # yscale
    ax[0, var].set_yscale('linear')
    ax[1, var].set_yscale('log')
    ax[2, var].set_yscale('log')
    ax[3, var].set_yscale('linear')
    # tau_rc
    for line in range(4):
        for system in systems:
            ax[line, var].axvline(x=tau_rc[system], color='black', linestyle=linestyles[system], clip_on=False)
    # CT/CL = 1
    ax[3, var].axhline(y=1, color='black', linestyle='-', clip_on=False)

    for system in systems:
        for dir in dirs[system]:
            # philocmax
            ax[0, var].scatter(1/dr[system][dir], philocmax[system][dir], color=colors[system][dr[system][dir]], marker=markers[time_step[system][dir]])
            # int_max
            ax[1, var].scatter(1/dr[system][dir], int_max[system][vars[var]][dir], color=colors[system][dr[system][dir]], marker=markers[time_step[system][dir]])
            # dr dt_max
            if dr[system][dir] in dr_list[system]:
                ax[2, var].scatter(1/dr[system][dir], (dr[system][dir] if not('DRDT' in os.environ and not(eval(os.environ['DRDT']))) else 1)*time_step[system][dir]*period_dump[system][dir]*dt_max[system][vars[var]][dir], color=colors[system][dr[system][dir]] if dr[system][dir] in dr_list[system] else 'white', marker=markers[time_step[system][dir]])
            # CLCT
            ax[3, var].scatter(1/dr[system][dir], CT[system][vars[var]][dir]/CL[system][vars[var]][dir], color=colors[system][dr[system][dir]], marker=markers[time_step[system][dir]])

fig.savefig('comparison_plot_Dk8000_Vj1000_Dl1000_Vk1000.eps')
print('Figure saved to \'comparison_plot_Dk8000_Vj1000_Dl1000_Vk1000.eps\'.')
# plt.show()
