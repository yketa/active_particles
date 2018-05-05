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

font_size = int(eval(os.environ['FONT_SIZE'])) if 'FONT_SIZE' in os.environ else 10
marker_size = int(eval(os.environ['MARKER_SIZE'])) if 'MARKER_SIZE' in os.environ else 20
mpl.rcParams.update({'font.size': font_size, 'lines.markersize': marker_size})

ratio_legend = int(eval(os.environ['RATIO_LEGEND'])) if 'RATIO_LEGEND' in os.environ else 7 # width ratio between legend and figure
ncol_legend = int(eval(os.environ['NCOL_LEGEND'])) if 'NCOL_LEGEND' in os.environ else 1 # number of columns for the legend

wspace = float(eval(os.environ['WSPACE'])) if 'WSPACE' in os.environ else 0.2
hspace = float(eval(os.environ['HSPACE'])) if 'HSPACE' in os.environ else 0.05

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else os.getcwd() # data directory
os.chdir(data_dir) # change working directory to data directory

r_cut = float(eval(os.environ['R_CUT'])) if 'R_CUT' in os.environ else 2 # cut-off radius for Css calculation
sigma = float(eval(os.environ['SIGMA'])) if 'SIGMA' in os.environ else 2 # gaussian standard deviation for Css calculation

Ncases = int(eval(os.environ['N_CASES'])) if 'N_CASES' in os.environ else 1000 # number of cases for Css calculation
init_frame = int(eval(os.environ['INITIAL_FRAME'])) if 'INITIAL_FRAME' in os.environ else 5000 # initial frame in Css calculation
int_max = int(eval(os.environ['INTERVAL_MAXIMUM'])) if 'INTERVAL_MAXIMUM' in os.environ else 1 # maximum number of interval in Css calculation

snip0 = str('Css' + ('b' if not('ENDPOINT' in os.environ and eval(os.environ['ENDPOINT'])) else ''))
snip1 = 'I%s' % float_to_letters(init_frame)
snip2 = 'M%s_C%s_RCUT%s_SIGM%s.pickle' % tuple(map(float_to_letters, [int_max, Ncases, r_cut, sigma]))
files = [file for file in os.listdir() if (snip0 in file) and (snip1 in file) and (snip2 in file)] # Css pickle files
dt = lambda file: eval(letters_to_float(file.split("_T")[1][:5])) # lag time from file name
dt_list = sorted(list(map(dt, files))) # list of lag times

slope0 = float(eval(os.environ['SLOPE'])) if 'SLOPE' in os.environ else -2 # default slope for slider
slope_min = float(eval(os.environ['SLOPE_MIN'])) if 'SLOPE_MIN' in os.environ else -5 # minimum slope for slider
slope_max = float(eval(os.environ['SLOPE_MAX'])) if 'SLOPE_MAX' in os.environ else 0 # maximum slope for slider

colormap = os.environ['COLORMAP'] if 'COLORMAP' in os.environ else 'jet' # colormap for curves
cm = plt.get_cmap(colormap)
cNorm  = colors.Normalize(vmin=0, vmax=len(dt_list) - 1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
colors = {dt_list[i]: scalarMap.to_rgba(i) for i in range(len(dt_list))} # curves colors

y_min = float(eval(os.environ['Y_MIN'])) if 'Y_MIN' in os.environ else 1e-3 # minimum C44 for plotting window
y_max = float(eval(os.environ['Y_MAX'])) if 'Y_MAX' in os.environ else 2e-1 # maximum C44 for plotting window
y0 = np.exp(np.mean([np.log(y_min), np.log(y_max)])) # y coordinate in the middle of the graph
r_min = float(eval(os.environ['R_MIN'])) if 'R_MIN' in os.environ else 10 # minimum raidus for plotting window
r_max = float(eval(os.environ['R_MAX'])) if 'R_MAX' in os.environ else 40 # maximum raidus for plotting window
x0 = np.exp(np.mean([np.log(r_min), np.log(r_max)])) # x coordinate in the middle of the graph

points_x = int(eval(os.environ['POINTS_THETA'])) if 'POINTS_THETA' in os.environ else 100 # number of radii to evaluate
points_theta = int(eval(os.environ['POINTS_THETA'])) if 'POINTS_THETA' in os.environ else 100 # number of angles to evaluate

with open('param.pickle', 'rb') as param_file:
	N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(param_file)
av_p_sep = box_size/np.sqrt(N) # average particle separation

# C44

Css = {} # strain correlation grids
for file in files:
	with open(file, 'rb') as Css_file:
		Css[file] = pickle.load(Css_file)[1]

sep = box_size/Ncases # distance corresponding to a grid box length
def value(r, ang, file):
    # return Css[t] value at (r, ang) in polar coordinates
    point = (r*av_p_sep)*np.array([np.cos(ang), np.sin(ang)])
    index = (point//sep + Ncases)%Ncases
    return Css[file][tuple(list(map(int, index)))]

X = np.linspace(r_min, r_max, points_x) # list of x coordinates
theta = np.linspace(0, 2*np.pi, points_theta) # list of angles
integ = lambda r, file: np.trapz(np.array(list(map(lambda ang: value(r, ang, file), theta)))*np.cos(4*theta), theta)/np.pi # C44 at point r for Css[file]

C44 = {} # C44 at X
for file in files:
	C44[file] = list(map(lambda r: integ(r, file), X))

# PLOT

law_func = lambda law, slope, x, y, X: {'Powerlaw': y/(x**slope)*(X**slope), 'Exponential': y*np.exp(slope*(X - x))}[law] # function corresponding to chosen law ('PWL' = powerlaw, 'EXP' = exponential)

class LineBuilder: # adjustable line
	def __init__(self, line, slope, x0, y0, leg):
		self.line = line
		self.slope = slope # slope of the line
		self.x0 = x0 # x coordinate through which line passes
		self.y0 = y0 # y coordinate through which line passes
		self.leg = leg # legend axis
		self.law = 'Powerlaw' # law of the line
		self.xs = np.array(line.get_xdata())
		self.ys = np.array(line.get_ydata())
		self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
	def __call__(self, event):
		if event.inaxes!=self.line.axes: return
		self.x0 = event.xdata
		self.y0 = event.ydata
		self.draw()
	def draw(self):
		self.ys = law_func(self.law, self.slope, self.x0, self.y0, self.xs) # line passes through clicked point according to law
		self.line.set_data(self.xs, self.ys)
		self.leg.legend(handles=legend0 + [legend_slope(self.slope, self.law)], loc='center', ncol=ncol_legend) # legend according to law
		ax.set_xscale({'Exponential': 'linear', 'Powerlaw': 'log'}[adjustable_line.law]) # change scale of the graph to always have a straight line
		self.line.figure.canvas.draw()

fig = plt.figure()
fig.set_size_inches(30, 30)
fig.subplots_adjust(wspace=wspace)
fig.subplots_adjust(hspace=hspace)

gs = GridSpec(1, 2, width_ratios=[1, 1/ratio_legend])
ax = plt.subplot(gs[0])
leg = plt.subplot(gs[1])
leg.axis('off')

ax.set_xlabel(r'$r/a$' + ' ' + r'$(a = L/\sqrt{N})$')
ax.set_ylabel(r'$C_4^4(r) = \frac{1}{\pi}\int_0^{2\pi}d\theta$' + ' ' + r'$C_{\epsilon_{xy}\epsilon_{xy}}(r, \theta)$' + ' ' + r'$\cos4\theta$')
ax.set_xlim([r_min, r_max])
ax.set_ylim([y_min, y_max])
ax.set_yscale('log')
ax.set_xscale('log')

for file in files:
	ax.loglog(X, C44[file], color=colors[dt(file)], label=r'$\Delta t = %.0e$' % (period_dump*time_step*dt(file)))

if not('TEMPERATURE' in os.environ and eval(os.environ['TEMPERATURE'])):
	legend0 = [mpatches.Patch(color='none', label=r'$nD_0 = \frac{Nv_0^2}{2\nu_r L^2} = %.2e$' % ((N*(vzero**2))/(2*dr*(box_size**2))))]
	legend0 += list(map(lambda t: Line2D([0], [0], color=colors[t], label=r'$nD_0\Delta t = %.2e$' % ((t*time_step*period_dump*N*(vzero**2))/(2*dr*(box_size**2)))), dt_list))
else:
	legend0 = [mpatches.Patch(color='none', label=r'$nD_0 = \frac{2 k_B T N}{\lambda a L^2} = %.2e$' % ((2*kT*N)/(damp_bro*a*(box_size**2))))]
	legend0 += list(map(lambda t: Line2D([0], [0], color=colors[t], label=r'$nD_0\Delta t = %.2e$' % ((2*kT*N*t*time_step*period_dump)/(damp_bro*a*(box_size**2)))), dt_list))
legend0 += [Line2D([0], [0], lw=0, label='')]
legend_slope = lambda sl, law: Line2D([0], [0], color='black', linestyle='--', label=r'$C_4^4(r) \propto$' + {'Powerlaw': r'$r^{%.2e}$' % sl, 'Exponential': r'$e^{%.2er}$' % sl}[law]) # legend for adjustable powerlaw line
leg.legend(handles=legend0 + [legend_slope(slope0, 'Powerlaw')], loc='center', ncol=ncol_legend)

line, = ax.plot(X, y0/(x0**slope0)*(X**slope0), color='black', linestyle='--') # adjustable powerlaw line
adjustable_line = LineBuilder(line, slope0, x0, y0, leg)

divider_ax = make_axes_locatable(ax)
cax = divider_ax.append_axes("bottom", size="5%", pad=0.6)
sr = Slider(cax, 'slope', slope_min, slope_max, valinit=slope0) # slider
def update_slope(val):
	# change line slope on slider update
	adjustable_line.slope = sr.val # new slope
	adjustable_line.draw()
sr.on_changed(update_slope)

divider_leg = make_axes_locatable(leg)
cleg = divider_leg.append_axes("bottom", size="20%", pad=0.6)
radio = RadioButtons(cleg, ('Powerlaw', 'Exponential'), active=0) # radio button
def update_law(law):
	# change fitting law
	adjustable_line.law = law # new law
	adjustable_line.draw()
radio.on_clicked(update_law)

title = r'$N=%.2e, \phi=%1.2f, $' % (N, density)
title += r'$\tilde{v}=%.2e, \tilde{\nu}_r=%.2e$' % (vzero, dr) if not('TEMPERATURE' in os.environ and eval(os.environ['TEMPERATURE'])) else r'$kT=%.2e, k=%.2e$' % (kT, k)
title += '\n'
title += r'$S_{init}=%.2e, S_{max}=%.2e, N_{cases}=%.2e, r_{cut}=%.2e, \sigma=%.2e$' % (init_frame, int_max, Ncases, r_cut, sigma)
title += '\n'
title += r'$N_r=%.2e, N_{\theta}=%.2e$' % (points_x, points_theta)
fig.suptitle(title)

plt.show()
