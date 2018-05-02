#! /home/yketa/miniconda3/bin/python3.6

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider

import numpy as np
import os
import pickle
import sys
sys.path.append('/home/yketa')
from exponents import *

font_size = int(eval(os.environ['FONT_SIZE'])) if 'FONT_SIZE' in os.environ else 10
marker_size = int(eval(os.environ['MARKER_SIZE'])) if 'MARKER_SIZE' in os.environ else 20
mpl.rcParams.update({'font.size': font_size, 'lines.markersize': marker_size})

data_dir = os.environ['DATA_DIRECTORY'] if 'DATA_DIRECTORY' in os.environ else os.getcwd() # data directory

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

with open(data_dir + '/param.pickle', 'rb') as param_file:
	N, a, pdi, N_sizes, density, box_size, kT, mu, k, vzero, dr, damp_bro, shear_rate, time_step, N_steps, period_dump, prep_steps = pickle.load(param_file)

# C44

Css = {} # strain correlation grids
for file in files:
	with open(file, 'rb') as Css_file:
		Css[file] = pickle.load(Css_file)[1]

sep = box_size/Ncases # distance corresponding to a grid box length
def value(r, ang, file):
    # return Css[t] value at (r, ang) in polar coordinates
    point = r*np.array([np.cos(ang), np.sin(ang)])
    index = (point//sep + Ncases)%Ncases
    return Css[file][tuple(list(map(int, index)))]

X = np.linspace(r_min, r_max, points_x) # list of x coordinates
theta = np.linspace(0, 2*np.pi, points_theta) # list of angles
integ = lambda r, file: np.trapz(np.array(list(map(lambda ang: value(r, ang, file), theta)))*np.cos(4*theta), theta)/np.pi # C44 at point r for Css[file]

C44 = {} # C44 at X
for file in files:
	C44[file] = list(map(lambda r: integ(r, file), X))

# PLOT

class LineBuilder: # adjustable powerlaw line
	def __init__(self, line, slope, x0, y0):
		self.line = line
		self.slope = slope # slope of the line
		self.x0 = x0 # x coordinate through which line passes
		self.y0 = y0 # y coordinate through which line passes
		self.xs = np.array(line.get_xdata())
		self.ys = np.array(line.get_ydata())
		self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
	def __call__(self, event):
		if event.inaxes!=self.line.axes: return
		self.x0 = event.xdata
		self.y0 = event.ydata
		self.ys = self.y0/(self.x0**self.slope)*(self.xs**self.slope) # line passes through clicked point
		self.line.set_data(self.xs, self.ys)
		self.line.figure.canvas.draw()

fig, ax = plt.subplots()
fig.set_size_inches(30, 30)

ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$C_4^4(r) = \frac{1}{\pi}\int_0^{2\pi}d\theta$' + ' ' + r'$C_{\epsilon_{xy}\epsilon_{xy}}(r, \theta)$' + ' ' + r'$\cos4\theta$')
ax.set_xlim([r_min, r_max])
ax.set_ylim([y_min, y_max])

for file in files:
	ax.loglog(X, C44[file], color=colors[dt(file)], label=r'$\Delta t = %.0e$' % (period_dump*time_step*dt(file)))

legend0 = list(map(lambda t: Line2D([0], [0], color=colors[t], label=r'$\Delta t = %.2e$' % (t*time_step*period_dump)), dt_list))
legend0 += [Line2D([0], [0], lw=0, label='')]
legend_slope = lambda sl: Line2D([0], [0], color='black', linestyle='--', label=r'$C_4^4(r) \propto r^{%.2e}$' % sl) # legend for adjustable powerlaw line
ax.legend(handles=legend0 + [legend_slope(slope0)])

line, = ax.loglog(X, y0/(x0**slope0)*(X**slope0), color='black', linestyle='--') # adjustable powerlaw line
adjustable_line = LineBuilder(line, slope0, x0, y0)

divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="5%", pad=0.6)
sr = Slider(cax, 'slope', slope_min, slope_max, valinit=slope0) # slider
def update(val):
	# change line slope on slider update
	adjustable_line.slope = sr.val # new slope
	adjustable_line.ys = adjustable_line.y0/(adjustable_line.x0**adjustable_line.slope)*(adjustable_line.xs**adjustable_line.slope) # y values of adjustable line
	adjustable_line.line.set_data(adjustable_line.xs, adjustable_line.ys)
	adjustable_line.line.figure.canvas.draw()
	ax.legend(handles=legend0 + [legend_slope(adjustable_line.slope)])
sr.on_changed(update)

title = r'$N=%.2e, \phi=%1.2f, $' % (N, density)
title += r'$\tilde{v}=%.2e, \tilde{\nu}_r=%.2e$' % (vzero, dr) if not('TEMPERATURE' in os.environ and eval(os.environ['TEMPERATURE'])) else r'$kT=%.2e, k=%.2e$' % (kT, k)
title += '\n'
title += r'$S_{init}=%.2e, S_{max}=%.2e, N_{cases}=%.2e, r_{cut}=%.2e, \sigma=%.2e$' % (init_frame, int_max, Ncases, r_cut, sigma)
title += '\n'
title += r'$N_r=%.2e, N_{\theta}=%.2e$' % (points_x, points_theta)
fig.suptitle(title)

plt.show()
