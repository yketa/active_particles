import active_particles.naming as naming
import active_particles.analysis.css as css

from active_particles.init import get_env

import pickle

import matplotlib as mpl
if get_env('SHOW', default=False, vartype=bool):
	mpl.use('Agg')	# avoids crash if launching without display
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
cmap = plt.cm.jet

if __name__ == '__main__':	# executing as script

	# VARIABLE DEFINITION

	css.define_variables()	# define variables used to plot

	# DATA

	with open(css.data_dir + '/' + css.Css_filename, 'rb') as Css_dump_file,\
		open(css.data_dir + '/' + css.Ccc_filename, 'rb') as Ccc_dump_file:
		Sgrid, Css2D = pickle.load(Css_dump_file)
		Cgrid, Ccc2D = pickle.load(Ccc_dump_file)

	# PLOT

	css.plot(Sgrid, Css2D, '\epsilon_{xy}')	# plotting shear strain map and correlation
	css.plot(Cgrid, Ccc2D, '\omega')		# plotting displacement vorticity map and correlation

	if get_env('SHOW', default=False, vartype=bool):
		plt.show()
