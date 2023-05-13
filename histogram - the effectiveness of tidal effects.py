import os
import glob
import numpy as np

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{color}"

# types = ['perpendicular-soft-hernquist', 'perpendicular-hard-hernquist']
# types = ['perpendicular-hard-plummer', 'perpendicular-hard-plummer-light']
# types = ['wide_range_1','wide_range_mtotal=1e5','wide_range_m=60','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e6_nokicks_plummer','wide_range_mtotal=1e6_nokicks','wide_range_mtotal=1e6']
# types = ['wide_range_1','wide_range_mtotal=1e6_nokicks']
# types = ['wide_range_mtotal=1e5','wide_range_mtotal=1e5_nokicks']
# types = ['wide_range_mtotal=1e6_nokicks_plummer']
# types = ['wide_range_mtotal=1e6']

for input_file in glob.glob('output/tidal_effects_strength/*.txt'):
	xi_1 = []
	xi_2 = []
	
	with open(input_file) as f:
		for line in f:
			data = line.split()
			if isfloat(data[1]):
				xi_1.append(float(data[1]))
				xi_2.append(float(data[3]))
	if len(xi_1) == 0: continue
	xi_1 = np.array(xi_1)
	xi_2 = np.array(xi_2)
	figure = pyplot.figure(figsize=(9/2, 7/2))
	figure.suptitle(fr'$\langle\xi_1\rangle = {xi_1.mean():.3f}$, $\langle\xi_2\rangle = {xi_2.mean():.3f}$')
	plot_de = figure.add_subplot(1,1,1)
	plot_de.hist (xi_1, histtype='step', bins=np.geomspace(xi_1.min(), xi_1.max(), 10), log=True, label=r'$\xi_1$')
	plot_de.hist (xi_2, histtype='step', bins=np.geomspace(xi_2.min(), xi_2.max(), 10), log=True, label=r'$\xi_2$')
	plot_de.legend(loc='lower left')
	pyplot.xscale ('log')
	pyplot.tight_layout()
	pyplot.savefig(input_file[:-3]+'pdf')
	pyplot.clf()