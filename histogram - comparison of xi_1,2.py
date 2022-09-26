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

figure = pyplot.figure(figsize=(9/2, 7/2))
plot_de = figure.add_subplot(1,1,1)

figure.suptitle(r'$M_{\rm cl}=10^6M_\odot$')
output_file = 'output/tidal_effects_strength/mtotal=1e6.pdf'
types = ['wide_range_1','wide_range_mtotal=1e6_nokicks','wide_range_mtotal=1e6_plummer','wide_range_mtotal=1e6_nokicks_plummer']
legends = ['Hernquist', 'Hernquist, no kicks', 'Plummer', 'Plummer, no kicks']
colors = ['black', 'black', 'red', 'red']
linestyles = ['solid', 'dashed', 'solid', 'dashed']

figure.suptitle(r'$M_{\rm cl}=10^5M_\odot$')
output_file = 'output/tidal_effects_strength/mtotal=1e5.pdf'
types = ['wide_range_mtotal=1e5','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e5_plummer','wide_range_mtotal=1e5_plummer_b=1']
legends = ['Hernquist', 'Hernquist, no kicks', 'Plummer', r'Plummer, $b=1\,\mathrm{pc}$']
colors = ['black', 'black', 'red', 'red']
linestyles = ['solid', 'dashed', 'solid', 'dashed']

figure.suptitle(r'Hernquist')
output_file = 'output/tidal_effects_strength/hernquist-norm.pdf'
types = ['wide_range_mtotal=1e5','wide_range_1','wide_range_mtotal=1e7_hernquist']
legends = [r'$M_{\rm cl}=10^5M_\odot$', r'$M_{\rm cl}=10^6M_\odot$', r'$M_{\rm cl}=10^7M_\odot$']
colors = ['black', 'blue', 'red']
linestyles = ['solid', 'solid', 'solid']

for i in range(0, len(types)):
	input_file = 'output/tidal_effects_strength/'+types[i]+'-tidal_effects_strength.txt'
	xi_1 = []
	with open(input_file) as f:
		for line in f:
			data = line.split()
			if isfloat(data[1]):
				xi_1.append(float(data[1]))
	if len(xi_1) == 0: continue
	xi_1 = np.array(xi_1)	
	# normalize
	hist_values, bins = np.histogram(xi_1, bins=np.geomspace(xi_1.min(), xi_1.max(), 10))
	# hist_values = np.array(hist_values)
	hist_values_norm = hist_values/hist_values.sum()
	widths = (bins[1:] - bins[:-1])
	plot_de.bar (bins[:-1], hist_values_norm, widths)
	# plot_de.hist (xi_1, histtype='step', bins=np.geomspace(xi_1.min(), xi_1.max(), 10), log=True, label=legends[i], color=colors[i], linestyle=linestyles[i])

plot_de.legend(loc='upper left')
pyplot.xscale('log')
pyplot.yscale('log')
pyplot.tight_layout()
pyplot.savefig(output_file)
pyplot.clf()