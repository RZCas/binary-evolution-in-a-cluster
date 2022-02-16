import os
import glob

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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
figure = pyplot.figure() #(figsize=(12, 12))
plot_Q = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$Q/a$', fontsize=16)
# pyplot.yscale('log')
# pyplot.text(1.5, 0.75, '$e='+str(0.999)+'$', fontsize=16)

root_dir = "output/cluster storage/a_dependence_7/"
for filepath in glob.iglob(root_dir + '**/*.txt', recursive=True):
	color = 'k'
	t_array = []
	r_p_array = []
	lineNumber = 0
	with open(filepath) as f:
		for line in f:
			lineNumber+=1
			data = line.split()
			if isfloat(data[0]) and isfloat(data[1]):
				t = float(data[0])
				a = float(data[7])
			elif data[0] == 'perturber:' and lineNumber>2:
				r_p = float(data[2])
				t_array.append(t)
				r_p_array.append(r_p/a)
			elif data[1] == 'destroyed': color = 'r'
			elif data[1] == 'merger': color = 'g'
	plot_Q.plot(t_array, r_p_array, color)

# plot_ecc = figure.add_subplot(3,2,2)
# ax = pyplot.gca()
# ax.minorticks_on()
# ax.tick_params(labelsize=14)
# ax.set_xlabel(r'$t$ [yr]', fontsize=16)
# ax.set_ylabel(r'$e$', fontsize=16)
# plot_ecc.plot(t, ecc, 'k')

# plot_inc = figure.add_subplot(3,2,3)
# ax = pyplot.gca()
# ax.minorticks_on()
# ax.tick_params(labelsize=14)
# ax.set_xlabel(r'$t$ [yr]', fontsize=16)
# ax.set_ylabel(r'$i$', fontsize=16)
# plot_inc.plot(t, inc, 'k')

# plot_long_asc = figure.add_subplot(3,2,4)
# ax = pyplot.gca()
# ax.minorticks_on()
# ax.tick_params(labelsize=14)
# ax.set_xlabel(r'$t$ [yr]', fontsize=16)
# ax.set_ylabel(r'$\Omega$', fontsize=16)
# plot_long_asc.plot(t, long_asc, 'k')

# plot_arg_peri = figure.add_subplot(3,2,5)
# ax = pyplot.gca()
# ax.minorticks_on()
# ax.tick_params(labelsize=14)
# ax.set_xlabel(r'$t$ [yr]', fontsize=16)
# ax.set_ylabel(r'$\omega$', fontsize=16)
# plot_arg_peri.plot(t, arg_peri, 'k')

pyplot.tight_layout()
pyplot.savefig(root_dir+"r_p.pdf")