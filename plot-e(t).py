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
plot_ecc = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$e$', fontsize=16)
# pyplot.yscale('log')
# pyplot.text(1.5, 0.75, '$e='+str(0.999)+'$', fontsize=16)

root_dir = "output/cluster storage/"
for filepath in glob.iglob(root_dir + '**/*.txt', recursive=True):
	color = 'k'
	t = []
	a = []
	ecc = []
	inc = []
	long_asc = []
	arg_peri = []
	R = []
	z = []
	t_previous = 0
	dt = []
	with open(filepath) as f:
		for line in f:
			data = line.split()
			if isfloat(data[0]) and isfloat(data[1]):
				t.append(float(data[0]))
				dt.append(t[-1]-t_previous)
				t_previous = t[-1]
				R.append(float(data[1]))
				z.append(float(data[2]))
				a.append(float(data[7]))
				ecc.append(float(data[10]))
				inc.append(float(data[11]))
				long_asc.append(float(data[12]))
				arg_peri.append(float(data[13]))
			elif data[1] == 'destroyed': color = 'r'
			elif data[1] == 'merger': color = 'g'
	if color == 'g': plot_ecc.plot(t, ecc, color)

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
pyplot.savefig("output/cluster storage/e(t)-mergers.pdf")