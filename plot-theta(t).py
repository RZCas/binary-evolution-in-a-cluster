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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
figure = pyplot.figure() #(figsize=(12, 12))
plot_theta = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$\Theta$', fontsize=16)
# pyplot.yscale('log')

t_max=1e8

root_dir = "output/perpendicular-noweak/"
for filepath in glob.iglob(root_dir + '**/*.txt', recursive=True):
	color = 'k'
	t = []
	theta = []
	t_0 = 0
	with open(filepath) as f:
		for line in f:
			data = line.split()
			if len(data) > 1:
				if isfloat(data[0]) and isfloat(data[1]):
					t_0 = float(data[0])
					if t_0 > t_max:	break
					# print(filepath, line)
					e = float(data[10])
					i = float(data[11])
					t.append(t_0)
					theta.append((1-e**2)*np.cos(i)**2)
				elif data[1] == 'destroyed': color = 'r'
				elif data[1] == 'merger': color = 'g'
	plot_theta.plot(t, theta, color)

pyplot.tight_layout()
pyplot.savefig(root_dir+"theta(t)-zoomin.pdf")