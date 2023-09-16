x1 = []
y1 = []
x2 = []
y2 = []
with open('../flybys-master/r3maxDependence3_sa.txt') as f:
	# lineNumber = 0
	for line in f:
		# lineNumber+=1
		# if lineNumber>1:
		data = line.split()
		x1.append(float(data[0]))
		y1.append(float(data[1]))
with open('../flybys-master/r3maxDependence4_sa.txt') as f:
	# lineNumber = 0
	for line in f:
		# lineNumber+=1
		# if lineNumber>1:
		data = line.split()
		x2.append(float(data[0]))
		y2.append(float(data[1]))

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
figure = pyplot.figure(figsize=(6, 4))
plot = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$r_{\rm 3,max}/a$', fontsize=16)
ax.set_ylabel(r'$\Delta e$', fontsize=16)
pyplot.xscale('log')
# pyplot.text(1.5, 0.75, '$e='+str(0.999)+'$', fontsize=16)
plot.plot(x1, y1, 'k')
plot.plot(x2, y2, 'r')
# plot.plot(x2, y2, 'r', label=r'$\Delta t=1$')
# plot.scatter(x1, y2, label=r'$\delta e$')
# plot.legend(fontsize=16)
# plot.set_xlim(5, 20)
# plot.set_ylim(-1e-5, 1e-5)
pyplot.tight_layout()
pyplot.savefig('output/for the paper/r3maxDependence.pdf')