x1 = []
y1 = []
x2 = []
y2 = []
fileName1 = '../flybys-master/de(t)1_Q=30_3body_meanAnomaly=0-test.txt'
fileName2 = '../flybys-master/de(t)1_Q=30_3body_meanAnomaly=1.5-test.txt'
save_file = 'output/for the paper/'
with open(fileName1) as f:
	for line in f:
		data = line.split()
		x1.append(float(data[0]))
		y1.append(float(data[1])) 
with open(fileName2) as f:
	for line in f:
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
figure = pyplot.figure() #(figsize=(12, 12))
plot = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.tick_params(labelsize=14)
# ax.set_xlabel(r'$r_3$', fontsize=16)
# ax.set_ylabel(r'$\delta e$', fontsize=16)
ax.set_xlabel(r'$t$', fontsize=16)
ax.set_ylabel(r'$\Delta e$', fontsize=16)
# pyplot.xscale('log')
# pyplot.yscale('log')
# pyplot.text(1.5, 0.75, '$e='+str(0.999)+'$', fontsize=16)
# plot.plot(x40, y40, 'k', label=r'3-body, $r_{3,\mathrm{max}}=40$')
plot.plot(x1, y1, 'k', label=r'initial mean anomaly = 0')
plot.plot(x2, y2, 'r', label=r'initial mean anomaly = 1.5')
# plot.scatter(time1, de1, label=r'$dt=1$'))
# plot.plot(x2, y2, 'k--', label=r'3-body, $\delta t=100$')
# plot.plot(x_sa, y_sa, 'r', label=r'SA')
plot.legend(fontsize=16, frameon=False)
plot.set_xlim(1200, 1400)
# plot.set_ylim(-5e-5, 5e-5)
pyplot.tight_layout()
pyplot.savefig(save_file)