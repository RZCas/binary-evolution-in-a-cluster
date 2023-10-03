t1 = []
de1 = []
r1 = []
t2 = []
de2 = []
r2 = []
fileName1 = '../flybys-master/de(t)1_Q=30_3body_meanAnomaly=0-test.txt'
fileName2 = '../flybys-master/de(t)1_Q=30_3body_meanAnomaly=1.5-test.txt'
save_file = 'output/for the paper/de(t)1_Q=30_3body.pdf'

with open(fileName1) as f:
	for line in f:
		data = line.split()
		t1.append(float(data[0]))
		de1.append(float(data[1])) 
		r1.append(float(data[2])) 
with open(fileName2) as f:
	for line in f:
		data = line.split()
		t2.append(float(data[0]))
		de2.append(float(data[1])) 
		r2.append(float(data[2]))  

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"

figure = pyplot.figure(figsize=(6, 8)) 
gs = figure.add_gridspec(3, 1, height_ratios=[0.7, 0.7, 1])
plot_r, plot_e, plot_e_zoomin = gs.subplots()

# fig, (plot_r, plot_e, plot_e_zoomin) = pyplot.subplots(nrows=3, figsize=(6, 8))

# plot_r = pyplot.subplot(311)
# plot_e = pyplot.subplot(312, sharex=plot_r)
# plot_e_zoomin = pyplot.subplot(313)

x_text = 0.02
y_text = 0.05

plot_r.minorticks_on() 
plot_r.tick_params(labelsize=14)
plot_r.set_xlabel(r'$t$', fontsize=16)
plot_r.set_ylabel(r'$r_3/a$', fontsize=16)
plot_r.plot(t1, r1, 'k')
plot_r.sharex(plot_e)
plot_r.text(x_text, y_text, '(a)', transform=plot_r.transAxes, fontsize=16)
# fig.draw_without_rendering() 
# plot_r.autoscale(enable=False)
xlim = plot_r.get_xlim()
ylim = plot_r.get_ylim()
plot_r.add_patch(pyplot.Rectangle([1200,0], 200, 1000, facecolor='grey', edgecolor='none'))
plot_r.set_xlim(xlim)
plot_r.set_ylim(ylim)

plot_e.minorticks_on() 
plot_e.tick_params(labelsize=14)
plot_e.set_xlabel(r'$t$', fontsize=16)
plot_e.set_ylabel(r'$\Delta e$', fontsize=16)
plot_e.plot(t1, de1, 'k')
plot_e.sharex(plot_r)
plot_e.text(x_text, y_text, '(b)', transform=plot_e.transAxes, fontsize=16)
xlim = plot_e.get_xlim()
ylim = plot_e.get_ylim()
plot_e.add_patch(pyplot.Rectangle([1200,-1], 200, 2, facecolor='grey', edgecolor='none'))
plot_e.set_xlim(xlim)
plot_e.set_ylim(ylim)

plot_e_zoomin.minorticks_on() 
plot_e_zoomin.tick_params(labelsize=14)
plot_e_zoomin.set_xlabel(r'$t$', fontsize=16)
plot_e_zoomin.set_ylabel(r'$\Delta e$', fontsize=16)
plot_e_zoomin.plot(t1, de1, 'k', label=r'initial mean anomaly = 0')
plot_e_zoomin.plot(t2, de2, 'r', label=r'initial mean anomaly = $85.9^\circ$')
plot_e_zoomin.legend(fontsize=16, frameon=False)
plot_e_zoomin.set_xlim(1200, 1400)
plot_e_zoomin.text(x_text, y_text, '(c)', transform=plot_e_zoomin.transAxes, fontsize=16)

pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
pyplot.savefig(save_file)