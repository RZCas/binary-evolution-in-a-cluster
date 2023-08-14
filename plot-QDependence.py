def plot_data(fileName1):
	x1 = []
	y1 = []
	x2 = []
	y2 = []
	x1neg = []
	y1neg = []
	x2neg = []
	y2neg = []
	x_sa = []
	y_sa = []
	with open(fileName1) as f:
		lineNumber = 0
		for line in f:
			lineNumber+=1
			if lineNumber>1:
				data = line.split()
				value1 = float(data[1])
				if value1>0: 
					x1.append(float(data[0]))
					y1.append(value1) 
				else:
					x1neg.append(float(data[0]))
					y1neg.append(-value1)
				value2 = float(data[2])
				if value2>0: 
					x2.append(float(data[0]))
					y2.append(value2) 
				else:
					x2neg.append(float(data[0]))
					y2neg.append(-value2)	
	return x1, y1, x1neg, y1neg, x2, y2, x2neg, y2neg

x1_1, y1_1, x1neg_1, y1neg_1, x2_1, y2_1, x2neg_1, y2neg_1 = plot_data('../flybys-master/QDependence4_eStar=1.01068346206_r3max=5000.txt')
x1_2, y1_2, x1neg_2, y1neg_2, x2_2, y2_2, x2neg_2, y2neg_2 = plot_data('../flybys-master/QDependence5_r3max=5000.txt')

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"

figure = pyplot.figure(figsize=(6, 7)) 
gs = figure.add_gridspec(2, 1, hspace=0)
plot1, plot2 = gs.subplots(sharex=True)

plot1.minorticks_on() 
plot1.tick_params(labelsize=14)
# plot1.set_xlabel(r'$Q/a$', fontsize=16)
plot1.set_ylabel(r'$|\Delta e|$', fontsize=16)
plot1.set_yscale('log')
plot1.plot(x1_1, y1_1, 'k', label='3-body')
plot1.plot(x1neg_1, y1neg_1, 'k--')
plot1.plot(x2_1, y2_1, 'r', label='orbit-averaged')
plot1.plot(x2neg_1, y2neg_1, 'r--')
plot1.legend(fontsize=16, frameon=False)

plot2.minorticks_on() 
plot2.tick_params(labelsize=14)
plot2.set_xlabel(r'$Q/a$', fontsize=16)
plot2.set_ylabel(r'$|\Delta e|$', fontsize=16)
plot2.set_yscale('log')
plot2.plot(x1_2, y1_2, 'k', label='3-body')
plot2.plot(x1neg_2, y1neg_2, 'k--')
plot2.plot(x2_2, y2_2, 'r', label='orbit-averaged')
plot2.plot(x2neg_2, y2neg_2, 'r--')
plot2.legend(fontsize=16, frameon=False)

pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
pyplot.savefig("output/for the paper/QDependence.pdf")
pyplot.clf()