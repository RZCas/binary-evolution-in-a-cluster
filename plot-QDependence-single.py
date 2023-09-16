import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"

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

filename = ['QDependence4-hamers-test','QDependence4,e=0.9-test','QDependence4-hamers-fig4a','QDependence4-test']
label = [r'$e=0.1$, $E=1.5$', r'$e=0.9$, $v=3$ km/s', r'$e=0.999$, $E=1.5$ (Hamers \& Samsing 2019, Fig. 4a)', r'$e=0.1$, $v=3$ km/s']

for i in range(len(filename)):
	x1_1, y1_1, x1neg_1, y1neg_1, x2_1, y2_1, x2neg_1, y2neg_1 = plot_data(f'../flybys-master/{filename[i]}.txt')

	figure = pyplot.figure(figsize=(6, 5)) 
	gs = figure.add_gridspec(1, 1, hspace=0)
	plot1 = gs.subplots(sharex=True)

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
	plot1.set_title(label[i], fontsize=16)

	pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
	pyplot.savefig(f"output/for the paper/{filename[i]}.pdf")
	pyplot.clf()