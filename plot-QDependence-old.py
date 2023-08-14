def makePlot(fileName1, save_file):
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
	# with open(fileName2) as f:
	# 	# lineNumber = 0
	# 	for line in f:
	# 		# lineNumber+=1
	# 		# if lineNumber>1:
	# 		data = line.split()
	# 		x2.append(float(data[0]))
	# 		y2.append(float(data[1]))
	# with open(fileName_sa) as f:
	# 	# lineNumber = 0
	# 	for line in f:
	# 		# lineNumber+=1
	# 		# if lineNumber>1:
	# 		data = line.split()
	# 		x_sa.append(float(data[0]))
	# 		y_sa.append(float(data[1]))		
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
	ax.set_xlabel(r'$Q/a$', fontsize=16)
	ax.set_ylabel(r'$\Delta e$', fontsize=16)
	# pyplot.xscale('log')
	pyplot.yscale('log')
	# pyplot.text(1.5, 0.75, '$e='+str(0.999)+'$', fontsize=16)
	# plot.plot(x40, y40, 'k', label=r'3-body, $r_{3,\mathrm{max}}=40$')
	plot.plot(x1, y1, 'k', label='3-body')
	plot.plot(x1neg, y1neg, 'k--')
	plot.plot(x2, y2, 'r', label='orbit-averaged')
	plot.plot(x2neg, y2neg, 'r--')
	# plot.plot(x2, y2, 'k--', label=r'3-body, $\delta t=100$')
	# plot.plot(x_sa, y_sa, 'r', label=r'SA')
	plot.legend(fontsize=16)
	# plot.set_xlim(500, 5000)
	# plot.set_ylim(1e-8, 1e-5)
	pyplot.tight_layout()
	pyplot.savefig(save_file)

makePlot('QDependence4_eStar=1.01068346206_r3max=5000.txt', 'QDependence4.pdf')
makePlot('QDependence5_r3max=5000.txt', 'QDependence5.pdf')