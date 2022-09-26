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

# types = ['perpendicular-soft-hernquist', 'perpendicular-hard-hernquist']
# types = ['perpendicular-hard-plummer', 'perpendicular-hard-plummer-light']
# types = ['wide_range_1','wide_range_mtotal=1e5','wide_range_m=60','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e6_nokicks_plummer','wide_range_mtotal=1e6_nokicks','wide_range_mtotal=1e6']
# types = ['wide_range_1','wide_range_mtotal=1e6_nokicks']
# types = ['wide_range_mtotal=1e5','wide_range_mtotal=1e5_nokicks']
# types = ['wide_range_mtotal=1e6_nokicks_plummer']
# types = ['wide_range_mtotal=1e6']
types = os.listdir('storage')
for type in types:
	if os.path.isfile('storage/'+type):
		continue
	root_dir = 'storage/'+type+'/'
	outfile = open (root_dir + type + '-tidal_effects_strength.txt', 'w+')
	print('index metric1 metric2 metric3', file=outfile)
	print('metric1 = sum(|de_t|)/(sum(|de_t|)+sum(|de_k|))', file=outfile)
	print('metric2 = sum(|de_t|/(|de_t|+|de_k|)*dt/T)', file=outfile)
	print('metric3 = sum(|de_t|/(|de_t|+|de_k|))', file=outfile)
	# for filepath in glob.glob(root_dir+'*.txt'):
	for index in range(0,1000):
		filepath = root_dir + str(index) + '.txt'
		if not os.path.isfile(filepath):
			continue
		print(type, index)
		lineNumber = 0
		shift = 0
		metric2 = 0
		metric3 = 0
		sum_de_t = 0
		sum_de_tk = 0
		epsilon = 0
		total_encounters = 0
		total_time = 0
		t = []
		e = []
		de_ratios = []
		de_k_array = []
		de_t_array = []
		with open(filepath) as f:
			for line in f:
				lineNumber+=1
				data = line.split()
				if len(data) > 1:
					if data[0]=="perturber:" and lineNumber%3==0:
						shift=1
					if isfloat(data[0]) and isfloat(data[1]):
						t.append(float(data[0])/1e9)
						if len(data)<11:
							break
						e.append(float(data[10]))
						if lineNumber%3==1+shift and lineNumber>1:	#before encounter
							if len(data)<18:
								break
							epsilon = float(data[17])
							if epsilon < 20:
								de_t = abs(e[-1]-e[-2])
								if de_t>0:
									de_t_array.append(de_t)
								sum_de_t += de_t
								sum_de_tk += de_t
								total_time += t[-1]-t[-2]
						if lineNumber%3==0+shift and lineNumber>3 and epsilon<20:	#after encounter
							de_k = abs(e[-1]-e[-2])
							if de_k>0:
								de_k_array.append(de_k)
							sum_de_tk += de_k
							metric2 += de_t / (de_t+de_k) * (t[-2]-t[-3])
							metric3 += de_t / (de_t+de_k)
							if de_t>0:
								de_ratios.append(de_t / (de_t+de_k))
							total_encounters += 1
		# index = os.path.split(filepath)[1][:-4]
		if total_encounters > 100:
			metric1 = sum_de_t / sum_de_tk
			metric2 /= total_time
			metric3 /= total_encounters
			print(index, metric1, metric2, metric3, file=outfile)
			# de_ratios = np.array(de_ratios)
			# de_t_array = np.array(de_t_array)
			# de_k_array = np.array(de_k_array)
			# figure = pyplot.figure(figsize=(9/2, 7/2))
			# figure.suptitle(fr'$\xi_1 = {metric1:.3f}$, $\xi_2 = {metric3:.3f}$')
			# plot_de = figure.add_subplot(1,1,1)
			# plot_de.hist (de_ratios, histtype='step', bins=np.geomspace(de_ratios.min(), de_ratios.max(), 100), log=True, label=r'$|\Delta e_{\rm t,i}|/(|\Delta e_{\rm t,i}|+|\Delta e_{\rm k,i})|)$')
			# plot_de.hist (de_t_array, histtype='step', bins=np.geomspace(de_t_array.min(), de_t_array.max(), 100), label=r'$|\Delta e_{\rm t,i}|$', log=True)
			# plot_de.hist (de_k_array, histtype='step', bins=np.geomspace(de_k_array.min(), de_k_array.max(), 100), label=r'$|\Delta e_{\rm k,i}|$', log=True)
			# plot_de.legend(loc='lower left')
			# pyplot.xscale ('log')
			# pyplot.tight_layout()
			# pyplot.savefig(root_dir+"de-distribution-"+type+'-'+str(index)+".pdf")
			# pyplot.clf()
		else:
			print(index, 'epsilon>20 at (almost) all times', file=outfile)