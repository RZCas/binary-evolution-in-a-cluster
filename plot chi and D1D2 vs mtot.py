import os
import glob
import shutil
import numpy as np
from pathlib import Path
from statistics import stdev

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
matplotlib.rcParams["errorbar.capsize"] = 10
# matplotlib.rcParams["figure.titlesize"] = 4

input_files = ['3e4','1e5','2e5','3e5','6e5','1e6','3e6','1e7']
mtot = [float(input_file) for input_file in input_files]
chi_mean = []
d1d2_mean = []
chi_sigma = []
d1d2_sigma = []
input_files = ['output/chi and D1D2/mtotal='+input_file+'-chiD1D2.txt' for input_file in input_files]

for input_file in input_files:
	chi = []
	D1D2 = []
	with open(input_file) as f:
		for line in f:
			data = line.split()
			if isfloat(data[1]):
				chi.append(float(data[1]))
				if len(data)>2:
					D1D2.append(float(data[2]))
	chi = np.array(chi)
	D1D2 = np.array(D1D2)
	chi_mean.append(chi.mean())
	d1d2_mean.append(D1D2.mean())
	chi_sigma.append(stdev(chi))
	d1d2_sigma.append(stdev(D1D2))

figure = pyplot.figure(figsize=(6, 5)) 
figure.suptitle(r'Hernquist, $b=a_{\rm out}=2$ pc, $m_1=m_2=10M_\odot$, $a_{\rm in}=100$ AU, $i=0\dots180^\circ$')
pyplot.xlabel(r'$M_{\rm tot}$ $[M_\odot]$', fontsize=16)
pyplot.xscale('log')
pyplot.ylim(0,1)
pyplot.errorbar(mtot, chi_mean, chi_sigma, label=r'$\chi$')
pyplot.errorbar(mtot, d1d2_mean, d1d2_sigma, label=r'$D_1/D_2$')
pyplot.legend(frameon=False)
pyplot.tight_layout()
pyplot.savefig('output/chi and D1D2/plot - chi and D1D2 vs mtot.pdf')
pyplot.savefig(str(Path.home())+'/Dropbox/CLUSTERS/chi and D1D2/plot - chi and D1D2 vs mtot.pdf')