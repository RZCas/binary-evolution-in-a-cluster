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
# matplotlib.rcParams["figure.titlesize"] = 4

input_files = ['1e5','3e5','1e6','3e6','1e7']
input_files = ['output/chi and D1D2/mtotal='+input_file+'-chiD1D2.txt' for input_file in input_files]
titles = [r'$M_{\rm tot}=10^5M_\odot$', r'$M_{\rm tot}=3\cdot10^5M_\odot$', r'$M_{\rm tot}=10^6M_\odot$', r'$M_{\rm tot}=3\cdot10^6M_\odot$', r'$M_{\rm tot}=10^7M_\odot$']

figure = pyplot.figure(figsize=(6, 15)) 
figure.suptitle(r'Hernquist, $b=a_{\rm out}=2$ pc, $m_1=m_2=10M_\odot$, $a_{\rm in}=100$ AU, $i=0\dots180^\circ$')
gs = figure.add_gridspec(5, 1)
(ax0, ax1, ax2, ax3, ax4) = gs.subplots(sharex=True,sharey=True)
ax = [ax0, ax1, ax2, ax3, ax4]

for n in range(len(ax)):
	input_file = input_files[n]
	chi = []
	D1D2 = []

	with open(input_file) as f:
		for line in f:
			data = line.split()
			if isfloat(data[1]):
				chi.append(float(data[1]))
				if len(data)>2:
					D1D2.append(float(data[2]))
					# if n==0 and D1D2[-1]<0.02:
					# 	print("simulation "+data[0]+" has D1/D2 = "+data[2])

	chi = np.array(chi)
	D1D2 = np.array(D1D2)
	# figure.suptitle(fr'$\langle\chi\rangle = {chi.mean():.3f}$')
	ax[n].set_title(titles[n]+'\n'+rf'$\chi={chi.mean():.2f}\pm{stdev(chi):.2f}$'+'\n'+rf'$D_1/D_2={D1D2.mean():.2f}\pm{stdev(D1D2):.2f}$', x=0.5, y=0.75)
	ax[n].hist (chi, histtype='step', density=True, bins=np.linspace(0, 1, 11), label=r'$\chi$')
	ax[n].hist (D1D2, histtype='step', density=True, bins=np.linspace(0, 1, 11), label=r'$D_1/D_2$')
	ax[n].legend(frameon=False)#, loc='lower left')

pyplot.tight_layout()
pyplot.savefig('output/chi and D1D2/chi and D1D2 vs mtot.pdf')
pyplot.savefig(str(Path.home())+'/Dropbox/CLUSTERS/chi and D1D2/chi and D1D2 vs mtot.pdf')