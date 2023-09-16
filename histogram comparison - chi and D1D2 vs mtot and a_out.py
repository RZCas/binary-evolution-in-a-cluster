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

input_files = ['1e5,aout=1','1e5','1e5,aout=3','1e5,aout=4','1e6,aout=1','1e6','1e6,aout=3','1e6,aout=4']
input_files = ['output/chi and D1D2/mtotal='+input_file+'-chiD1D2.txt' for input_file in input_files]
titles = [r'$M_{\rm tot}=10^5M_\odot$, $a_{\rm out}=1$ pc', r'$M_{\rm tot}=10^5M_\odot$, $a_{\rm out}=2$ pc', r'$M_{\rm tot}=10^5M_\odot$, $a_{\rm out}=3$ pc', r'$M_{\rm tot}=10^5M_\odot$, $a_{\rm out}=4$ pc',r'$M_{\rm tot}=10^6M_\odot$, $a_{\rm out}=1$ pc', r'$M_{\rm tot}=10^6M_\odot$, $a_{\rm out}=2$ pc', r'$M_{\rm tot}=10^6M_\odot$, $a_{\rm out}=3$ pc', r'$M_{\rm tot}=10^6M_\odot$, $a_{\rm out}=4$ pc']

figure = pyplot.figure(figsize=(12, 14)) 
figure.suptitle(r'Hernquist, $b=2$ pc, $m_1=m_2=10M_\odot$, $a_{\rm in}=100$ AU, $i=0\dots180^\circ$')
gs = figure.add_gridspec(4, 2)
(ax0, ax4), (ax1, ax5), (ax2, ax6), (ax3, ax7) = gs.subplots(sharex=True,sharey=True)
ax = [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7]

for n in range(8):
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

	chi = np.array(chi)
	D1D2 = np.array(D1D2)
	# figure.suptitle(fr'$\langle\chi\rangle = {chi.mean():.3f}$')
	ax[n].set_title(titles[n]+'\n'+rf'$\chi={chi.mean():.2f}\pm{stdev(chi):.2f}$'+'\n'+rf'$\xi={D1D2.mean():.2f}\pm{stdev(D1D2):.2f}$', x=0.5, y=0.75)
	ax[n].hist (chi, histtype='step', density=True, bins=np.linspace(0, 1, 11), label=r'$\chi$')
	ax[n].hist (D1D2, histtype='step', density=True, bins=np.linspace(0, 1, 11), label=r'$\xi$')
	ax[n].legend(frameon=False)#, loc='lower left')

pyplot.tight_layout()
pyplot.savefig('output/for the paper 2/chi and D1D2 vs mtot and a_out.pdf')
pyplot.savefig(str(Path.home())+'/Dropbox/CLUSTERS/chi and D1D2/chi and D1D2 vs mtot and a_out.pdf')