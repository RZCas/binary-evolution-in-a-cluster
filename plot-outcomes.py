import os
import glob
import numpy as np
import astropy.units as u
from galpy.potential import evaluaterforces, evaluatePotentials, PotentialError, KeplerPotential, TwoPowerTriaxialPotential, PlummerPotential, HernquistPotential
from binary_evolution_with_flybys import a_h, sigma
from amuse.lab import units, constants 
from binary_evolution.tools import rarp
_pc = 8000
_kms = 220
G = constants.G
c = constants.c
t_H = 1.4e10|units.yr
H = 15
potential = "Plummer"

m_per = 1|units.MSun
def tau_0_factor (a, m_bin, r, Q_max_a=50, type="Plummer", m_total=4e6, b=1, V=0|units.kms):
	Q_max = Q_max_a * a
	v0 = np.sqrt(G*(m_bin+m_per)/Q_max)
	sigma_rel = np.sqrt(sigma(r, type, m_total, b)**2 + V**2)
	return (Q_max**2*(1+(v0/sigma_rel)**2))**-1

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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
# pyplot.yscale('log')

t_fin = 1e10
a_unfinished = []
a_survived = []
a_merged = []		#withot exchange
a_destroyed = []	#includes both evaporated and a>1000AU
a_exchange = []

lineNumber = 0
with open('output/outcomes.txt') as f:
	for line in f:
		lineNumber+=1
		data = line.split()
		if lineNumber % 4 == 1:	#initial conditions
			skip = False
			a = float(data[7])
			m_0 = float(data[8])
		if lineNumber % 4 == 2:	#conditions before the finale (to check for exchanges)
			if data[0] == "perturber:":
				a_unfinished.append(a)
				skip = True
			elif float(data[8]) != m_0:
				a_exchange.append(a)	
				skip = True
		if lineNumber % 4 == 3 and not skip:	#outcome
			if not isfloat(data[0]):	#data[0]=="perturber", i.e. stopped during a 3-body interaction
				a_unfinished.append(a)
			else:
				if isfloat(data[1]):
					if float(data[0]) >= t_fin:
						a_survived.append(a)
					else:
						a_unfinished.append(a)
				else:
					if data[1] == 'destroyed': 
						a_destroyed.append(a)
					if data[1] == 'merger': 
						a_merged.append(a)
					if data[1] == 'maximum' and data[2] == 'semimajor': 
						a_destroyed.append(a)

a_min = 2
a_max = 200
n = 10	#number of bins
bins = np.logspace(np.log10(a_min), np.log10(a_max), n+1)
unfinished, bin_edges = np.histogram (a_unfinished, bins=bins)
survived, bin_edges = np.histogram (a_survived, bins=bins)
merged, bin_edges = np.histogram (a_merged, bins=bins)
destroyed, bin_edges = np.histogram (a_destroyed, bins=bins)
exchange, bin_edges = np.histogram (a_exchange, bins=bins)

total = unfinished + survived + merged + destroyed + exchange
unfinished_fraction = unfinished / total
survived_fraction = survived / total
merged_fraction = merged / total
destroyed_fraction = destroyed / total
exchange_fraction = exchange / total

bin_centres = []
for i in range(n):
	bin_centres.append((bin_edges[i]+bin_edges[i+1])/2)

figure = pyplot.figure(figsize=(6, 5))
plot = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$a_0$ [AU]', fontsize=16)
ax.set_ylabel(r'outcome fractions', fontsize=16)
plot.plot(bin_centres, unfinished_fraction, label='unfinished')
plot.plot(bin_centres, survived_fraction, label='survived')
plot.plot(bin_centres, merged_fraction, label='merged')
plot.plot(bin_centres, destroyed_fraction, label='destroyed')
plot.plot(bin_centres, exchange_fraction, label='exchange')
pyplot.xscale("log")
plot.legend()
pyplot.tight_layout()
pyplot.savefig("output/outcomes.pdf")
pyplot.show()