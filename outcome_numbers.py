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
root_dir = "output/m1=m2=10/"
for filepath in glob.iglob(root_dir + 'outcomes-*.txt'):
	output_file = open(root_dir + 'outcome_numbers' + os.path.basename(filepath)[8:], 'w+')
	lineNumber = 0
	n_total = 0
	n_unfinished = 0
	n_survived = 0
	n_merged = 0	#withot exchange
	n_destroyed = 0	#includes both evaporated and a>1000AU
	n_exchange = 0
	n_ejected = 0
	with open(filepath) as f:
		for line in f:
			lineNumber+=1
			data = line.split()
			if lineNumber % 5 == 2:	#initial conditions
				n_total += 1
				skip = False
				a = float(data[7])
				m_0 = float(data[8])
			if lineNumber % 5 == 3:	#conditions before the finale (to check for exchanges)
				if data[0] == "perturber:":
					n_unfinished += 1
					skip = True
				elif isfloat(data[8]) and float(data[8]) != m_0:
					n_exchange += 1	
					skip = True
			if lineNumber % 5 == 4 and not skip:	#outcome
				if not isfloat(data[0]):	#data[0]=="perturber", i.e. stopped during a 3-body interaction
					n_unfinished += 1
				else:
					if isfloat(data[1]):
						if float(data[0]) >= t_fin:
							n_survived += 1
						else:
							n_unfinished += 1
					else:
						if data[1] == 'destroyed': 
							n_destroyed += 1
						if data[1] == 'merger': 
							n_merged += 1
						if data[1] == 'maximum' and data[2] == 'semimajor': 
							n_destroyed += 1
						if data[1] == 'ejected': 
							n_ejected += 1
	print(f'{n_total} simulations:', file=output_file)
	print(f'unfinished: {int(100*n_unfinished/n_total)}%', file=output_file)
	print(f'survived: {int(100*n_survived/n_total)}%', file=output_file)
	print(f'merged: {int(100*n_merged/n_total)}%', file=output_file)
	print(f'destroyed: {int(100*n_destroyed/n_total)}%', file=output_file)
	print(f'exchange: {int(100*n_exchange/n_total)}%', file=output_file)
	print(f'ejected: {int(100*n_ejected/n_total)}%', file=output_file)
	print(f'\n{n_total-n_unfinished} finished simulations:', file=output_file)
	print(f'survived: {int(100*n_survived/(n_total-n_unfinished))}%', file=output_file)
	print(f'merged: {int(100*n_merged/(n_total-n_unfinished))}%', file=output_file)
	print(f'destroyed: {int(100*n_destroyed/(n_total-n_unfinished))}%', file=output_file)
	print(f'exchange: {int(100*n_exchange/(n_total-n_unfinished))}%', file=output_file)
	print(f'ejected: {int(100*n_ejected/(n_total-n_unfinished))}%', file=output_file)
	output_file.close()
