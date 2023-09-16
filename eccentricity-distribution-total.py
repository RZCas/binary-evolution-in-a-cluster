import numpy as np
import glob, os, subprocess
import matplotlib
from pathlib import Path
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{color}"

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

root_folder = 'output/m1=m2=10/'
simulation_sets = ['mtotal=1e5/','mtotal=1e5,enc_only/']
for simulation_set in simulation_sets:
	e = []
	folder = root_folder + simulation_set[:-1] + '-e-distribution/'
	for file in glob.glob (folder + "*.txt"):
		with open(file) as f:
			for line in f:
				e.append(float(line))

	fig = pyplot.figure(figsize=(6, 4)) 
	pyplot.hist (e, histtype='step', density=True, bins=np.linspace(0, 1, 31))
	pyplot.xlabel (r'$e$', fontsize=16)
	pyplot.tight_layout()
	pyplot.savefig(folder + simulation_set[:-1] + '-e-distribution.pdf')
	pyplot.close()