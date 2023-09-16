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

e_tidal = []
e_notidal = []

for file in glob.glob ('output/m1=m2=10/mtotal=1e5-e-distribution/*-e-distribution.txt'):
	with open(file) as f:
		for line in f:
			e_tidal.append(float(line))
for file in glob.glob ('output/m1=m2=10/mtotal=1e5,enc_only-e-distribution/*-e-distribution.txt'):
	with open(file) as f:
		for line in f:
			e_notidal.append(float(line))

fig = pyplot.figure(figsize=(6, 4)) 
pyplot.hist (e_tidal, histtype='step', density=True, bins=np.linspace(0, 1, 31), label='flybys+tides')
pyplot.hist (e_notidal, histtype='step', density=True, bins=np.linspace(0, 1, 31), label='flybys')
pyplot.xlabel (r'$e$', fontsize=16)
pyplot.legend()
pyplot.tight_layout()
pyplot.savefig('output/m1=m2=10/mtotal=1e5-e-distribution-comparison.pdf')
