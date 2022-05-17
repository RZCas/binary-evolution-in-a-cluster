import os
import glob
import numpy as np
import astropy.units as u
from galpy.potential import evaluaterforces, evaluatePotentials, PotentialError, KeplerPotential, TwoPowerTriaxialPotential, PlummerPotential, HernquistPotential
from binary_evolution_with_flybys import a_h
from amuse.lab import units, constants 
_pc = 8000
_kms = 220
G = constants.G
c = constants.c
t_H = 1.4e10|units.yr
H = 15
potential = "Hernquist"

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

def f(x): return 12*x*(x+1)**3*np.log(1+1/x) - x/(x+1)*(25+52*x+42*x**2+12*x**3)

x_values = np.linspace (2000, 10000, 1000)
f_values = f(x_values)

from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
figure = pyplot.figure(figsize=(6, 4))
plot_theta = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$r/a$', fontsize=16)
ax.set_ylabel(r'$\sigma / (GM/12a)$', fontsize=16)
plot_theta.plot(x_values, f_values, 'k')
pyplot.tight_layout()
pyplot.savefig("output/Hernquist-sigma-3.pdf")