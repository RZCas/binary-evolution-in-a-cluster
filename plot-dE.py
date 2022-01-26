from galpy.potential import evaluaterforces, evaluatePotentials, PotentialError, KeplerPotential, TwoPowerTriaxialPotential, PlummerPotential
import astropy.units as u
import numpy as np
import os
_pc = 8000
_kms = 220
m_total = 4e+6
b=1
pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc)

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

# def makePlot(fileName1, save_file):
# follow the time evolution of a[AU] ecc inc long_asc arg_peri
t = []
a = []
ecc = []
inc = []
long_asc = []
arg_peri = []
t_previous = 0
dt = []
dE_total = 0
dE_total_dE_0 = []
# filename = os.path.dirname(os.path.abspath(__file__))+"/history-1000n"
filename = "output/t_i-test/t_i-test-noweak"

# R1 = 5.56036948138
# z1 = 1.23961836864 
# phi1 = -2.3542439268
# v1 = np.linalg.norm([-35.6092489981, -7.76059469061, 23.2840513106])
# R2 = 0.924152548492
# z2 = -0.087550416546 
# phi2 = 1.95381235911
# v2 = np.linalg.norm([-23.4603470506, 23.0163896367, 122.284219355])
# E1 = (v1/_kms)**2/2 + evaluatePotentials(pot, R1/_pc, z1/_pc, phi=phi1, use_physical=False) 
# E2 = (v2/_kms)**2/2 + evaluatePotentials(pot, R2/_pc, z2/_pc, phi=phi2, use_physical=False) 
# print(E1, E2)

with open(filename+".txt") as f:
	lineNumber = 0
	for line in f:
		lineNumber+=1
		data = line.split()
		if isfloat(data[0]) and isfloat(data[1]):
			# dt.append(t[-1]-t_previous)
			# t_previous = t[-1]
			R = float(data[1])
			z = float(data[2])
			phi = float(data[3])
			v_R = float(data[4])
			v_z = float(data[5])
			v_phi = float(data[6])
			v = np.sqrt(v_R**2+v_z**2+v_phi**2)
			a.append(float(data[7]))
			ecc.append(float(data[10]))
			inc.append(float(data[11]))
			long_asc.append(float(data[12]))
			arg_peri.append(float(data[13]))
			E = (v/_kms)**2/2 + evaluatePotentials(pot, R/_pc, z/_pc, phi=phi, use_physical=False) 
			if lineNumber == 3:
				E_0 = E
				# print(E_0)
			if lineNumber % 3 == 0: 
				E_prev = E
			if lineNumber % 3 == 1:
				t.append(float(data[0]))
				dE_total += E - E_prev
				dE_total_dE_0.append(dE_total/E_0)
				# if (E - E_prev)/E_0>0.05: 
				# 	print(t[-1]) 
				# 	print(lineNumber)

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
figure = pyplot.figure() #(figsize=(12, 12))
plot_a = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$dE/E_0$', fontsize=16)
# pyplot.xscale('log')
# pyplot.text(1.5, 0.75, '$e='+str(0.999)+'$', fontsize=16)
plot_a.plot(t, dE_total_dE_0, 'k')

pyplot.tight_layout()
pyplot.savefig(filename+"_dE.pdf")