import os
import glob
from amuse.lab import units, constants 
import numpy as np

G = constants.G
c = constants.c
t_H = 1.4e10|units.yr
H = 15
# assuming m1=m2, m=m1+m2
m_i = 1
m_f = 100

def a_h (m, m_cl, b): 
	return b|units.pc*(m|units.MSun)/(m_cl|units.MSun)
def a_3body (m, m_cl, b):
	return 4*np.pi / (3*H*t_H) * sqrt((b|units.pc)**5/G/(m_cl|units.MSun))
def a_gw (m, m_cl, b):
	return (16*t_H*G**3*(m|units.MSun)**3/c**5)**0.25
def a_tidal (m, m_cl, b):
	A = 0.5*G*(m_cl|units.MSun)/(b|units.pc)**3
	return (24*G**2*(m|units.MSun)**2/c**2/A)**0.25
def a_tsec01tH (m, m_cl, b):
	A = 0.5*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return ((8/(3*A*t1))**2*G*m|units.MSun)**(1/3)

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
figure = pyplot.figure() #(figsize=(12, 12))
b_values = (1, 3)
m_cl_values = (1e4, 1e5, 1e6)

index=0
for n_b in range(len(b_values)):
	for n_m_cl in range(len(m_cl_values)):	
		m_cl = m_cl_values[n_m_cl]
		b = b_values[n_b]
		a_h_i = a_h(m_i, m_cl, b).value_in(units.pc)
		a_h_f = a_h(m_f, m_cl, b).value_in(units.pc)
		index += 1
		plot = figure.add_subplot(len(m_cl_values),len(b_values),index)
		ax = pyplot.gca()
		ax.minorticks_on() 
		# ax.xaxis.set_major_locator(MultipleLocator(1))
		# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
		# ax.yaxis.set_major_locator(MultipleLocator(0.1))
		# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
		ax.tick_params(labelsize=14)
		ax.set_xlabel(r'$m_{\rm bin}$ [MSun]', fontsize=16)
		ax.set_ylabel(r'$a$ [pc]', fontsize=16)
		# pyplot.yscale('log')
		# pyplot.text(1.5, 0.75, '$e='+str(0.999)+'$', fontsize=16)
		plot.plot ((m_i, m_f), (a_h_i, a_h_f), 'k')

pyplot.tight_layout()
pyplot.savefig("a_h(m).pdf")