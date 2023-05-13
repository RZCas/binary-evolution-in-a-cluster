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
	return 4*np.pi / (3*H*t_H) * np.sqrt((b|units.pc)**5/G/(m_cl|units.MSun))
def a_gw (m, m_cl, b):
	return (16*t_H*G**3*(m|units.MSun)**3/c**5)**0.25
def a_tidal (m, m_cl, b):
	epsilon = 20
	A = 0.5*G*(m_cl|units.MSun)/(b|units.pc)**3
	return (24*G**2*(m|units.MSun)**2/c**2/A/epsilon)**0.25
def a_tsec01tH (m, m_cl, b):
	A = 0.5*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return ((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
b_values = (1, 3)
m_cl_values = (1e4, 1e5, 1e6)

index=0
for n_b in range(len(b_values)):
	for n_m_cl in range(len(m_cl_values)):	
		index+=1
		m_cl = m_cl_values[n_m_cl]
		b = b_values[n_b]
		a_h_i = a_h(m_i, m_cl, b).value_in(units.AU)
		a_h_f = a_h(m_f, m_cl, b).value_in(units.AU)
		a_3body_i = a_3body(m_i, m_cl, b).value_in(units.AU)
		a_3body_f = a_3body(m_f, m_cl, b).value_in(units.AU)
		a_gw_i = a_gw(m_i, m_cl, b).value_in(units.AU)
		a_gw_f = a_gw(m_f, m_cl, b).value_in(units.AU)
		a_tidal_i = a_tidal(m_i, m_cl, b).value_in(units.AU)
		a_tidal_f = a_tidal(m_f, m_cl, b).value_in(units.AU)
		a_tsec01tH_i = a_tsec01tH(m_i, m_cl, b).value_in(units.AU)
		a_tsec01tH_f = a_tsec01tH(m_f, m_cl, b).value_in(units.AU)		
		figure = pyplot.figure() #(figsize=(12, 12))
		# figure.suptitle(r'$M_{\rm cl}=10^'+str(int(np.log10(m_cl)))+r'$ $M_\odot$, $b=$ '+str(b)+r' pc', fontsize=16)
		plot = figure.add_subplot(1,1,1)
		ax = pyplot.gca()
		ax.minorticks_on() 
		ax.tick_params(labelsize=14)
		ax.set_xlabel(r'$m_{\rm bin}$ [$M_\odot$]', fontsize=16)
		ax.set_ylabel(r'$a_{\rm in}$ [AU]', fontsize=16)
		pyplot.xscale('log')
		pyplot.yscale('log')
		# pyplot.text(1, 1, r'$M_{\rm cl}=10^'+str(int(np.log10(m_cl)))+r'$ $M_\odot$, $b=$ '+str(b)+r' pc', fontsize=16)
		ax.text(x=0.02, y=0.93, s=r'$M_{\rm cl}=10^'+str(int(np.log10(m_cl)))+r'$ $M_\odot$, $b=$ '+str(b)+r' pc', fontsize=16, transform=ax.transAxes)
		plot.plot ((m_i, m_f), (a_h_i, a_h_f), 'k', label=r'$a_h$')
		plot.plot ((m_i, m_f), (a_3body_i, a_3body_f), 'k--', label=r'$t_{\rm 3body}=t_H$')
		plot.plot ((m_i, m_f), (a_gw_i, a_gw_f), 'r', label=r'$t_{\rm GW}=t_H$')
		plot.plot ((m_i, m_f), (a_tidal_i, a_tidal_f), 'g', label=r'$\epsilon_{\rm GR}=1$')
		plot.plot ((m_i, m_f), (a_tsec01tH_i, a_tsec01tH_f), 'g--', label=r'$t_{\rm sec}=0.1t_H$')	
		plot.legend(loc='lower right', fontsize=16)	
		pyplot.tight_layout()
		pyplot.savefig("output/a_h(m)/a_h(m)_"+str(index)+".pdf")