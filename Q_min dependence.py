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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{color}"

t_max=1e80
a_out = 3
m_total = 1e6
b = 2
A_ast = 0.3 #A_* for Hernquist potential
potential = "Hernquist"
pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc)
# potential = "Plummer"
# pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

Q_max_a = 10	#maximum value of Q_min we consider
n = 50			#number of Q_min values
Q_a_n = [i/n*Q_max_a for i in range(n)]

types = ['hernquist,m_total=1e5,b=1,a_out=4,i=89.9,nokicks,a_in=300', 'hernquist,m_total=1e5,b=1,a_out=4,i=89.9,nokicks,a_in=300,ns','hernquist,m_total=1e5,b=1,a_out=4,i=89.9,nokicks','hernquist,m_total=1e5,b=1,a_out=3,i=89.9,nokicks','perpendicular-soft-hernquist-nokicks','perpendicular-soft-hernquist-nokicks-weakonly']
for type in types[1:]:
	root_dir = "output/"+type+"/"
	for index in range(20):
		shift = 0
		filepath = root_dir + str(index) + '.txt'
		e = []
		de_flybys_max = [0]*n 	#max value of |de| of all time, corresponding to the n values of Q_min
		de_flybys = [0]*n 	#de at the current time, corresponding to the n values of Q_min
		de_tidal = 0
		de_tidal_max = 0
		lineNumber = 0
		with open(filepath) as f:
			for line in f:
				lineNumber+=1
				data = line.split()
				if len(data) > 1:
					if data[0]=="perturber:":
						if lineNumber>2: 
							Q_a = float(data[2]) / a
						if lineNumber%3==0:
							shift=1
					if isfloat(data[0]) and isfloat(data[1]):
						a = float(data[7])
						e.append(float(data[10]))
						if lineNumber%3==1+shift and lineNumber>1:
							de_tidal += e[-1]-e[-2]
							if abs(de_tidal) > de_tidal_max:
								de_tidal_max = abs(de_tidal)
						if lineNumber%3==0+shift and lineNumber>3:
							n_current = int(Q_a / Q_max_a * n)
							if n_current >= n:  
								n_current = n-1
							for j in range(0,n_current+1):
								de_flybys[j] += e[-1]-e[-2]
								if abs(de_flybys[j]) > de_flybys_max[j]:
									de_flybys_max[j] = abs(de_flybys[j])
		figure = pyplot.figure(figsize=(9, 7))
		plot = figure.add_subplot(1,1,1)
		ax = pyplot.gca()
		ax.minorticks_on() 
		ax.tick_params(labelsize=14)
		ax.set_xlabel(r'$Q_{\rm min}/a$', fontsize=16)
		ax.set_ylabel(r'$|\Delta e|_{\rm max}$', fontsize=16)
		plot.plot(Q_a_n, de_flybys_max, 'k', label='encounters')
		plot.plot([0,Q_max_a], [de_tidal_max,de_tidal_max], 'r', label='tidal effects')
		plot.legend()
		pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
		pyplot.savefig(root_dir+"de_max-"+type+"-"+str(index)+".pdf")
		pyplot.close()