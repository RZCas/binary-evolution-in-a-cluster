import os
import glob
import numpy as np
from binary_evolution_with_flybys import a_h
from amuse.lab import units, constants 

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

t_max=1e80
a_out = 2
m_total = 1e6
b = 2
A_ast = 0.3 #A_* for Hernquist potential
# print(a_h(10, 10, a_out, type="Hernquist", m_total=1e6, b=1))

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

# types = ['wide_range_1','wide_range_mtotal=1e5','wide_range_m=60','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e6_nokicks_plummer','wide_range_mtotal=1e6_nokicks']
types = ['wide_range_mtotal=1e6_nokicks']
for type in types:
	root_dir = "output/"+type+"/"
	for index in range(0,10):
		notYet = True
		filepath = root_dir + str(index) + '.txt'
		color = 'k'
		t = []
		theta = []
		e = []
		cosi = []
		a = []
		t_0 = 0

		n = 30	#number of Q/a bins
		Qa_max = 50
		de = [0]*n
		binEdges = [Qa_max * i/n for i in range(1,n+1)] 
		lineNumber=0
		with open(filepath) as f:
			for line in f:
				lineNumber+=1
				data = line.split()
				if len(data) > 1:
					if lineNumber == 3:	#initial parameters
						a_i = float(data[7])
						m = float(data[8])
						q = float(data[9])
						m1 = m/(1+q)
						m2 = m*q/(1+q)
					if lineNumber%3==1 and lineNumber>=4 and isfloat(data[1]):	#before the encounter
						a = float(data[7])
						e_i = float(data[10])
					if data[0]=="perturber:" and lineNumber>=5:	#encounter parameters
						Q = float(data[2])
					if lineNumber%3==0 and lineNumber>=6 and isfloat(data[1]):	#after the encounter
						e_f = float(data[10])
						binNumber = int(Q/a / Qa_max * n) 
						for i in range(binNumber, n):
							de[i] += e_f-e_i
		de.insert(0,0)
		binEdges.insert(0,0)
		figure = pyplot.figure(figsize=(9/2, 7/2))
		figure.suptitle(r'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_0$ = {a_0:.1f} AU'.format(m1=m1, m2=m2, a_0=a_i, fontsize=12))
		plot_de = figure.add_subplot(1,1,1)
		pyplot.xlabel(r'$Q_{\rm max}/a$')
		pyplot.ylabel(r'$(\Delta e)_{\rm total}(Q<Q_{\rm max})$')
		plot_de.plot (binEdges, de, 'k')
		pyplot.tight_layout()
		pyplot.savefig(root_dir+"de_total(Qa_max)-"+type+"-"+str(index)+".pdf")
		pyplot.clf()