# plot the histogram of |de_tidal|/|de_encounters| for a certain cluster parameters
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
# matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
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

types = ['wide_range_1','wide_range_mtotal=1e5','wide_range_m=60','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e6_nokicks_plummer','wide_range_mtotal=1e6_nokicks','wide_range_mtotal=1e5_plummer','wide_range_mtotal=1e6_plummer','wide_range_mtotal=1e5_plummer_b=1','wide_range_mtotal=1e7_hernquist']
# types = ['wide_range_1']
for type in types:
	de_ratios = []
	root_dir = "output/"+type+"/"
	index = 0
	while True:	#for every file in the folder
		filepath = root_dir + str(index) + '.txt'
		index += 1
		if not os.path.isfile(filepath): 
			break
		lineNumber = 0
		de_encounters_total = 0
		de_tidal_total = 0
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
						de_tidal_total += float(data[19])
					if data[0]=="perturber:" and lineNumber>=5:	#encounter parameters
						Q = float(data[2])
					if lineNumber%3==0 and lineNumber>=6 and isfloat(data[1]):	#after the encounter
						e_f = float(data[10])
						de_encounters_total += abs(e_f-e_i)
		de_ratios.append(de_tidal_total/de_encounters_total)
	
	de_ratios = np.array(de_ratios)
	de_ratios = np.log10(de_ratios)

	figure = pyplot.figure(figsize=(9/2, 7/2))
	# figure.suptitle(r'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_0$ = {a_0:.1f} AU'.format(m1=m1, m2=m2, a_0=a_i, fontsize=12))

	plot_de = figure.add_subplot(1,1,1)
	pyplot.xlabel(r'$\log_{10}(|\Delta e|_{\rm tidal}/|\Delta e|_{\rm enc})$')
	plot_de.hist (de_ratios, histtype='step', bins=10)

	pyplot.tight_layout()
	pyplot.savefig(root_dir+"de-ratios-"+type+".pdf")
	pyplot.clf()