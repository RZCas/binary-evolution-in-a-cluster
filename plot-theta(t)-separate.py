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
a_out = 1.6
m_total = 1e6
b = 1
A_ast = 0.3 #A_* for Hernquist potential
# print(a_h(10, 10, a_out, type="Hernquist", m_total=1e6, b=1))

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

type = 'spherical-noweak'
root_dir = "output/"+type+"/"
for index in range(20):
	# index=0
	filepath = root_dir + str(index) + '.txt'
	color = 'k'
	t = []
	theta = []
	e = []
	cosi = []
	a = []
	t_0 = 0
	theta_previous = 0
	lineNumber=0
	with open(filepath) as f:
		for line in f:
			lineNumber+=1
			data = line.split()
			if len(data) > 1:
				if lineNumber == 3:
					m = float(data[8])
					q = float(data[9])
					m1 = m/(1+q)
					m2 = m*q/(1+q)
				if isfloat(data[0]) and isfloat(data[1]):
					t_0 = float(data[0])
					if t_0 > t_max:	break
					# print(filepath, line)
					a_0 = float(data[7])
					e_0 = float(data[10])
					i_0 = float(data[11])
					a.append(a_0)
					t.append(t_0)
					theta.append((1-e_0**2)*np.cos(i_0)**2)
					e.append(e_0)
					cosi.append(np.cos(i_0))
					# if (theta[-1]-theta_previous < -0.15):
					# 	print(lineNumber)
					# theta_previous = theta[-1]
				elif data[1] == 'destroyed': color = 'r'
				elif data[1] == 'merger': color = 'g'

	figure = pyplot.figure(figsize=(18, 12))
	figure.suptitle(r'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_\mathrm{{h}}$ = {a_h:.1f} AU, $a(\epsilon_\mathrm{{GR}}=1)$ = {a_tidal:.1f} AU, $a(t_\mathrm{{sec}}=0.1t_\mathrm{{H}})$ = {a_tsec01tH:.1f} AU'.format(m1=m1, m2=m2, a_h=a_h(m1, m2, a_out, type=potential, m_total=m_total, b=b), a_tidal=a_tidal(m1+m2, m_total, b), a_tsec01tH=a_tsec01tH(m1+m2, m_total, b)), fontsize=24)

	plot_theta = figure.add_subplot(2,2,1)
	ax = pyplot.gca()
	ax.minorticks_on() 
	ax.tick_params(labelsize=14)
	ax.set_xlabel(r'$t$ [yr]', fontsize=16)
	ax.set_ylabel(r'$\Theta$', fontsize=16)
	plot_theta.plot(t, theta, color)

	plot_e = figure.add_subplot(2,2,2)
	ax = pyplot.gca()
	ax.minorticks_on() 
	ax.tick_params(labelsize=14)
	ax.set_xlabel(r'$t$ [yr]', fontsize=16)
	ax.set_ylabel(r'$e$', fontsize=16)
	plot_e.plot(t, e, color)

	plot_cosi = figure.add_subplot(2,2,3)
	ax = pyplot.gca()
	ax.minorticks_on() 
	ax.tick_params(labelsize=14)
	ax.set_xlabel(r'$t$ [yr]', fontsize=16)
	ax.set_ylabel(r'$\cos{i}$', fontsize=16)
	plot_cosi.plot(t, cosi, color)

	plot_a = figure.add_subplot(2,2,4)
	ax = pyplot.gca()
	ax.minorticks_on() 
	ax.tick_params(labelsize=14)
	ax.set_xlabel(r'$t$ [yr]', fontsize=16)
	ax.set_ylabel(r'$a$ [AU]', fontsize=16)
	plot_a.plot(t, a, color)

	pyplot.tight_layout()
	pyplot.savefig(root_dir+"evolution-"+type+"-"+str(index)+".pdf")