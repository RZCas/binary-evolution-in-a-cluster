import os
import glob
import numpy as np
import astropy.units as u
from galpy.potential import evaluaterforces, evaluatePotentials, PotentialError, KeplerPotential, TwoPowerTriaxialPotential, PlummerPotential, HernquistPotential
from binary_evolution_with_flybys import a_h
from amuse.lab import units, constants 
from binary_evolution.tools import rarp
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

t_max=1e80
a_out = 1.6
m_total = 1e6
b = 1
A_ast = 0.3 #A_* for Hernquist potential
pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc)

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

# out = open('output/rarp.txt', 'w')

# types = ['perpendicular-hard', 'perpendicular-veryhard', 'perpendicular-noweak-hard', 'perpendicular-noweak-veryhard', 'perpendicular', 'perpendicular-noweak']
types = ['perpendicular-hard']
for type in types:
	root_dir = "output/"+type+"/"
	for index in range(19):
		# if type=='perpendicular-hard' and index==4:
		if True:#index<3:
			filepath = root_dir + str(index) + '.txt'
			color = 'k'
			t = []
			theta = []
			e = []
			cosi = []
			a = []
			r = []
			dE_total_dE_0 = []
			t_dE = []
			E_array = []
			v_array = []
			logepsilon = []
			t_logepsilon = []
			r_p = []
			exchange = []
			ra_array = []
			rp_array = []
			t_rarp = []
			dE_total = 0
			t_0 = 0
			theta_previous = 0
			lineNumber = 0
			with open(filepath) as f:
				for line in f:
					lineNumber+=1
					data = line.split()
					if len(data) > 1:
						if data[0]=="perturber:":
							startLineNumber = lineNumber + 1
						if isfloat(data[0]) and isfloat(data[1]):
							t_0 = float(data[0])/1e9
							if t_0 > t_max/1e9:	break
							R = float(data[1])
							z = float(data[2])
							phi = float(data[3])
							r_0 = np.sqrt(R**2+z**2)
							v_R = float(data[4])
							v_z = float(data[5])
							v_phi = float(data[6])
							v = np.sqrt(v_R**2+v_z**2+v_phi**2)
							a_0 = float(data[7])
							e_0 = float(data[10])
							i_0 = float(data[11])
							r.append(r_0)
							v_array.append(v)
							a.append(a_0)
							r_p.append(a_0*(1-e_0))
							t.append(t_0)
							theta.append((1-e_0**2)*np.cos(i_0)**2)
							e.append(e_0)
							cosi.append(np.cos(i_0))
							E = (v/_kms)**2/2 + evaluatePotentials(pot, R/_pc, z/_pc, phi=phi, use_physical=False) 
							# print((v/_kms)**2/2, evaluatePotentials(pot, R/_pc, z/_pc, phi=phi, use_physical=False) )
							if E<0:
								ra, rp = rarp(pot, [R, z, phi], [v_R, v_z, v_phi])
							else:
								ra, rp = 0, 0
							# print(t_0, ra, rp, r_0, file=out)
							if ra>0 and rp>0:
								ra_array.append(ra)
								rp_array.append(rp)
								t_rarp.append(t_0)
							m = float(data[8])
							q = float(data[9])
							if lineNumber > 3 and abs(m - m_prev)>1e-5:
								color = 'm'
								exchange.append(t_0)
							m_prev = m
							if lineNumber == 3:
								m_prev = m
								E_0 = E
								a_i = float(data[7])
								m1 = m/(1+q)
								m2 = m*q/(1+q)
							if lineNumber == startLineNumber: 
								E_prev = E
							if lineNumber == startLineNumber + 1:
								t_dE.append(t_0)
								dE_total += E - E_prev
								dE_total_dE_0.append(np.log10(abs(dE_total/E_0)))
							# if abs(E/E_prev-1)>0.1:
								# print(lineNumber)
							E_array.append(E)
							if len(data)>17:
								epsilon = float(data[17])
								if epsilon>0 and E<0:
									t_logepsilon.append(t_0)
									logepsilon.append(np.log10(epsilon))
						elif data[1] == 'destroyed': color = 'r'
						elif data[1] == 'merger': color = 'g'

			figure = pyplot.figure(figsize=(18, 15))
			figure.suptitle(r'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_0$ = {a_0:.1f} AU, \\$a_\mathrm{{h}}$ = {a_h:.1f} AU, $a(\epsilon_\mathrm{{GR}}=1)$ = {a_tidal:.1f} AU, $a(t_\mathrm{{sec}}=0.1t_\mathrm{{H}})$ = {a_tsec01tH:.1f} AU'.format(m1=m1, m2=m2, a_0=a_i, a_h=a_h(m1, m2, a_out, type=potential, m_total=m_total, b=b), a_tidal=a_tidal(m1+m2, m_total, b), a_tsec01tH=a_tsec01tH(m1+m2, m_total, b)), fontsize=24)

			plot_theta = figure.add_subplot(3,2,1)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$\Theta$', fontsize=16)
			plot_theta.plot(t, theta, color)

			plot_e = figure.add_subplot(3,2,2)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$e$', fontsize=16)
			plot_e.plot(t, e, color)
			for exchange_time in exchange:
				plot_e.plot([exchange_time,exchange_time], [0,1], 'k--')

			plot_cosi = figure.add_subplot(3,2,3)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$\cos{i}$', fontsize=16)
			plot_cosi.plot(t, cosi, color)

			plot_a = figure.add_subplot(3,2,4)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$a$ [AU], $r_p$ [AU]', fontsize=16)
			ax.set_yscale('log')
			plot_a.plot(t, a, color)
			plot_a.plot(t, r_p, color+'--')

			plot_r = figure.add_subplot(3,2,5)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$r$ [pc]', fontsize=16)
			# plot_r.set_xlim(5,6)
			# plot_rs.set_ylim(0,10)
			plot_r.plot(t_rarp, ra_array, color)
			plot_r.plot(t_rarp, rp_array, color)

			# plot_dE = figure.add_subplot(3,2,6)
			# ax = pyplot.gca()
			# ax.minorticks_on() 
			# ax.tick_params(labelsize=14)
			# ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			# # ax.set_ylabel(r'$\log_{10}{|\Delta E/E_0|}$', fontsize=16)
			# ax.set_ylabel(r'$E$', fontsize=16)
			# # plot_dE.plot(t_dE, dE_total_dE_0, color)
			# # plot_dE.set_ylim(bottom=-10)
			# plot_dE.plot(t, E_array, color)

			# plot_v = figure.add_subplot(3,2,6)
			# ax = pyplot.gca()
			# ax.minorticks_on() 
			# ax.tick_params(labelsize=14)
			# ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			# ax.set_ylabel(r'$v$ [km/s]', fontsize=16)
			# plot_v.plot(t, v_array, color)

			plot_epsilon = figure.add_subplot(3,2,6)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$\log_{10}\epsilon$', fontsize=16)
			plot_epsilon.plot(t_logepsilon, logepsilon, color)
			plot_epsilon.plot([0, t[-1]], [np.log10(20), np.log10(20)], 'r')

			pyplot.tight_layout()
			pyplot.savefig(root_dir+"evolution-"+type+"-"+str(index)+".pdf")
			pyplot.clf()