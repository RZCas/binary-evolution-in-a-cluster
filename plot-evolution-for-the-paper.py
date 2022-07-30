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
potential = "Plummer"

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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
# pyplot.yscale('log')

t_max=1e80
a_out = 2
m_total = 1e6
b = 2
A_ast = 0.3 #A_* for Hernquist potential
# pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc)
pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

# out = open('output/rarp.txt', 'w')

# types = ['perpendicular-soft-hernquist', 'perpendicular-hard-hernquist']
# types = ['perpendicular-hard-plummer', 'perpendicular-hard-plummer-light']
# types = ['wide_range_1','wide_range_mtotal=1e5','wide_range_m=60','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e6_nokicks_plummer','wide_range_mtotal=1e6_nokicks']
types = ['wide_range_mtotal=1e6_nokicks']
for type in types:
	root_dir = "output/"+type+"/"
	for index in range(0,1):
		if True:
		# if type == 'perpendicular-hard-plummer-light':
			filepath = root_dir + str(index) + '.txt'
			color = 'k'
			result = 'binary survived'
			t = []
			theta = []
			e = []
			cosi = []
			a = []
			r = []
			V = []
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

			t_no_doubles = [0]
			de_abs_tidal = [0]
			de_abs_flybys = [0]
			t_strong_flybys = [0]
			t_weak_flybys = [0]
			de_abs_strong_flybys = [0]
			de_abs_weak_flybys = [0]

			dedt_tidal = []
			dedt_flybys = []

			with open(filepath) as f:
				for line in f:
					lineNumber+=1
					data = line.split()
					if len(data) > 1:
						if data[0]=="perturber:":
							startLineNumber = lineNumber + 1
							if lineNumber>2: Q = float(data[2])
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
							V.append(v)
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
							# if ra>0 and rp>0:
							ra_array.append(ra)
							rp_array.append(rp)
							t_rarp.append(t_0)
							m = float(data[8])
							q = float(data[9])
							if lineNumber > 3 and abs(m - m_prev)>1e-5:
								# color = 'm'
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
							if lineNumber%3==0:
								t_previous = t_0
							if lineNumber%3==1 and lineNumber>1:
								t_no_doubles.append(t_0)
								de_abs_tidal.append(de_abs_tidal[-1]+float(data[19]))
								if t_0==t_previous:
									print('hmm, that\'s bad...', index, lineNumber)
									dedt_tidal.append(dedt_tidal[-1])
								else:
									dedt_tidal.append(float(data[19])/(t_0-t_previous))
									# if lineNumber>1950 and lineNumber<1970:
									# 	print(lineNumber, "de =", de_abs_tidal[-1], "dt =", t_0-t_previous)
									# if len(dedt_tidal)>=2 and dedt_tidal[-1]>1e-4 and dedt_tidal[-2]<1e-4:
									# 	print("something happened here:", lineNumber)
							if lineNumber%3==0 and lineNumber>3:
								de_abs_flybys.append(de_abs_flybys[-1]+abs(e[-1]-e[-2]))
								if Q/a_0 < 25:
									t_strong_flybys.append(t_0)
									de_abs_strong_flybys.append(de_abs_strong_flybys[-1]+abs(e[-1]-e[-2]))
							E_array.append(E)
							if len(data)>17:
								epsilon = float(data[17])
								if epsilon>0 and E<0:
									t_logepsilon.append(t_0)
									logepsilon.append(np.log10(epsilon))
						elif data[1] == 'calculation':	#N-body calculation abandoned
							t_prev = float(data[0])
							de_abs_flybys.append(de_abs_flybys[-1])
						elif data[1] == 'destroyed': 
							# color = 'r'
							result = 'binary destroyed'
						elif data[1] == 'merger': 
							# color = 'g'
							result = 'binary merged'
						elif data[1] == 'maximum' and data[2] == 'semimajor': 
							# color = 'b'
							result = 'calculation adandoned (semimajor axis too large)'

			figure = pyplot.figure(figsize=(18, 20))
			figure.suptitle(r'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_0$ = {a_0:.1f} AU, \\$a_\mathrm{{h}}$ = {a_h:.1f} AU, $a(\epsilon_\mathrm{{GR}}=1)$ = {a_tidal:.1f} AU, $a(t_\mathrm{{sec}}=0.1t_\mathrm{{H}})$ = {a_tsec01tH:.1f} AU\\Outcome: {result}'.format(m1=m1, m2=m2, a_0=a_i, a_h=a_h(m1, m2, a_out, type=potential, m_total=m_total, b=b), a_tidal=a_tidal(m1+m2, m_total, b), a_tsec01tH=a_tsec01tH(m1+m2, m_total, b), result=result), fontsize=24)

			# normalize the encounter rate to the initial semimajor axis
			i = 0
			for i_t in range(len(t)):
				if t[i_t]>bin_centres[i]:
					encounter_rate[i] *= tau_0_factor (a[i_t]|units.AU, (m1+m2)|units.MSun, r[i_t]|units.pc, Q_max_a=50, type=potential, m_total=m_total, b=b, V=V[i_t]|units.kms) / tau_0_factor (a[0]|units.AU, (m1+m2)|units.MSun, r[0]|units.pc, Q_max_a=50, type=potential, m_total=m_total, b=b, V=V[0]|units.kms)
					# print(i, encounter_rate[i])
					i += 1
					if i>=n: break

			plot_theta = figure.add_subplot(4,2,1)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$\Theta$', fontsize=16)
			plot_theta.plot(t, theta, color)

			plot_e = figure.add_subplot(4,2,2)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$e$', fontsize=16)
			plot_e.plot(t, e, color)
			for exchange_time in exchange:
				plot_e.plot([exchange_time,exchange_time], [0,1], 'k--')

			plot_cosi = figure.add_subplot(4,2,3)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$\cos{i}$', fontsize=16)
			plot_cosi.plot(t, cosi, color)

			plot_a = figure.add_subplot(4,2,4)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$a$ [AU], $r_p$ [AU]', fontsize=16)
			ax.set_yscale('log')
			plot_a.plot(t, a, color)
			plot_a.plot(t, r_p, color+'--')

			plot_r = figure.add_subplot(4,2,5)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$r_{a,out}$ [pc], $r_{p,out}$ [pc]', fontsize=16)
			# plot_r.set_xlim(5,6)
			# plot_r.set_ylim(0,10)
			ax.plot(t_rarp, ra_array, color)
			ax.plot(t_rarp, rp_array, color)
			# ax.plot(t, r, 'r')

			plot_epsilon = figure.add_subplot(4,2,6)
			ax = pyplot.gca()
			ax.minorticks_on() 
			ax.tick_params(labelsize=14)
			ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			ax.set_ylabel(r'$\log_{10}\epsilon$', fontsize=16)
			plot_epsilon.plot(t_logepsilon, logepsilon, color)
			plot_epsilon.plot([0, t[-1]], [np.log10(20), np.log10(20)], 'r')

			pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
			pyplot.savefig(root_dir+"evolution-"+type+"-"+str(index)+".pdf")
			pyplot.clf()