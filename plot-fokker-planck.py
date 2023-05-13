import os
import glob
import numpy as np
import astropy.units as u
from galpy.potential import evaluaterforces, evaluatePotentials, PotentialError, KeplerPotential, TwoPowerTriaxialPotential, PlummerPotential, HernquistPotential
from binary_evolution_with_flybys import a_h, sigma
from amuse.lab import units, constants 
from binary_evolution.tools import rarp
from statistics import mean
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

def calculate_chi (t, de_tidal, de_flybys, n):
	# cut the timespan into n equal pieces and calculate chi for each one
	chi_n = []
	t_chi_n = [(i+0.5)*t[-1]/n for i in range(n)]
	chi_n_bin = 0
	i_begin = 0
	for i in range(len(t)):
		chi_n_bin_new = int(t[i]/t[-1]*n)
		if chi_n_bin_new > chi_n_bin:
			i_end = i
			de_t_max = max(abs(de_tidal[i_begin:i_end]-de_tidal[i_begin]))
			de_f_max = max(abs(de_flybys[i_begin:i_end]-de_flybys[i_begin]))
			# if de_t_max == 0:
			# 	print(de_tidal[i_begin:i_end], i_begin, i_end, len(t))
			chi_n.append(de_t_max / (de_t_max + de_f_max))
			i_begin = i_end
			chi_n_bin += 1
	return t_chi_n, chi_n

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
m_total = 1e5
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

# out = open('output/rarp.txt', 'w')

# types = ['perpendicular-soft-hernquist', 'perpendicular-hard-hernquist']
# types = ['perpendicular-hard-plummer', 'perpendicular-hard-plummer-light']
# types = ['wide_range_1','wide_range_mtotal=1e5','wide_range_m=60','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e6_nokicks_plummer','wide_range_mtotal=1e6_nokicks']
# types = ['wide_range_1','wide_range_mtotal=1e6_nokicks']
# types = ['wide_range_mtotal=1e5','wide_range_mtotal=1e5_nokicks']
# types = ['wide_range_mtotal=1e6_nokicks_plummer']
nokicks = True
# types = ['mtotal=1e6_nokicks', 'uniform_mtotal=1e5_plummer', 'uniform_mtotal=1e5_hernquist', 'uniform_mtotal=1e6_hernquist']
types = ['hernquist,m_total=1e5,b=1,a_out=3,i=89.9,nokicks', 'hernquist,m_total=1e5,b=1,a_out=4,i=89.9,nokicks,a_in=300,ns', 'hernquist,m_total=1e5,b=1,a_out=4,i=89.9,nokicks,a_in=300', 'uniform_mtotal=1e5_hernquist']
a_fixed = [3, 4, 4, 0]
potentials = ['Hernquist', 'Hernquist', 'Hernquist', 'Hernquist']
m_totals = [1e5, 1e5, 1e5, 1e5]
for i in [3]:
# for type in types:
	type = types[i]
	potential = potentials[i]
	m_total = m_totals[i]
	root_dir = "output/"+type+"/"
	# for filepath in glob.glob(root_dir+'*.txt'):
	for index in range(20):
		if True:
			shift = 0
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

			de_abs_tidal = [0]
			de2_flybys = [0]

			with open(filepath) as f:
				for line in f:
					lineNumber+=1
					data = line.split()
					if len(data) > 1:
						if data[0]=="perturber:":
							startLineNumber = lineNumber + 1
							if lineNumber>2: 
								Q = float(data[2])
							if lineNumber%3==0:
								shift=1
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
							if lineNumber == 3:
								a_initial = a_0
								e_initial = e_0
							i_0 = float(data[11])
							r.append(r_0)
							v_array.append(v)
							a.append(a_0)
							r_p.append(a_0*(1-e_0))
							t.append(t_0)
							theta.append((1-e_0**2)*np.cos(i_0)**2)
							e.append(e_0)
							cosi.append(np.cos(i_0))
							# if nokicks:
							# 	E = -1
							# 	ra, rp = a_fixed[i], a_fixed[i]
							# else:
							# 	E = (v/_kms)**2/2 + evaluatePotentials(pot, R/_pc, z/_pc, phi=phi, use_physical=False) 
							# 	if E<0:
							# 		ra, rp = rarp(pot, [R, z, phi], [v_R, v_z, v_phi])
							# 	else:
							# 		ra, rp = 0, 0
							# ra_array.append(ra)
							# rp_array.append(rp)
							# t_rarp.append(t_0)
							m = float(data[8])
							q = float(data[9])
							if lineNumber > 3 and abs(m - m_prev)>1e-5:
								color = 'm'
								exchange.append(t_0)
							m_prev = m
							if lineNumber == 3:
								m_prev = m
								a_i = float(data[7])
								m1 = m/(1+q)
								m2 = m*q/(1+q)
							if lineNumber%3==0+shift:
								t_previous = t_0
							if lineNumber%3==1+shift and lineNumber>1:
								de_abs_tidal.append(de_abs_tidal[-1]+abs(e[-1]-e[-2]))
								de_abs_total.append(de_abs_total[-1]+abs(e[-1]-e[-2]))
								de_abs_flybys.append(de_abs_flybys[-1])
								de_tidal.append(de_tidal[-1]+e[-1]-e[-2])
								de_total.append(e[-1]-e_initial)
								de_flybys.append(de_flybys[-1])
								de2_flybys.append(de2_flybys[-1])

								de_tidal_max.append(max(abs(de_tidal[-1]),de_tidal_max[-1]))
								de_flybys_max.append(max(abs(de_flybys[-1]),de_flybys_max[-1]))

								if t_0==t_previous:
									print('hmm, that\'s bad...', index, lineNumber)
							if lineNumber%3==0+shift and lineNumber>3:
								de_abs_flybys.append(de_abs_flybys[-1]+abs(e[-1]-e[-2]))
								de_abs_total.append(de_abs_total[-1]+abs(e[-1]-e[-2]))
								de_abs_tidal.append(de_abs_tidal[-1])
								de_flybys.append(de_flybys[-1]+e[-1]-e[-2])
								de2_flybys.append(de2_flybys[-1]+(e[-1]-e[-2])**2)
								de_total.append(e[-1]-e_initial)
								de_tidal.append(de_tidal[-1])
								
								de_tidal_max.append(max(abs(de_tidal[-1]),de_tidal_max[-1]))
								de_flybys_max.append(max(abs(de_flybys[-1]),de_flybys_max[-1]))

								de_t = abs(e[-2]-e[-3])
								de_f = abs(e[-1]-e[-2])
								t_xi_2.append(t_0)
								# de_ratios.append(de_t/(de_t+de_f))
								# if lineNumber<100 or lineNumber>21600:
								# 	print(lineNumber, de_t/(de_t+de_f), de_t)
								if xi_2 == []:
									xi_2.append(de_t/(de_t+de_f))
								else:
									xi_2.append((xi_2[-1]*len(xi_2)+de_t/(de_t+de_f))/(len(xi_2)+1))
							if len(data)>17:
								epsilon = float(data[17])
								if epsilon>0:# and E<0:
									t_logepsilon.append(t_0)
									logepsilon.append(np.log10(epsilon))
						elif data[1] == 'calculation':	#N-body calculation abandoned
							t_prev = float(data[0])
							# de_abs_flybys.append(de_abs_flybys[-1])
							# de_flybys.append(de_flybys[-1])
							# de2_flybys.append(de2_flybys[-1])
						elif data[1] == 'destroyed': 
							color = 'r'
							result = 'binary destroyed'
						elif data[1] == 'merger': 
							color = 'g'
							result = 'binary merged'
						elif data[1] == 'maximum' and data[2] == 'semimajor': 
							color = 'b'
							result = 'calculation adandoned (semimajor axis too large)'
						elif data[1] == 'ejected':
							color = 'tab:purple'
							result = 'binary ejected from the cluster'

			de_abs_tidal = np.array(de_abs_tidal)
			de_abs_total = np.array(de_abs_total)
			xi_1 = de_abs_tidal[1:] / de_abs_total[1:]

			chi = np.array(de_tidal_max[1:]) / (np.array(de_tidal_max[1:]) + np.array(de_flybys_max[1:]))		

			t_chi_100, chi_100 = calculate_chi (t, np.array(de_tidal), np.array(de_flybys), 100)
			t_chi_30, chi_30 = calculate_chi (t, np.array(de_tidal), np.array(de_flybys), 30)

			figure = pyplot.figure(figsize=(12, 12))
			figure.suptitle(rf'$\xi_1$ = {xi_1[-1]:.3f}, $\xi_2$ = {xi_2[-1]:.3f}, $\chi$ = {chi[-1]:.3f}', fontsize=16)
			gs = figure.add_gridspec(4, 1)#, hspace=0, wspace=0)
			ax1, ax2, ax3, ax4 = gs.subplots(sharex=True)

			ax1.minorticks_on() 
			ax1.tick_params(labelsize=14)
			ax1.set_ylabel(r'$a$ [AU], $r_p$ [AU]', fontsize=16)
			ax1.set_yscale('log')
			ax1.plot(t, a, 'k')
			ax1.plot(t, r_p, 'k--')

			ax2.minorticks_on()     
			ax2.tick_params(labelsize=14)
			ax2.set_ylabel(r'$\Delta e$', fontsize=16)
			ax2.plot([0,t[-1]], [0,0], 'k--')
			ax2.plot(t, de_flybys, 'r', label='flybys')
			ax2.plot(t, de_tidal, 'b', label='tidal')
			ax2.plot(t, de_total, 'k', label='total')
			ax2.legend(fontsize=16)

			ax3.minorticks_on()     
			ax3.tick_params(labelsize=14)
			ax3.set_ylabel(r'$\sum|\delta e|$', fontsize=16)
			ax3.plot(t, de_abs_flybys, 'r', label='flybys')
			ax3.plot(t, de_abs_tidal, 'b', label='tidal')
			ax3.plot(t, de_abs_total, 'k', label='total')
			ax3.legend(fontsize=16)
			
			ax4.minorticks_on()     
			ax4.tick_params(labelsize=14)
			ax4.set_ylabel(r'$\xi_1$, $\xi_2$, $\chi$', fontsize=16)
			ax4.plot([0,t[-1]], [0,0], 'k--')
			# ax4.plot(t[1:], xi_1, 'r', label=r'$\xi_1$')
			# ax4.plot(t_xi_2, xi_2, 'b', label=r'$\xi_2$')
			# ax4.plot(t_chi_100, chi_100, 'r', label=r'$\chi_{100}$')
			ax4.plot(t_chi_30, chi_30, 'b', label=r'$\chi_{30}$')
			ax4.plot(t[1:], chi, 'k', label=r'$\chi$')
			ax4.legend(fontsize=16)

			# ax5.minorticks_on()     
			# ax5.tick_params(labelsize=14)
			# ax5.set_ylabel(r'$\sum(\delta e_{\rm f})^2$', fontsize=16)
			# ax5.plot(t, de2_flybys, 'k')

			pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
			pyplot.savefig(root_dir+"de_tidal_vs_de_flybys-"+type+"-"+str(index)+".pdf")
			pyplot.clf()