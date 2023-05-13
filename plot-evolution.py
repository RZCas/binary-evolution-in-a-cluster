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

# a_out = 3
# m_total = 1e6
# b = 1
# A_ast = 0.3 #A_* for Hernquist potential
# potential = "Hernquist"
# pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc)
# potential = "Plummer"
# pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

types = ['mtotal=1e5', 'mtotal=1e5,aout=3', 'perpendicular-hard-hernquist-light', 'wide_range_mtotal=1e7_hernquist', 'perpendicular-hard', 'wide_range_mtotal=1e6_plummer', 'perpendicular-soft-hernquist']
m_totals = [1e5, 1e5, 1e6, 1e6, 1e6, 1e6, 1e6]
nokicks = [False, False, False, False, False, False]
# if 'nokick' in types[0]:
# 	nokicks = True
# else:
# 	nokicks = False
a_fixed = [2]
potential_types = ['Hernquist', 'Hernquist', 'Hernquist', 'Hernquist', 'Hernquist', 'Plummer', 'Hernquist']
b = [2, 2, 1, 2, 1, 2, 1]
# if 'b=1' in types[0]: 
# 	b = 1
# else:
# 	b = 2
def potential (type, m_tot, b):
	if type=='Plummer':
		return PlummerPotential(amp=m_tot*u.solMass, b=b*u.pc) 
	else:
		return HernquistPotential(amp=2*m_tot*u.solMass, a=b*u.pc)
potentials = [potential (potential_types[i], m_totals[i], b[i]) for i in range(len(types))]

for i in [0]:#range(len(types)):
	type = types[i]
	potential = potentials[i]
	m_total = m_totals[i]
	subfolder = 'm1=m2=10/'
	root_dir = "output/"+subfolder+type+"/"
	# for filepath in glob.glob(root_dir+'*.txt'):
	for index in [12,32,83,147,187,248,270,297,343,472,584]:#range(10):
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

			t_no_doubles = [0]
			de_abs_tidal = [0]
			de_abs_flybys = [0]
			t_strong_flybys = [0]
			t_weak_flybys = [0]
			de_abs_strong_flybys = [0]
			de_abs_weak_flybys = [0]

			# first- and second-order terms:
			de_tidal = [0]
			de_flybys = [0]
			de_strong_flybys = [0]
			de2_tidal = [0]
			de2_flybys = [0]
			de2_strong_flybys = [0]

			dedt_tidal = []
			dedt_flybys = []

			with open(filepath) as f:
				for line in f:
					lineNumber+=1
					data = line.split()
					if len(data) > 1:
						# if data[0]=="perturber:":
							# startLineNumber = lineNumber + 1
							# if lineNumber>2: 
							# 	Q = float(data[2])
							# if lineNumber%3==0:
							# 	shift=1
						if isfloat(data[0]) and isfloat(data[1]):# and (len(data)==15 or len(data)==21 or lineNumber==3):
							# if len(de2_tidal)-len(de2_flybys)>2:
								# print(lineNumber) 
							try:
								t_0 = float(data[0])/1e9
								if len(t)>0 and t_0<t[-1]:								
									print("time reversed", index, lineNumber)
									continue
								R = float(data[1])
								z = float(data[2])
								phi = float(data[3])
								v_R = float(data[4])
								v_z = float(data[5])
								v_phi = float(data[6])
								a_0 = float(data[7])
								m = float(data[8])
								q = float(data[9])
								e_0 = float(data[10])
								if e_0>1:
									print("e =", e_0, index, lineNumber)
								i_0 = float(data[11])
								if len(data)>17:
									epsilon = float(data[17])
									if epsilon>0:# and E<0:
										t_logepsilon.append(t_0)
										logepsilon.append(np.log10(epsilon))
							except:
								print("bad data at", index, lineNumber)
								continue
							# if len(t)>0 and t_0<t[-1]:	
							# 	print('how\'s that possible?')
							if t_0 > t_max/1e9:	break
							if lineNumber == 3:
								a_initial = a_0
							v = np.sqrt(v_R**2+v_z**2+v_phi**2)
							V.append(v)
							v_array.append(v)
							a.append(a_0)
							r_p.append(a_0*(1-e_0))
							t.append(t_0)
							theta.append((1-e_0**2)*np.cos(i_0)**2)
							e.append(e_0)
							cosi.append(np.cos(i_0))
							if nokicks[i]:
								E = -1
								ra, rp = a_fixed[i], a_fixed[i]
							else:
								E = (v/_kms)**2/2 + evaluatePotentials(potential, R/_pc, z/_pc, phi=phi, use_physical=False) 
								if E<0:
									ra, rp = rarp(potential, [R, z, phi], [v_R, v_z, v_phi])
								else:
									ra, rp = 0, 0
								# ra, rp = np.sqrt(R**2+z**2), np.sqrt(R**2+z**2)
							ra_array.append(ra)
							rp_array.append(rp)
							t_rarp.append(t_0)
							if lineNumber > 3 and abs(m - m_prev)>1e-5:
								color = 'm'
								exchange.append(t_0)
							m_prev = m
							if lineNumber == 3:
								m_prev = m
								# E_0 = E
								a_i = float(data[7])
								m1 = m/(1+q)
								m2 = m*q/(1+q)
							# if lineNumber == startLineNumber: 
							# 	E_prev = E
							# if lineNumber == startLineNumber + 1:
							# 	t_dE.append(t_0)
							# 	dE_total += E - E_prev
							# 	dE_total_dE_0.append(np.log10(abs(dE_total/E_0)))
							# if lineNumber%3==0+shift:
							# 	t_previous = t_0
							# if lineNumber%3==1+shift and lineNumber>1:
							# 	t_no_doubles.append(t_0)
							# 	# de_abs_tidal.append(de_abs_tidal[-1]+float(data[19]))
							# 	de_tidal.append(de_tidal[-1]+e[-1]-e[-2])
							# 	de2_tidal.append(de2_tidal[-1]+(e[-1]-e[-2])**2)
							# 	if t_0==t_previous:
							# 		print('hmm, that\'s bad...', index, lineNumber)
									# dedt_tidal.append(dedt_tidal[-1])
								# else:
									# dedt_tidal.append(float(data[19])/(t_0-t_previous))
									# if lineNumber>1950 and lineNumber<1970:
									# 	print(lineNumber, "de =", de_abs_tidal[-1], "dt =", t_0-t_previous)
									# if len(dedt_tidal)>=2 and dedt_tidal[-1]>1e-4 and dedt_tidal[-2]<1e-4:
									# 	print("something happened here:", lineNumber)
							# if lineNumber%3==0+shift and lineNumber>3:
							# 	# de_abs_flybys.append(de_abs_flybys[-1]+abs(e[-1]-e[-2]))
							# 	de_flybys.append(de_flybys[-1]+e[-1]-e[-2])
							# 	de2_flybys.append(de2_flybys[-1]+(e[-1]-e[-2])**2)
							# 	if Q/a_0 < 25:
							# 		t_strong_flybys.append(t_0)
							# 		# de_abs_strong_flybys.append(de_abs_strong_flybys[-1]+abs(e[-1]-e[-2]))
							# 		de_strong_flybys.append(de_strong_flybys[-1]+e[-1]-e[-2])
							# 		de2_strong_flybys.append(de2_strong_flybys[-1]+(e[-1]-e[-2])**2)
							# E_array.append(E)
						elif data[1] == 'calculation':	#N-body calculation abandoned
							t_prev = float(data[0])
							# de_abs_flybys.append(de_abs_flybys[-1])
							de_flybys.append(de_flybys[-1])
							de2_flybys.append(de2_flybys[-1])
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

			figure = pyplot.figure(figsize=(18, 16))
			# figure.suptitle(fr'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_0$ = {a_initial:.1f} AU, \\$a_\mathrm{{h}}$ = {a_h(m1, m2, a_fixed[0], type=potential, m_total=m_total, b=b):.1f} AU, $a(\epsilon_\mathrm{{GR}}=1)$ = {a_tidal(m1+m2, m_total, b):.1f} AU, $a(t_\mathrm{{sec}}=0.1t_\mathrm{{H}})$ = {a_tsec01tH(m1+m2, m_total, b):.1f} AU', fontsize=24)
			figure.suptitle(fr'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_0$ = {a_initial:.1f} AU' + f'\n{result}', fontsize=24)
			# pyplot.figtext(0.3, 0.93, f'{result}', fontsize=24, color=color)

			# figure.suptitle(r'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_0$ = {a_0:.1f} AU, \\$a_\mathrm{{h}}$ = {a_h:.1f} AU, $a(\epsilon_\mathrm{{GR}}=1)$ = {a_tidal:.1f} AU, $a(t_\mathrm{{sec}}=0.1t_\mathrm{{H}})$ = {a_tsec01tH:.1f} AU\\ {{\textcolor{red}{result}}}'.format(m1=m1, m2=m2, a_0=a_i, a_h=a_h(m1, m2, a_out, type=potential, m_total=m_total, b=b), a_tidal=a_tidal(m1+m2, m_total, b), a_tsec01tH=a_tsec01tH(m1+m2, m_total, b), result=result, color=color), fontsize=24)

			# n = 30
			# encounter_rate, bin_edges = np.histogram (t, bins=n)
			# encounter_rate = encounter_rate.astype(np.double)
			# bin_width = bin_edges[1]-bin_edges[0]
			# bin_centres = []
			# for i in range(n):
			# 	bin_centres.append((bin_edges[i]+bin_edges[i+1])/2)
			# 	# encounter rate in Myr^-1 (factor of 2 because there are 2 entries for every encounter)
			# 	# print(i, encounter_rate[i])
			# 	encounter_rate[i] = encounter_rate[i] / (2*1000*bin_width)
			# 	# print(i, encounter_rate[i])

			# normalize the encounter rate to the initial semimajor axis
			# i = 0
			# for i_t in range(len(t)):
			# 	if t[i_t]>bin_centres[i]:
			# 		encounter_rate[i] *= tau_0_factor (a[i_t]|units.AU, (m1+m2)|units.MSun, r[i_t]|units.pc, Q_max_a=50, type=potential, m_total=m_total, b=b, V=V[i_t]|units.kms) / tau_0_factor (a[0]|units.AU, (m1+m2)|units.MSun, r[0]|units.pc, Q_max_a=50, type=potential, m_total=m_total, b=b, V=V[0]|units.kms)
			# 		# print(i, encounter_rate[i])
			# 		i += 1
			# 		if i>=n: break

			color = 'k'

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
			# plot_r.set_ylim(0,10)
			ax.plot(t_rarp, ra_array, color)
			ax.plot(t_rarp, rp_array, color)
			# ax.plot(t, r, 'r')

			# ax1 = ax.twinx()
			# ax1.minorticks_on() 
			# ax1.tick_params(labelsize=14)
			# ax1.set_ylabel(r'normalized encounter rate [Myr$^{-1}$]', fontsize=16)
			# ax1.plot(bin_centres, encounter_rate)

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

			# plot_de = figure.add_subplot(4,2,7)
			# ax = pyplot.gca()
			# ax.minorticks_on() 
			# ax.tick_params(labelsize=14)
			# ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			# ax.set_ylabel(r'$|\Delta e|$', fontsize=16)
			# plot_de.plot(t_no_doubles, de_abs_tidal, color, label=r'total $|\Delta e|$ due to tidal effects')
			# if len(t_no_doubles)==len(de_abs_flybys):
			# 	plot_de.plot(t_no_doubles, de_abs_flybys, color+'--', label=r'total $|\Delta e|$ due to flybys')
			# else:
			# 	plot_de.plot(t_no_doubles[:-1], de_abs_flybys, color+'--', label=r'total $|\Delta e|$ due to flybys')
			# plot_de.plot(t_strong_flybys, de_abs_strong_flybys, color+':', label=r'$Q/a<25$ only')
			# pyplot.yscale('log')

			# de_flybys = np.array(de_flybys)
			# de_tidal = np.array(de_tidal)
			# de_strong_flybys = np.array(de_strong_flybys)
			# plot_de = figure.add_subplot(4,2,7)
			# ax = pyplot.gca()
			# ax.minorticks_on() 
			# ax.tick_params(labelsize=14)
			# ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			# ax.set_ylabel(r'$\Delta e$', fontsize=16)
			# if len(t_no_doubles) != len(de_flybys):
			# 	t_no_doubles_new = t_no_doubles[:-1]
			# else:
			# 	t_no_doubles_new = t_no_doubles

			# # plot the negative values with a different linestyle
			# start_index = 1
			# end_index = 1
			# current_sign = np.sign(de_flybys[1])
			# while end_index<len(de_flybys):
			# 	while current_sign == np.sign(de_flybys[end_index]):
			# 		end_index += 1
			# 		if end_index==len(de_flybys):
			# 			break
			# 	if current_sign == 1:
			# 		plot_de.plot(t_no_doubles_new[start_index:end_index], abs(de_flybys[start_index:end_index]), 'r', label=r'total $\Delta e$ due to flybys')
			# 	else: 
			# 		plot_de.plot(t_no_doubles_new[start_index:end_index], abs(de_flybys[start_index:end_index]), 'r--', label=r'total $\Delta e$ due to flybys, $\Delta e<0$')
			# 	if end_index<len(de_flybys):
			# 		start_index = end_index
			# 		current_sign = np.sign(de_flybys[end_index])

			# start_index = 1
			# end_index = 1
			# current_sign = np.sign(de_tidal[1])
			# while end_index<len(de_tidal):
			# 	while current_sign == np.sign(de_tidal[end_index]):
			# 		end_index += 1
			# 		if end_index==len(de_tidal):
			# 			break
			# 	if current_sign == 1:
			# 		plot_de.plot(t_no_doubles[start_index:end_index], de_tidal[start_index:end_index], 'k', label=r'total $\Delta e$ due to tidal effects')
			# 	else: 
			# 		plot_de.plot(t_no_doubles[start_index:end_index], -de_tidal[start_index:end_index], 'k--', label=r'total $\Delta e$ due to tidal effects, $\Delta e<0$')	
			# 	if end_index<len(de_tidal):
			# 		start_index = end_index
			# 		current_sign = np.sign(de_tidal[end_index])

			# # plot_de.plot(t_strong_flybys_new, abs(de_strong_flybys), color+':', label=r'$Q/a<25$ only')
			# handles, labels = pyplot.gca().get_legend_handles_labels()
			# by_label = dict(zip(labels, handles))
			# pyplot.legend(by_label.values(), by_label.keys())
			# # plot_de.legend()
			# pyplot.yscale('log')

			# plot_de2 = figure.add_subplot(4,2,8)
			# ax = pyplot.gca()
			# ax.minorticks_on() 
			# ax.tick_params(labelsize=14)
			# ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			# ax.set_ylabel(r'$(\Delta e)^2$', fontsize=16)
			# # print(len(t_no_doubles), len(de2_tidal), len(t_no_doubles_new), len(de2_flybys))
			# plot_de2.plot(t_no_doubles, de2_tidal, 'k', label=r'total $(\Delta e)^2$ due to tidal effects')
			# plot_de2.plot(t_no_doubles_new, de2_flybys, 'r', label=r'total $(\Delta e)^2$ due to flybys')
			# plot_de2.plot(t_strong_flybys, de2_strong_flybys, 'r--', label=r'$Q/a<25$ only')
			# plot_de2.legend()
			# pyplot.yscale('log')

			# for i in range(0,len(t_no_doubles)):
			# 	if t_no_doubles[i] > t_no_doubles[-1]/10:
			# 		cutoffNumber = i
			# 		break

			# plot_de.set_ylim(bottom=de_abs_tidal[cutoffNumber]/10)
			# plot_de.legend()

			# plot_dedt = figure.add_subplot(4,2,8)
			# ax = pyplot.gca()
			# ax.minorticks_on() 
			# ax.tick_params(labelsize=14)
			# ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
			# ax.set_ylabel(r'$|de/dt|_{\rm tidal}$', fontsize=16)
			# if len(t_no_doubles)==len(dedt_tidal):
			# 	plot_dedt.plot(t_no_doubles, dedt_tidal, color)
			# else:
			# 	plot_dedt.plot(t_no_doubles[:-1], dedt_tidal, color)
			# pyplot.yscale('log')

			pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
			# index = os.path.split(filepath)[1][:-4]
			pyplot.savefig(root_dir+"evolution-"+type+"-"+str(index)+".pdf")
			pyplot.clf()