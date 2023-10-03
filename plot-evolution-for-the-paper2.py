import os
import glob
import numpy as np
import astropy.units as u
from galpy.potential import evaluaterforces, evaluatePotentials, PotentialError, KeplerPotential, TwoPowerTriaxialPotential, PlummerPotential, HernquistPotential
# from binary_evolution_with_flybys import a_h, sigma
from amuse.lab import units, constants 
from binary_evolution.tools import rarp
from binary_evolution import KeplerRing
from A_averaged import A_averaged
_pc = 8000
_kms = 220
G = constants.G
c = constants.c
t_H = 1.4e10|units.yr
H = 15
potential = "Hernquist"

m_per = 1|units.MSun
# def tau_0_factor (a, m_bin, r, Q_max_a=50, type="Plummer", m_total=4e6, b=1, V=0|units.kms):
# 	Q_max = Q_max_a * a
# 	v0 = np.sqrt(G*(m_bin+m_per)/Q_max)
# 	sigma_rel = np.sqrt(sigma(r, type, m_total, b)**2 + V**2)
# 	return (Q_max**2*(1+(v0/sigma_rel)**2))**-1

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
# a_out = 2
# m_total = 1e6
# b = 2
# A_ast = 0.3 #A_* for Hernquist potential
# pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc)
# pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

# old:
# perpendicular-soft-hernquist
# 6 - merged
# 16 - abandoned
# 15 - destroyed
# 11 - merged after an exchange

# new:
# m1=m2=10, mtotal=1e6, hernquist
# 0 - abandoned
# 2 - merged after an exchange
# 7 - merged
# 19 - destroyed
# nokicks, 0 - no kicks
# evolution-hernquist,m_total=1e5,b=1,a_out=4,i=89.9,nokicks,a_in=300-5 - the case where tidal effects dominate
# uniform_mtotal=1e5_hernquist-9 - ejected after an exchange
# ejected/wide_range_mtotal=1e5_plummer_b=1-145 - ejected

root_dir = ["output/m1=m2=10/mtotal=1e6/",
			"output/m1=m2=10/mtotal=1e7/"
			]
potential_type = ["Hernquist", "Hernquist"]
mtotal = [1e6, 1e7]
b = [2, 2]
pot = [HernquistPotential(amp=2*1e6*u.solMass, a=2*u.pc),
		HernquistPotential(amp=2*1e7*u.solMass, a=2*u.pc)]
nokicks = [False,False]
indices = [4,4]
chi = [0.14, 0.33]
D1D2 = [0.32, 0.87]
# de_min = [1e-4, 1e-3]
fileNames = ['example1','example2']
for i in range(len(indices)):
	shift = 0
	index = indices[i]
	filepath = root_dir[i] + str(index) + '.txt'
	color = 'k'
	result = 'Binary survived'
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
	de_tidal = [0]
	de_flybys = [0]

	need_to_calculate_R = False
	need_to_calculate_epsilon = False
	t_logepsilon_real = []
	logepsilon_real = []

	filepath_R = "output/for the paper 2/"+fileNames[i]+"-R.txt"
	if os.path.exists(filepath_R):
		with open(filepath_R) as f:
			for line in f:
				data = line.split()
				t_rarp.append(float(data[0]))
				rp_array.append(float(data[1]))
				ra_array.append(float(data[2]))
	else:
		file_R = open(filepath_R, 'w+')
		need_to_calculate_R = True

	with open(filepath) as f:
		for line in f:
			lineNumber+=1
			data = line.split()
			if len(data) > 1:
				if data[0]=="perturber:":
					startLineNumber = lineNumber + 1
					if lineNumber>2: Q = float(data[2])
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
					i_0 = float(data[11])
					Omega_0 = float(data[12])
					omega_0 = float(data[13])
					r.append(r_0)
					v_array.append(v)
					a.append(a_0)
					r_p.append(a_0*(1-e_0))
					t.append(t_0)
					theta.append((1-e_0**2)*np.cos(i_0)**2)
					e.append(e_0)
					cosi.append(np.cos(i_0))
					if lineNumber == 3:
						R_initial = R
						a_initial = a_0
						e_initial = e_0
						i_initial = i_0 * 180/np.pi
						omega_initial = omega_0 * 180/np.pi
						Omega_initial = Omega_0 * 180/np.pi
					if nokicks[i]:
						E = -1
						ra, rp = R_initial, R_initial
						t_rarp.append(t_0)
						ra_array.append(ra)
						rp_array.append(rp)
					elif need_to_calculate_R:
						E = (v/_kms)**2/2 + evaluatePotentials(pot[i], R/_pc, z/_pc, phi=phi, use_physical=False) 
						if E<0:
							ra, rp = rarp(pot[i], [R, z, phi], [v_R, v_z, v_phi])
							if (ra == 0 or rp == 0) and lineNumber<100:
								ra, rp = R_initial, R_initial
						else:
							ra, rp = 0, 0
						t_rarp.append(t_0)
						ra_array.append(ra)
						rp_array.append(rp)
						print(t_0, rp, ra, file=file_R)
					m = float(data[8])
					q = float(data[9])

					if lineNumber > 3 and abs(m - m_prev)>1e-5:
						# color = 'm'
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
						t_no_doubles.append(t_0)
						# de_abs_tidal.append(de_abs_tidal[-1]+float(data[19]))
						de_tidal.append(de_tidal[-1]+e[-1]-e[-2])
						if t_0==t_previous:
							print('hmm, that\'s bad...', index, lineNumber)
							# dedt_tidal.append(dedt_tidal[-1])
						# else:
							# dedt_tidal.append(float(data[19])/(t_0-t_previous))
							# if lineNumber>1950 and lineNumber<1970:
							# 	print(lineNumber, "de =", de_abs_tidal[-1], "dt =", t_0-t_previous)
							# if len(dedt_tidal)>=2 and dedt_tidal[-1]>1e-4 and dedt_tidal[-2]<1e-4:
							# 	print("something happened here:", lineNumber)
					if lineNumber%3==0+shift and lineNumber>3:
						# de_abs_flybys.append(de_abs_flybys[-1]+abs(e[-1]-e[-2]))
						de_flybys.append(de_flybys[-1]+e[-1]-e[-2])
						if Q/a_0 < 25:
							t_strong_flybys.append(t_0)
							# de_abs_strong_flybys.append(de_abs_strong_flybys[-1]+abs(e[-1]-e[-2]))
				elif data[1] == 'calculation':	#N-body calculation abandoned
					t_prev = float(data[0])
					de_flybys.append(de_flybys[-1])
				elif data[1] == 'destroyed': 
					result = 'Binary destroyed (ionized)'
				elif data[1] == 'merger': 
					if len(exchange)>0:
						result = 'Binary merged after an exchange'
					else:
						result = 'Binary merged'
				elif data[1] == 'maximum' and data[2] == 'semimajor': 
					result = 'Calculation adandoned (semimajor axis too large)'
				elif data[1] == 'ejected':
					result = 'Binary ejected from the cluster'
					t_rarp.pop()
					rp_array.pop()
					ra_array.pop()
					rp_array[-1] = max(ra_array)
					ra_array[-1] = max(ra_array)

	e = np.array(e)
	de_flybys = np.array(de_flybys)
	de_tidal = np.array(de_tidal)
	e_flybys = de_flybys + e[0]
	e_tidal = de_tidal + e[0]

	figure = pyplot.figure(figsize=(6, 9)) 				
	figure.suptitle(fr'\noindent $m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a_0$ = {a_initial:.1f} AU, $e$ = {e_initial:.1f}, \\ $i_0 = {i_initial:.1f}^\circ$, $\omega_0$ = ${omega_initial:.1f}^\circ$, $\Omega_0$ = {Omega_initial:.0f} \\ $\chi={chi[i]}$, $\xi={D1D2[i]}$ \\ {result}', fontsize=16)

	gs = figure.add_gridspec(4, 1, hspace=0, wspace=0)
	r_plot, a_plot, e_plot, de_plot = gs.subplots(sharex=True)

	r_plot.minorticks_on() 
	r_plot.tick_params(labelsize=14)
	r_plot.set_ylabel(r'$R_a$, $R_p$ [pc]', fontsize=16)
	r_plot.plot(t_rarp, rp_array, 'r', label=r'$R_p$')
	r_plot.plot(t_rarp, ra_array, color, label=r'$R_a$')
	r_plot.legend(fontsize=16, frameon=False)

	a_plot.minorticks_on() 
	a_plot.tick_params(labelsize=14)
	a_plot.set_ylabel(r'$a$, $r_p$ [AU]', fontsize=16)
	a_plot.set_yscale('log')
	a_plot.plot(t, r_p, 'r', label=r'$r_p$')
	a_plot.plot(t, a, color, label=r'$a$')
	a_plot.legend(fontsize=16, frameon=False)

	start_index = 1
	end_index = 1
	if len(t_no_doubles) != len(de_flybys):
		t_no_doubles_new = t_no_doubles[:-1]
	else:
		t_no_doubles_new = t_no_doubles

	e_plot.minorticks_on() 
	e_plot.tick_params(labelsize=14)
	e_plot.set_xlabel(r'$t$ [Gyr]', fontsize=16)
	e_plot.set_ylabel(r'$1-e$', fontsize=16)
	e_plot.set_yscale('log')
	e_plot.plot(t, 1-e, 'k', label='actual eccentricity') 
	e_plot.plot(t_no_doubles_new, 1-e_flybys, 'r--', label='flybys only')
	e_plot.plot(t_no_doubles_new, 1-e_tidal, 'b:', label='tidal effects only')
	e_plot.legend(fontsize=16, frameon=False)

	de_plot.minorticks_on() 
	de_plot.tick_params(labelsize=14)
	de_plot.set_xlabel(r'$t$ [Gyr]', fontsize=16)
	de_plot.set_ylabel(r'$|\Delta e|$', fontsize=16)
	de_plot.set_yscale('log')
	de_plot.set_ylim(1e-3, 1)

	current_sign = np.sign(de_flybys[1])
	while end_index<len(de_flybys):
		while current_sign == np.sign(de_flybys[end_index]):
			end_index += 1
			if end_index==len(de_flybys):
				break
		if current_sign == 1:
			de_plot.plot(t_no_doubles_new[start_index:end_index], abs(de_flybys[start_index:end_index]), 'r', label=r'total $\Delta e$ due to flybys')
		else: 
			de_plot.plot(t_no_doubles_new[start_index:end_index], abs(de_flybys[start_index:end_index]), 'r:', label=r'total $\Delta e$ due to flybys, $\Delta e<0$')
		if end_index<len(de_flybys):
			start_index = end_index
			current_sign = np.sign(de_flybys[end_index])

	start_index = 1
	end_index = 1
	current_sign = np.sign(de_tidal[1])
	while end_index<len(de_tidal):
		while current_sign == np.sign(de_tidal[end_index]):
			end_index += 1
			if end_index==len(de_tidal):
				break
		if current_sign == 1:
			de_plot.plot(t_no_doubles[start_index:end_index], de_tidal[start_index:end_index], 'k', label=r'total $\Delta e$ due to tidal effects')
		else: 
			de_plot.plot(t_no_doubles[start_index:end_index], -de_tidal[start_index:end_index], 'k:', label=r'total $\Delta e$ due to tidal effects, $\Delta e<0$')	
		if end_index<len(de_tidal):
			start_index = end_index
			current_sign = np.sign(de_tidal[end_index])

	handles, labels = pyplot.gca().get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	pyplot.legend(by_label.values(), by_label.keys())
	pyplot.yscale('log')

	pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
	pyplot.savefig("output/for the paper 2/"+fileNames[i]+".pdf")
	pyplot.clf()