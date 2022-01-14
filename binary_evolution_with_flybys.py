import os
import numpy as np
import time
import astropy.units as u
from galpy.potential import evaluatePotentials, KeplerPotential, HernquistPotential, PlummerPotential
from galpy.util import conversion
from binary_evolution import KeplerRing, PointMass
from binary_evolution.tools import ecc_to_vel
from flybys3body import scattering_hybrid, scattering_SA
from amuse.lab import *
from numpy.random import default_rng

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

G = constants.G
c = constants.c
HubbleTime = 1.4e10|units.yr
_pc = 8000
_kms = 220

class inputParameters:
	def __init__(self, t=1e4, a_out=0.5, e_out=0, inc_out=np.pi/6, m1=5, m2=5, a=1, e=0.05, i=1, Omega=1.5, omega=0, output_file='output.txt', output_file_2='output2.txt', forcePrecise=False, potential="Plummer", rtol=1e-11, tmax=1e20, resume=False, includeWeakEncounters=True):
		self.t = t # Integration time [yr] 
		self.a_out = a_out # Outer orbit semi-major axis [pc]
		self.e_out = e_out # Outer orbit eccentricity
		self.inc_out = inc_out # Outer orbit inclination
		self.m1 = m1 # Primary mass [MSun]
		self.m2 = m2 # Secondatu mass [MSun]
		self.a = a # Inner orbit semimajor axis [AU]
		self.e = e # Inner orbit eccentricity
		self.i = i # Inner orbit inclination
		self.Omega = Omega # Inner orbit longitude os ascending node
		self.omega = omega # Inner orbit argument of periapsis
		self.output_file = output_file # Output file name
		self.output_file_2 = output_file_2 # Output file name
		self.forcePrecise = forcePrecise
		self.potential = potential
		self.rtol = rtol
		self.tmax = tmax
		self.resume = resume
		self.includeWeakEncounters = includeWeakEncounters

def m_final(m):
	stellar = SSE()
	stellar.parameters.metallicity = 0.001
	stellar.particles.add_particle(Particle(mass=m))
	stellar.evolve_model(5|units.Gyr)	
	result = stellar.particles.mass.in_(units.MSun)
	stellar.stop()
	return result

def sample_v_icdf (x, sigma_rel, n=10):
# x = GM/Q_max/sigma_rel^2
# u = v/sigma_rel
	rng = default_rng()
	random_number = rng.random()
	cdf1 = 0
	i = -1
	while cdf1 < random_number:
		i += 1
		u1 = (i+1)/n
		cdf1 = 1-np.exp(-u1**2/2) - u1**2*np.exp(-u1**2/2)/2/(1+x)
	u0 = i/n
	cdf0 = 1-np.exp(-u0**2/2) - u0**2*np.exp(-u0**2/2)/2/(1+x)	
	u = u0 + (u1-u0)*(random_number-cdf0)/(cdf1-cdf0)
	return u*sigma_rel

def sample_v_hamers (sigma_rel, v0):
# v0 = sqrt(GM/Q_max)
	rng = default_rng()
	while True:
		vz = np.sqrt(rng.exponential(2*sigma_rel.value_in(units.kms)**2))|units.kms
		vx = rng.normal(0, sigma_rel.value_in(units.kms))|units.kms
		vy = rng.normal(0, sigma_rel.value_in(units.kms))|units.kms
		v2 = vx**2+vy**2+vz**2
		if (v2 > 2*v0**2): break
	return np.sqrt(v2 - 2*v0**2)

m_total = 4e+6
b=1
# pot = TwoPowerTriaxialPotential(amp=16*m_bh*u.solMass, a=4*u.pc, alpha=1, beta=4, c=0.7)
pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 

Q_max_a_default = 50
Q_hybrid_a = 10
m_per = 1|units.MSun
m_per_max = m_per
# sigma_rel = 3|units.kms
# n = 1e6|units.pc**-3

def sigma_rel (r, type="Plummer"):
	# sqrt(2) * one-dimensional velocity dispersion
	if type=="Plummer": return np.sqrt(2)*np.sqrt(G*(m_total|units.MSun)/6/np.sqrt(r**2+(b|units.pc)**2))
	elif type=="Hernquist": 
		x = r/(b|units.pc)
		return np.sqrt(2)*G*(m_total|units.MSun)/12/(b|units.pc) * (12*x*(x+1)**3*log(1+x) - x/(x+1)*(25+52*x+42*x**2+12*x**3))
	else: return 0|units.kms

def rho (r, type="Plummer"):
	if type=="Plummer": return 3*(m_total|units.MSun)/4/np.pi/(b|units.pc)**3*(1+(r/(b|units.pc))**2)**-2.5
	elif type=="Hernquist": return (m_total|units.MSun)/2/np.pi/(b|units.pc)**3/(r/(b|units.pc))*(1+r/(b|units.pc))**-3
	else: return 0|units.kg/units.m**3

def tau_0 (a, m_bin, r, Q_max_a=50):
	Q_max = Q_max_a * a
	v0 = np.sqrt(G*(m_bin+m_per)/Q_max)
	return (2*np.sqrt(2*np.pi)*Q_max**2*sigma_rel(r)*(rho(r)/m_per)*(1+(v0/sigma_rel(r))**2))**-1

def a_h (m1, m2, r):
	return (G*(m1*m2/(m1+m2)|units.MSun)/4/sigma_rel(r|units.pc)**2).value_in(units.AU)

def sample_encounter_parameters (a, m_bin, r, Q_max_a=50):
	Q_max = Q_max_a * a
	rng = default_rng()
	v0 = np.sqrt(G*(m_bin+m_per)/Q_max)
	# time until the encounter
	# tau = rng.exponential(tau_0(a,m_bin,n).value_in(units.yr))|units.yr
	# perturber mass
	# m_per = 1|units.MSun
	# perturber orbital parameters
	# v = sample_v_hamers (sigma_rel(r), v0)
	x = (v0/sigma_rel(r))**2
	v = sample_v_icdf (x, sigma_rel(r))
	p_max2 = Q_max**2*(1+2*(v0/v)**2)
	p = np.sqrt(p_max2*rng.random())
	aStar = -G*(m_bin+m_per)/v**2
	eStar = np.sqrt(1+p**2/aStar**2)
	iStar = np.arccos(rng.random())
	OmegaStar = rng.random()*2*np.pi
	omegaStar = rng.random()*2*np.pi
	return m_per, aStar, eStar, iStar, OmegaStar, omegaStar

# Inner binary parameters
a_in = 0.01              # Semi-major axis in AU
ecc = 0.05            	# Eccentricity
inc = np.pi/3           # Inclination with respect to the z-axis
long_asc = 0            # Longitude of the ascending node
arg_peri = np.pi / 2    # Arugment of pericentre
m1 = 5
m2 = 5
m_bin = m1+m2                  # Total mass in solar masses
m_bin_init = m_bin|units.MSun
# q = 1

M_max = max(m_bin_init+m_per_max, 3*m_per_max) 

# Outer binary parameters
ecc_out = 0.0         # Outer orbit eccentricity
inc_out = np.pi/6             # Outer orbit inclination
a_out = 0.5        # Outer semi-major axis in pc

# Start at pericentre
r = a_out * (1 - ecc_out)   # Spherical radius
R = r * np.cos(inc_out)     # Cylindrical radius
z = r * np.sin(inc_out)     # Cylindrical height

# Compute the correct v_phi at pericentre for the selected eccentricity and potential
v_phi = ecc_to_vel(pot, ecc_out, [R, z, 0])

# Define the KeplerRing
k = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=1)

# output_file_name = os.path.dirname(os.path.abspath(__file__))+'/history-rtol7.txt'

# print(np.sum([[[0, 0 ,0], [2, 2, 2]], [[1, 2, 3], 2]], axis=0))
# def func1 (x, y):
# 	return x, y
# def func2 (x, y):
# 	return x
# funcs=[]
# funcs.append(func1)
# funcs.append(lambda x, y: (func2(x, y), 0))
# def derivatives(x, y):
#     return np.sum([f(x, y) for f in funcs], axis=0)
# x = np.ndarray([0, 0, 0])
# y = np.ndarray([1, 1, 1])
# print(derivatives(x, y))

# k.integrate(ts, pot=pot, relativity=False, gw=False)
# print(k.r(ts[-1]), k.v(ts[-1]))
# k = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=1)
# k.integrate(ts, pot=pot, relativity=False, gw=False, rtol=1e-3, atol=1e-6)
# print(k.r(ts[-1]), k.v(ts[-1]))

# test_file_name = 'test.txt'
# test_file = open(test_file_name, 'w+')
# t_max = 1e7
# ts = np.linspace(0, t_max, 1000)
# k.integrate(ts, pot=pot, relativity=True, gw=True, tau_0=lambda *args: tau_0(args[0]|units.pc, m_bin|units.MSun, args[1]|units.pc).value_in(units.yr), random_number=100)#, rtol=1e-7, atol=1e-10)
# t1 = k.t_fin
# print(k.r(t1))
# k.integrate(ts, pot=pot, relativity=True, gw=True, tau_0=lambda *args: tau_0(args[0]|units.pc, m_bin|units.MSun, args[1]|units.pc).value_in(units.yr), random_number=1e10)#, rtol=1e-7, atol=1e-10)
# print(k.r(t1))

# for t1 in ts:
# 	print(t1, k.r(t1), k.v(t1), file=test_file)

# is the binary hard?
# a_h = G*(m1*m2/(m1+m2)|units.MSun)/4/sigma_rel(r|units.pc)**2
# # print(a_in/a_h.value_in(units.AU))
# Q=0.25
# print(((a_in|units.AU)/(64/5 * Q * G**3 * (k.m()|units.MSun)**3 / c**5 / (a_in|units.AU)**3)).value_in(units.yr))
# print(tau_0 (a_in|units.AU, m_bin|units.MSun, r|units.pc).value_in(units.yr))

# t1 = 481719441.126
# R1 = 0.146515758464
# z1 = -0.199396101429 
# phi1 = 1.35525468694
# v1 = [-71.0051074216, -0.220052592881, 3.17866255331]
# t2 = 481721453.32
# R2 = 0.00708097862563
# z2 = -0.193131002586 
# phi2 = -2.93773350624
# v2 = [19.1350160868, 6.31164318845, 49.9642453788]
# random_number = 1.5219247855019409
# dt = 4187.69708645
# n = 20

# k = KeplerRing(ecc, inc, long_asc, arg_peri, [R1, z1, phi1], v1, a=a_in, m=m_bin, q=1)
# ts = np.linspace(0, 2*(t2-t1), 1000)
# k.integrate(ts, pot=pot, relativity=True, gw=True, tau_0=lambda *args: tau_0(args[0]|units.pc, m_bin|units.MSun, args[1]|units.pc).value_in(units.yr), random_number=1.4, rtol=1e-7, atol=1e-10)
# R2real, z2real, phi2real = k.r(t2-t1)
# v2real = k.v(t2-t1)

# E1 = (np.linalg.norm(v1)/_kms)**2/2 + evaluatePotentials(pot, R1/_pc, z1/_pc, phi=phi1, use_physical=False) 
# E2 = (np.linalg.norm(v2)/_kms)**2/2 + evaluatePotentials(pot, R2/_pc, z2/_pc, phi=phi2, use_physical=False) 
# E2real = (np.linalg.norm(v2real)/_kms)**2/2 + evaluatePotentials(pot, R2real/_pc, z2real/_pc, phi=phi2real, use_physical=False) 

def approximation_test (input):

	t = input.t
	if input.includeWeakEncounters: Q_max_a = Q_max_a_default
	else: Q_max_a = Q_hybrid_a

	# Outer binary parameters
	a_out = input.a_out        # Outer semi-major axis in pc
	ecc_out = input.e_out         # Outer orbit eccentricity
	inc_out = input.inc_out        # Outer orbit inclination

	# Inner binary parameters
	m1 = input.m1
	m2 = input.m2
	a_in = input.a              # Semi-major axis in AU
	ecc = input.e           	# Eccentricity
	inc = input.i           # Inclination with respect to the z-axis
	arg_peri = input.omega     # Arugment of pericentre
	long_asc = input.Omega             # Longitude of the ascending node
	m_bin = m1+m2                  # Total mass in solar masses

	# Start at pericentre
	r = a_out * (1 - ecc_out)   # Spherical radius
	R = r * np.cos(inc_out)     # Cylindrical radius
	z = r * np.sin(inc_out)     # Cylindrical height

	if input.potential=="Plummer": pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 
	elif input.potential=="Hernquist": pot = HernquistPotential(amp=m_total*u.solMass, a=b*u.pc) 

	# Compute the correct v_phi at pericentre for the selected eccentricity and potential
	v_phi = ecc_to_vel(pot, ecc_out, [R, z, 0])

	# Define the KeplerRing
	k = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=m2/m1)
	k1 = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=m2/m1)

	output_file = open(input.output_file, 'w+')
	print('t[yr] R[pc] z phi v_R[km/s] v_z v_phi a[AU] m[MSun] q ecc inc long_asc arg_peri', file=output_file)
	print(0, R, z, 0, 0, 0, v_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), file=output_file, flush=True)

	ts = np.linspace(0, t, 1000)
	rtol=input.rtol #1e-11
	atol= rtol*1e-3 #1e-14

	k.integrate(ts, pot=pot, relativity=True, gw=True, tau_0=lambda *args: tau_0(args[0]|units.pc, k.m()|units.MSun, args[1]|units.pc, Q_max_a=Q_max_a).value_in(units.yr), random_number=1e10, rtol=rtol, atol=atol, forcePrecise=False)
	print('gr_ratio =', k.gr_ratio, ', t =', t, file=output_file)
	print('da de di dOmega domega', file=output_file)
	if k.merger: print('approximate: merger at', k.t_fin, file=output_file, flush=True)
	else: print('approximate: ', k.a_fin-a_in, k.ecc_fin-ecc, k.inc_fin-inc, k.long_asc_fin-long_asc, k.arg_peri_fin-arg_peri, file=output_file, flush=True)

	k1.integrate(ts, pot=pot, relativity=True, gw=True, tau_0=lambda *args: tau_0(args[0]|units.pc, k.m()|units.MSun, args[1]|units.pc, Q_max_a=Q_max_a).value_in(units.yr), random_number=1e10, rtol=rtol, atol=atol, forcePrecise=True)
	if k1.merger: print('precise: merger at', k1.t_fin, file=output_file, flush=True)
	else: print('precise: ', k1.a_fin-a_in, k1.ecc_fin-ecc, k1.inc_fin-inc, k1.long_asc_fin-long_asc, k1.arg_peri_fin-arg_peri, file=output_file, flush=True)

	matplotlib.rcParams['text.usetex'] = True
	matplotlib.rcParams['mathtext.fontset'] = 'stix'
	matplotlib.rcParams['font.family'] = 'STIXGeneral'
	matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
	figure = pyplot.figure()

	plot_a = figure.add_subplot(2,1,1)
	ax = pyplot.gca()
	ax.minorticks_on() 
	# ax.xaxis.set_major_locator(MultipleLocator(1))
	# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
	# ax.yaxis.set_major_locator(MultipleLocator(0.1))
	# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
	ax.tick_params(labelsize=14)
	ax.set_xlabel(r'$t$ [yr]', fontsize=16)
	ax.set_ylabel(r'$\Delta a$ [AU]', fontsize=16)
	# pyplot.xscale('log')
	# pyplot.yscale('log')
	# pyplot.text(0.3, 0.9, '$e='+str(e)+'$', fontsize=16, transform=ax.transAxes)
	plot_a.plot(ts, k.a_array-a_in, 'k', label=r'GR only')
	plot_a.plot(ts, k1.a_array-a_in, 'r', label=r'GR + tidal')
	plot_a.legend(fontsize=16, frameon=False)
	# plot.set_xlim(500, 5000)
	# plot.set_ylim(1e-8, 1e-5)

	plot_e = figure.add_subplot(2,1,2)
	ax = pyplot.gca()
	ax.minorticks_on() 
	# ax.xaxis.set_major_locator(MultipleLocator(1))
	# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
	# ax.yaxis.set_major_locator(MultipleLocator(0.1))
	# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
	ax.tick_params(labelsize=14)
	ax.set_xlabel(r'$t$ [yr]', fontsize=16)
	ax.set_ylabel(r'$\Delta e$', fontsize=16)
	# pyplot.xscale('log')
	# pyplot.yscale('log')
	# pyplot.text(0.3, 0.9, '$e='+str(e)+'$', fontsize=16, transform=ax.transAxes)
	plot_e.plot(ts, k.e_array-ecc, 'k', label=r'GR only')
	plot_e.plot(ts, k1.e_array-ecc, 'r', label=r'GR + tidal')
	plot_e.legend(fontsize=16, frameon=False)
	# plot.set_xlim(500, 5000)
	# plot.set_ylim(1e-8, 1e-5)

	pyplot.tight_layout()
	pyplot.savefig(input.output_file_2)

def evolve_binary (input):
	# 0 - binary has survived until t_final
	# 1 - binary has merged
	# 2 - binary has been destroyed
	# 3 - maximum calculation time exceeded

	t_final = input.t|units.yr
	
	if input.includeWeakEncounters: Q_max_a = Q_max_a_default
	else: Q_max_a = Q_hybrid_a
	
	if input.potential=="Plummer": pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 
	elif input.potential=="Hernquist": pot = HernquistPotential(amp=m_total*u.solMass, a=b*u.pc) 

	if input.resume and os.path.isfile(input.output_file):
		with open(input.output_file) as f:
			for line in f: pass
			data = line.split()
			if data[-1] == 'merger': return 1
			if data[-1] == 'destroyed': return 2
			t = float(data[0])|units.yr
			if t>t_final: return 0
			R = float(data[1])
			z = float(data[2])
			phi = float(data[3])
			v_R = float(data[4])
			v_z = float(data[5])
			v_phi = float(data[6])
			a_in = float(data[7])
			m_bin = float(data[8])
			q = float(data[9])
			ecc = float(data[10])
			inc = float(data[11])
			long_asc = float(data[12])
			arg_peri = float(data[13])
		k = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, phi], [v_R, v_z, v_phi], a=a_in, m=m_bin, q=q)
		output_file = open(input.output_file, 'a')
	else:
		t = 0|units.yr

		# Outer binary parameters
		a_out = input.a_out        # Outer semi-major axis in pc
		ecc_out = input.e_out         # Outer orbit eccentricity
		inc_out = input.inc_out        # Outer orbit inclination,

		# Inner binary parameters
		m1 = input.m1
		m2 = input.m2
		a_in = input.a              # Semi-major axis in AU
		ecc = input.e           	# Eccentricity
		inc = input.i           # Inclination with respect to the z-axis
		arg_peri = input.omega     # Arugment of pericentre
		long_asc = input.Omega             # Longitude of the ascending node
		m_bin = m1+m2                  # Total mass in solar masses
		m_bin_init = m_bin|units.MSun

		# Start at pericentre
		r = a_out * (1 - ecc_out)   # Spherical radius
		R = r * np.cos(inc_out)     # Cylindrical radius
		z = r * np.sin(inc_out)     # Cylindrical height

		# Compute the correct v_phi at pericentre for the selected eccentricity and potential
		v_phi = ecc_to_vel(pot, ecc_out, [R, z, 0])

		# Define the KeplerRing
		k = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=m2/m1)

		output_file = open(input.output_file, 'w+')
		print('t[yr] R[pc] z phi v_R[km/s] v_z v_phi a[AU] m[MSun] q ecc inc long_asc arg_peri random_number_0 dt[yr] n gr_ratio', file=output_file)
		print('perturber: m_per[MSun] Q[AU] eStar iStar OmegaStar omegaStar', file=output_file)
		print(0, R, z, 0, 0, 0, v_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), file=output_file)
		output_file.flush()

	rtol=input.rtol #1e-11
	atol= rtol*1e-3 #1e-14

	timeTotal1 = time.time()
	timeClose = 0
	timeDistant = 0
	timeOrbit = 0
	timeLoop = 0
	while t<t_final:
		# bad_orbit_file = open('bad_orbit.txt', 'w+')
		# print(t.value_in(units.yr), k.r(), k.v(), k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), file=output_file)
		# output_file = open(output_file_name, 'w+')
		# print("heyyyyyy", file=output_file)
		# output_file.flush()
		# integrate the orbit until the next flyby
		rng = default_rng()
		random_number = rng.exponential()
		random_number_0 = random_number
		r = np.sqrt(k.r()[0]**2+k.r()[1]**2)
		tau_0_value = tau_0 (k.a()|units.AU, k.m()|units.MSun, r|units.pc, Q_max_a=Q_max_a)
		T = 2*np.pi*(r|units.pc)/sigma_rel(r|units.pc)	# approximate outer period
		# integration_time = tau_0_value*random_number
		timeOrbit1 = time.time()
		# if integration_time < 0.01*T:
		# 	# print("simple")
		# 	t += integration_time
		# 	if t > HubbleTime: reachedHubbleTime = True
		# 	ts = np.linspace(0, integration_time.value_in(units.yr), 10)
		# 	k.integrate(ts, pot=pot, relativity=False, gw=False)
		# 	# assign new orbital parameters to the binary
		# 	k = KeplerRing(k.ecc(ts[-1]), k.inc(ts[-1]), k.long_asc(ts[-1]), k.arg_peri(ts[-1]), k.r(ts[-1]), k.v(ts[-1]), a=k.a(), m=k._m, q=k._q)
		# else:
		# print("complicated")
		# integral = 0
		# dt = min(5*tau_0_value*random_number, 1e6|units.yr)
		Q = k._q / (1+k._q)**2
		t_gw = (k.a()|units.AU)/(64/5 * Q * G**3 * (k.m()|units.MSun)**3 / c**5 / (k.a()|units.AU)**3)
		# print(t_gw.value_in(units.yr))
		dt = 1.1*min(tau_0_value*random_number, t_gw, (t_final-t))
		# print(dt.value_in(units.yr))
		n = max(int(dt/(0.01*T)), 10)
		# previous_tau_0_value = 0
		while (random_number>0):
			# print(dt.value_in(units.yr))
			ts = np.linspace(0, dt.value_in(units.yr), 100*n+1)
			# ts = np.linspace(0, 1e4, n+1)
			k.integrate(ts, pot=pot, relativity=True, gw=True, tau_0=lambda *args: tau_0(args[0]|units.pc, k.m()|units.MSun, args[1]|units.pc, Q_max_a=Q_max_a).value_in(units.yr), random_number=random_number, rtol=rtol, atol=atol, forcePrecise=input.forcePrecise, debug_file=input.output_file_2) #, rtol=1e-3, atol=1e-6)
			# timeLoop1 = time.time()
			# for i in range(n):
			# 	R, z, phi = k.r(ts[i])
			# 	r = np.sqrt(R**2+z**2)
			# 	tau_0_value = tau_0 (k.a()|units.AU, k.m()|units.MSun, r|units.pc)
			# 	if i>0 and abs(previous_tau_0_value/tau_0_value-1)>0.2: print("big error") 
			# 	previous_tau_0_value = tau_0_value
			# 	integral += dt/n/tau_0_value
			# 	t += dt/n
			# 	if integral >= random_number: break
			t += k.t_fin|units.yr
			# print(t.value_in(units.yr), file=output_file)
			# output_file.flush()
			if k.merger: break
			# if t > HubbleTime: 
			# 	reachedHubbleTime = True
			# 	break
			# timeLoop2 = time.time()
			# timeLoop += timeLoop2 - timeLoop1
			# assign new orbital parameters to the binary
			# print(k.a())
			# k = KeplerRing(k.ecc(ts[i+1]), k.inc(ts[i+1]), k.long_asc(ts[i+1]), k.arg_peri(ts[i+1]), k.r(ts[i+1]), k.v(ts[i+1]), a=k.a_array[i+1], m=k._m, q=k._q)
			random_number = k.probability
			gr_ratio = k.gr_ratio
			# print(random_number)
			k = KeplerRing(k.ecc_fin, k.inc_fin, k.long_asc_fin, k.arg_peri_fin, k.r(k.t_fin), k.v(k.t_fin), a=k.a_fin, m=k._m, q=k._q)
			if t>=t_final: break
		# t_list.add(t.value_in(units.yr))
		timeOrbit2 = time.time()
		timeOrbit += timeOrbit2 - timeOrbit1
		# print('t[yr] R[pc] z phi v_R[km/s] v_z v_phi a[AU] m[MSun] q ecc inc long_asc arg_peri')
		R, z, phi = k.r()
		v_R, v_z, v_phi = k.v()
		# print(t.value_in(units.yr), k.ecc())
		if k.merger:
			print(t.value_in(units.yr), "merger", file=output_file)
			return 1
		print(t.value_in(units.yr), R, z, phi, v_R, v_z, v_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), random_number_0, dt.value_in(units.yr), n, gr_ratio, file=output_file)
		output_file.flush()
		if t>=t_final: return 0

		# t = 1e4|units.yr
		# Q=0.25
		# print(ecc/((a_in|units.AU)**2.5 * c**2 * (1 - ecc**2)**1.5 / 3 / (G * (k.m()|units.MSun))**1.5) * t)
		# print((-64/5 * Q * G**3 * (k.m()|units.MSun)**3 / c**5 / (a_in|units.AU)**3 * t).value_in(units.AU) / a_in)
		# ts = np.linspace(0, t.value_in(units.yr), 10)
		# time1 = time.time()
		# k.integrate(ts, pot=pot, relativity=False, gw=True, rtol=1e-5, atol=1e-8)
		# time2 = time.time()
		# print("dt = ", time2-time1, " s")

		# print(k.long_asc(ts[-1]))
		# # assign new orbital parameters to the binary
		# k = KeplerRing(k.ecc(ts[-1]), k.inc(ts[-1]), k.long_asc(ts[-1]), k.arg_peri(ts[-1]), k.r(ts[-1]), k.v(ts[-1]), a=k.a(), m=k._m, q=k._q)
		# print(k.long_asc())
		# k.integrate(ts, pot=pot, relativity=True)
		# print(k.long_asc(ts[-1]))
		# # assign new orbital parameters to the binary
		# k = KeplerRing(k.ecc(ts[-1]), k.inc(ts[-1]), k.long_asc(ts[-1]), k.arg_peri(ts[-1]), k.r(ts[-1]), k.v(ts[-1]), a=k.a(), m=k._m, q=k._q)
		# print(k.long_asc())

		if 1==1:
			# sample the perturber parameters
			m_per, aStar, eStar, iStar, OmegaStar, omegaStar = sample_encounter_parameters (k.a()|units.AU, k.m()|units.MSun, np.sqrt(R**2+z**2)|units.pc, Q_max_a=Q_max_a)
			Q = aStar*(1-eStar)
			# print('perturber: m_per[MSun] aStar[AU] eStar iStar OmegaStar omegaStar', file=output_file)
			# print('perturber: ', m_per.value_in(units.MSun), aStar.value_in(units.AU), eStar, iStar, OmegaStar, omegaStar, file=output_file)
			print('perturber: ', m_per.value_in(units.MSun), Q.value_in(units.AU), eStar, iStar, OmegaStar, omegaStar, file=output_file)
			output_file.flush()

			# perform the scattering
			q = k._q
			m1 = k.m()/(1+q)
			m2 = k.m()*q/(1+q)
			if Q<=Q_hybrid_a*(k.a()|units.AU):
				timeClose1 = time.time()
				result, third_body_final, dv_binary, a_fin, e_fin, i_fin, Omega_fin, omega_fin, n_3body = scattering_hybrid (m1|units.MSun, m2|units.MSun, k.a()|units.AU, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), m_per, aStar, eStar, iStar, OmegaStar, omegaStar)
				timeClose2 = time.time()
				timeClose += timeClose2 - timeClose1 
			else:
				if 1==1:
					result = 0
					third_body_final = 2
					timeDistant1 = time.time()
					dv_binary, a_fin, e_fin, i_fin, Omega_fin, omega_fin = scattering_SA (m1|units.MSun, m2|units.MSun, k.a()|units.AU, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), m_per, aStar, eStar, iStar, OmegaStar, omegaStar)
					timeDistant2 = time.time()
					timeDistant += timeDistant2 - timeDistant1 
				else:
					result = 4	#ignore weak interactions

			if result == 2:
				print(t.value_in(units.yr), "destroyed", file=output_file)
				return 2 
			elif result == 1: 
				print(t.value_in(units.yr), "calculation abandoned after more than n_orbits_max bound orbits", file=output_file)
				output_file.flush()
			elif result == 3: 
				print(t.value_in(units.yr), "calculation abandoned after spending too much time in a 3-body phase", file=output_file)
				output_file.flush()
			elif result == 0:
				# assign new orbital parameters to the binary
				if third_body_final == 0: m1 = m_per.value_in(units.MSun)
				if third_body_final == 1: m2 = m_per.value_in(units.MSun)
				v_R, v_z, v_phi = k.v()	#velocity before the scattering
				x = R * np.cos(phi)
				y = R * np.sin(phi)
				dv_x = dv_binary[0].value_in(units.kms)	#velocity change during the scattering
				dv_y = dv_binary[1].value_in(units.kms)
				dv_z = dv_binary[2].value_in(units.kms)
				dv_R = (x*dv_x+y*dv_y)/R
				phi_unit_vector = np.cross([0, 0, 1], [x/R, y/R, 0])
				dv_phi = np.dot(phi_unit_vector, [dv_x, dv_y, dv_z])
				# print(e_fin, i_fin, Omega_fin, omega_fin, [R, z, phi], [v_R+dv_R, v_z+dv_z, v_phi+dv_phi], a_fin.value_in(units.AU), m1+m2)
				# print(a_fin.value_in(units.AU)-a_in)
				k = KeplerRing(e_fin, i_fin.value_in(units.rad), Omega_fin.value_in(units.rad), omega_fin.value_in(units.rad), [R, z, phi], [v_R+dv_R, v_z+dv_z, v_phi+dv_phi], a=a_fin.value_in(units.AU), m=m1+m2, q=min(m1/m2, m2/m1))
				m_bin = m1+m2
				print(t.value_in(units.yr), R, z, phi, v_R+dv_R, v_z+dv_z, v_phi+dv_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), file=output_file)
				output_file.flush()
		
		timeTotal = time.time()-timeTotal1
		if timeTotal>input.tmax: 
			print('Maximum calculation time (', str(input.tmax), ' s) exceeded', file=output_file)
			output_file.flush()
			return 3

	# test_file_name = 'test.txt'
	# test_file = open(test_file_name, 'w+')
	# t_max = 1e7
	# ts = np.linspace(0, t_max, 1000)
	# k.integrate(ts, pot=pot, relativity=True, gw=True, tau_0=lambda *args: tau_0(args[0]|units.pc, m_bin|units.MSun, args[1]|units.pc).value_in(units.yr), random_number=100)#, rtol=1e-7, atol=1e-10)
	# for t1 in t_list:

	timeTotal2 = time.time()
	print("total time", timeTotal2-timeTotal1, "s", file=output_file)
	print("close interaction time", timeClose, "s", file=output_file)
	print("distant interaction time", timeDistant, "s", file=output_file)
	print("outer orbit integration time", timeOrbit, "s", file=output_file)
	output_file.close()

	return 0

# R, z, phi = k.r(ts).T
# x = R * np.cos(phi)
# y = R * np.sin(phi)

# perturber parameters
# m_per = 1|units.MSun
# v = 3 | units.kms
# iStar = 0.
# OmegaStar = 0.1
# omegaStar = np.pi/2
# aStar = -constants.G*(m_per+(m_bin|units.MSun))/v**2
# Q = 20*a_in|units.AU
# eStar = 1 - Q/aStar

# print(np.sqrt(x**2+y**2))
# plt.plot(x, y)
# plt.xlabel('x (pc)')
# plt.ylabel('y (pc)')
# plt.show()

# e = k.ecc(ts)
# arg_peri = k.arg_peri(ts)
# long_asc = k.long_asc(ts)
# print(long_asc)
# i = k.inc(ts)

# plt.plot(ts, e)
# plt.xlabel('t (years)')
# plt.ylabel('e')
# plt.show()

# plt.plot(ts, i)
# plt.xlabel('t (years)')
# plt.ylabel('i')
# plt.show()
