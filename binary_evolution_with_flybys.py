import os
import numpy as np
import time
import astropy.units as u
from galpy.potential import evaluatePotentials, KeplerPotential, HernquistPotential, PlummerPotential
from galpy.util import conversion
from binary_evolution import KeplerRing, PointMass
from binary_evolution.tools import ecc_to_vel
# from flybys3body import scattering_hybrid, scattering_SA
from fortran.flybys3body_fortran import scattering_hybrid, scattering_SA
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

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

class inputParameters:
	def __init__(self, t=1e4, a_out=0.5, e_out=0, inc_out=np.pi/6, m1=5, m2=5, a=1, e=0.05, i=1, Omega=1.5, omega=0, output_file='output.txt', output_file_2='', approximation=0, potential="Plummer", m_total=4e6, b=1, rtol=1e-11, tmax=1e20, relativity=True, tidal_effects=True, gw=True, resume=False, includeEncounters=True, includeWeakEncounters=True, Q_max_a=50, Q_min_a=0, n=10, a_max=1000, sameParameters='', disableKicks=False, t0=0, m_per=1):
		self.t = t # Integration time [yr] 
		self.a_out = a_out # Outer orbit semi-major axis [pc]
		self.e_out = e_out # Outer orbit eccentricity
		self.inc_out = inc_out # Outer orbit inclination
		self.m1 = m1 # Primary mass [MSun]
		self.m2 = m2 # Secondatu mass [MSun]
		self.a = a # Inner orbit semimajor axis [AU]
		self.e = e # Inner orbit eccentricity
		self.i = i # Inner orbit inclination
		self.Omega = Omega # Inner orbit longitude of ascending node
		self.omega = omega # Inner orbit argument of periapsis
		self.output_file = output_file # Output file name
		self.output_file_2 = output_file_2 # Additional output file name
		self.approximation = approximation # 0 - use precise if epsilon_gr<20 and GR-only otherwise; 1 - always precise; 2 - always GR-only
		self.potential = potential # Cluster potential
		self.b = b
		self.m_total = m_total
		self.rtol = rtol # Inner orbit integration accuracy
		self.tmax = tmax # Maximum calculation time [s]
		self.resume = resume # Resume the integration from the last line in self.output_file (ignores the provided initial conditions)
		self.includeWeakEncounters = includeWeakEncounters  
		self.Q_max_a = Q_max_a # The maximum pericenter of the encounters to include
		self.Q_min_a = Q_min_a # The minimum pericenter of the encounters to include
		self.includeEncounters = includeEncounters 
		self.n = n # The number of points per (approximate) outer orbital period used to interpolate the outer orbit 
		self.relativity = relativity #include GR effects
		self.tidal_effects = tidal_effects #include tidal terms
		self.gw = gw #include GW emission
		self.a_max = a_max #Stop the integration if the inner binary semimajor axis exceeds this value in AU
		self.sameParameters = sameParameters #if not empty, take the initial conditions from that file (overwritten by resume)
		self.disableKicks = disableKicks #if True, encounters don't change the binary CM velocity 
		self.t0 = t0 	#initial time
		self.m_per = m_per 	#Perturber mass in M_Sun

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

# m_total = 4e+6
# b=1
# pot = TwoPowerTriaxialPotential(amp=16*m_bh*u.solMass, a=4*u.pc, alpha=1, beta=4, c=0.7)
# pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 

# Q_max_a_default = 50
Q_hybrid_a = 10
# m_per = 1|units.MSun
# m_per_max = m_per
# sigma = 3|units.kms
# n = 1e6|units.pc**-3

def sigma (r, type="Plummer", m_total=4e6, b=1):
	# sqrt(2) * one-dimensional velocity dispersion
	if type=="Plummer": return np.sqrt(G*(m_total|units.MSun)/6/np.sqrt(r**2+(b|units.pc)**2))
	elif type=="Hernquist": 
		x = r/(b|units.pc)
		# print(x, 12*x*(x+1)**3*np.log(1+1/x) - x/(x+1)*(25+52*x+42*x**2+12*x**3))
		return np.sqrt(G*(m_total|units.MSun)/12/(b|units.pc) * (12*x*(x+1)**3*np.log(1+1/x) - x/(x+1)*(25+52*x+42*x**2+12*x**3)))
	else: 
		return 0|units.kms

def rho (r, type="Plummer", m_total=4e6, b=1):
	if type=="Plummer": return 3*(m_total|units.MSun)/4/np.pi/(b|units.pc)**3*(1+(r/(b|units.pc))**2)**-2.5
	elif type=="Hernquist": return (m_total|units.MSun)/2/np.pi/(b|units.pc)**3/(r/(b|units.pc))*(1+r/(b|units.pc))**-3
	else: return 0|units.kg/units.m**3

def tau_0 (a, m_bin, r, Q_max_a=50, type="Plummer", m_total=4e6, b=1, V=0|units.kms, m_per=1|units.MSun):
	Q_max = Q_max_a * a
	v0 = np.sqrt(G*(m_bin+m_per)/Q_max)
	sigma_rel = np.sqrt(sigma(r, type, m_total, b)**2 + V**2)
	return (2*np.sqrt(2*np.pi)*Q_max**2*sigma_rel*(rho(r, type, m_total, b)/m_per)*(1+(v0/sigma_rel)**2))**-1

def a_h (m1, m2, r, type="Plummer", m_total=4e6, b=1):
	return (G*(m1*m2/(m1+m2)|units.MSun)/4/sigma(r|units.pc, type, m_total, b)**2).value_in(units.AU)

def sample_encounter_parameters_old (a, m_bin, r, Q_max_a=50, type="Plummer", m_total=4e6, b=1, m_per=1|units.MSun):
	Q_max = Q_max_a * a
	rng = default_rng()
	v0 = np.sqrt(G*(m_bin+m_per)/Q_max)
	# time until the encounter
	# tau = rng.exponential(tau_0(a,m_bin,n).value_in(units.yr))|units.yr
	# perturber mass
	# m_per = 1|units.MSun
	# perturber orbital parameters
	# v = sample_v_hamers (sigma(r), v0)

	# sigma -> sigma_rel!
	x = (v0/sigma(r, type, m_total, b))**2
	v = sample_v_icdf (x, sigma(r, type, m_total, b))
	p_max2 = Q_max**2*(1+2*(v0/v)**2)
	p = np.sqrt(p_max2*rng.random())
	aStar = -G*(m_bin+m_per)/v**2
	eStar = np.sqrt(1+p**2/aStar**2)
	iStar = np.arccos(rng.random()*2-1)
	OmegaStar = rng.random()*2*np.pi
	omegaStar = rng.random()*2*np.pi
	return m_per, aStar, eStar, iStar, OmegaStar, omegaStar

def normalize (vector):
	return vector / np.linalg.norm(vector)

def sample_encounter_parameters (a, m_bin, r, phi, Q_max_a=50, type="Plummer", m_total=4e6, b=1, v_bin=[0,0,0], m_per=1|units.MSun):
	# returns the parameters in physical units
	# v_bin is the binary CM velocity [v_R, v_phi, v_z]=[v_x, v_y, v_z]  in km/s 
	rng = default_rng()

	Q_max = Q_max_a * a.value_in(units.m)
	mu = (G*(m_bin+m_per)).value_in(units.m**3/units.s**2)
	v0 = np.sqrt(mu/Q_max)
	x = (v0/sigma(r, type, m_total, b).value_in(units.m/units.s))**2
	v = sample_v_icdf (x, sigma(r, type, m_total, b)).value_in(units.m/units.s)

	theta_v = np.arccos(rng.random()*2-1)
	phi_v = rng.random()*2*np.pi
	v_x = v*np.sin(theta_v)*np.cos(phi_v)
	v_y = v*np.sin(theta_v)*np.sin(phi_v)
	v_z = v*np.cos(theta_v)
	#initial velocity vector in the binary reference frame
	v1 = [v_x-v_bin[0]*1000, v_y-v_bin[1]*1000, v_z-v_bin[2]*1000]

	p_max2 = Q_max**2*(1+2*(v0/np.linalg.norm(v1))**2)
	p = np.sqrt(p_max2*rng.random())
	aStar = -mu/np.linalg.norm(v1)**2
	eStar = np.sqrt(1+p**2/aStar**2)

	# unit vectors perpendicular to velocity
	e1 = normalize([v1[1], -v1[0], 0])
	e2 = normalize(np.cross(v1, e1))
	angle = rng.random()*2*np.pi
	# angular momentum
	h = (e1*np.cos(angle) + e2*np.sin(angle)) * np.sqrt(mu*aStar*(1-eStar**2))
	# eccentricity vector
	ecc_vector = np.cross(v1, h)/mu + normalize(v1)

	# now we convert all the vectors to the cluster reference frame
	# the unut vectors for that reference frame in cylindrical coordinates
	cluster_x = [np.cos(phi), -np.sin(phi), 0]
	cluster_y = [np.sin(phi), np.cos(phi), 0]
	cluster_z = [0, 0, 1]
	rotation_matrix = [cluster_x, cluster_y, cluster_z]
	# transforming the vectors to the cluster reference frame
	h = np.matmul(rotation_matrix, h)
	ecc_vector = np.matmul(rotation_matrix, ecc_vector)

	# calculating the orbital angles
	iStar = np.arccos(normalize(h)[2])
	# vector pointitng towards the ascending node
	n = normalize([-h[1], h[0], 0])
	if n[0]==0 and n[1]==0:
		OmegaStar = 0
	else:
		OmegaStar = np.arccos(n[0])
		if n[1]<0: OmegaStar = 2*np.pi - OmegaStar
	if h[0]==0 and h[1]==0:
		omegaStar = np.arctan2(ecc_vector[1], ecc_vector[0])
		if h[2]<0: omegaStar = 2*np.pi - omegaStar
	else:
		omegaStar = np.arccos(np.dot(n, normalize(ecc_vector)))
		if ecc_vector[2]<0: omegaStar = 2*np.pi - omegaStar
	return m_per, aStar|units.m, eStar, iStar, OmegaStar, omegaStar

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

# M_max = max(m_bin_init+m_per_max, 3*m_per_max) 

# Outer binary parameters
ecc_out = 0.0         # Outer orbit eccentricity
inc_out = np.pi/6             # Outer orbit inclination
a_out = 0.5        # Outer semi-major axis in pc

# Start at pericentre
r = a_out * (1 - ecc_out)   # Spherical radius
R = r * np.cos(inc_out)     # Cylindrical radius
z = r * np.sin(inc_out)     # Cylindrical height

# output_file_name = os.path.dirname(os.path.abspath(__file__))+'/history-rtol7.txt'

# is the binary hard?
# a_h = G*(m1*m2/(m1+m2)|units.MSun)/4/sigma_rel(r|units.pc)**2
# # print(a_in/a_h.value_in(units.AU))
# Q=0.25
# print(((a_in|units.AU)/(64/5 * Q * G**3 * (k.m()|units.MSun)**3 / c**5 / (a_in|units.AU)**3)).value_in(units.yr))
# print(tau_0 (a_in|units.AU, m_bin|units.MSun, r|units.pc).value_in(units.yr))

# k = KeplerRing(ecc, inc, long_asc, arg_peri, [R1, z1, phi1], v1, a=a_in, m=m_bin, q=1)
# ts = np.linspace(0, 2*(t2-t1), 1000)
# k.integrate(ts, pot=pot, relativity=True, gw=True, tau_0=lambda *args: tau_0(args[0]|units.pc, m_bin|units.MSun, args[1]|units.pc).value_in(units.yr), random_number=1.4, rtol=1e-7, atol=1e-10)
# R2real, z2real, phi2real = k.r(t2-t1)
# v2real = k.v(t2-t1)

# E1 = (np.linalg.norm(v1)/_kms)**2/2 + evaluatePotentials(pot, R1/_pc, z1/_pc, phi=phi1, use_physical=False) 
# E2 = (np.linalg.norm(v2)/_kms)**2/2 + evaluatePotentials(pot, R2/_pc, z2/_pc, phi=phi2, use_physical=False) 
# E2real = (np.linalg.norm(v2real)/_kms)**2/2 + evaluatePotentials(pot, R2real/_pc, z2real/_pc, phi=phi2real, use_physical=False) 

def evolve_binary_noenc (input):

	t = input.t

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

	# Potential
	b = input.b
	m_total = input.m_total
	type = input.potential
	if type=="Plummer": pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 
	elif type=="Hernquist": pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc) 
	elif type=='Kepler': pot = KeplerPotential(amp=m_total*u.solMass)

	# Compute the correct v_phi at pericentre for the selected eccentricity and potential
	v_phi = ecc_to_vel(pot, ecc_out, [R, z, 0])

	# Define the KeplerRing
	k = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=m2/m1)
	k1 = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=m2/m1)

	# output_file = open(input.output_file, 'w+')
	# print('t[yr] R[pc] z phi v_R[km/s] v_z v_phi a[AU] m[MSun] q ecc inc long_asc arg_peri, outer_integration_time, tidal_time, inner_integration_time', file=output_file)
	# print(0, R, z, 0, 0, 0, v_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), file=output_file, flush=True)

	if type=='Kepler':
		T = 2*np.pi*np.sqrt((r|units.pc)**3/G/(m_total|units.MSun))
	else:
		T = 2*np.pi*(r|units.pc)/sigma(r|units.pc, type, m_total, b)	# approximate outer period
	n = max(int(input.n*t/(T.value_in(units.yr))), 100)	#number of points used to approximate the outer orbit
	# n=10
	ts = np.linspace(0, t, n)+input.t0
	rtol=input.rtol #1e-11
	atol= rtol*1e-3 #1e-14
	
	k.integrate(ts, pot=pot, relativity=input.relativity, gw=input.gw, tau_0=lambda *args: tau_0(args[0]|units.pc, k.m()|units.MSun, args[1]|units.pc, 50, type, m_total, b, m_per=m_per).value_in(units.yr), random_number=1e10, rtol=rtol, atol=atol, approximation=input.approximation, debug_file=input.output_file_2, points_per_period=input.n, ej_method='LSODA')
	if k.switch_to_gr:
		ts += k.t_fin
		k = KeplerRing(k.ecc_fin, k.inc_fin, k.long_asc_fin, k.arg_peri_fin, k.r(k.t_fin), k.v(k.t_fin), a=k.a_fin, m=k._m, q=k._q)
		k.integrate(ts, pot=pot, relativity=input.relativity, gw=input.gw, tau_0=lambda *args: tau_0(args[0]|units.pc, k.m()|units.MSun, args[1]|units.pc, 50, type, m_total, b, m_per=m_per).value_in(units.yr), random_number=1e10, rtol=rtol, atol=atol, approximation=2, debug_file=input.output_file_2, points_per_period=input.n)

	# print('da de di dOmega domega', file=output_file)
	# if k.merger: print('merger at', k.t_fin, file=output_file, flush=True)
	# else: print(k.t_fin, k.a_fin-a_in, k.ecc_fin-ecc, k.inc_fin-inc, k.long_asc_fin-long_asc, k.arg_peri_fin-arg_peri, k.outer_integration_time, k.tidal_time, k.inner_integration_time, file=output_file, flush=True)
	# print('epsilon =', k.epsilon_gr, file=output_file, flush=True)

def evolve_binary (input):
	# 0 - binary has survived until t_final
	# 1 - binary has merged
	# 2 - binary has been destroyed
	# 3 - maximum calculation time exceeded
	# 4 - maximum semimajor axis exceeded
	# 5 - binary has been ejected from the cluster

	t_final = input.t|units.yr
	m_per = input.m_per|units.MSun

	if input.includeWeakEncounters: Q_max_a = input.Q_max_a
	else: Q_max_a = Q_hybrid_a
	Q_min_a = input.Q_min_a
	
	# Potential
	b = input.b
	m_total = input.m_total
	type = input.potential
	if type=="Plummer": pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 
	elif type=="Hernquist": pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc) 

	if os.path.isfile(input.sameParameters):
		with open(input.sameParameters) as f:
			for line in f:
				data = line.split()
				if isfloat(data[0]): break
			t = float(data[0])|units.yr
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
		output_file = open(input.output_file, 'w+')
		print('t[yr] R[pc] z phi v_R[km/s] v_z v_phi a[AU] m[MSun] q ecc inc long_asc arg_peri random_number_0 dt[yr] n epsilon_gr t_orbit[s] |de|', file=output_file)	# |de| is the integral of |de/dt|_tidal
		print('perturber: m_per[MSun] Q[AU] eStar iStar OmegaStar omegaStar t_3body[s]', file=output_file)
		print(0, R, z, 0, 0, 0, v_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), file=output_file)
		output_file.flush()
	elif input.resume and os.path.isfile(input.output_file):
		with open(input.output_file) as f:
			for line in f: pass
			data = line.split()
			if data[1] == 'merger': return 1
			if data[1] == 'destroyed': return 2
			if data[1] == 'maximum' and data[2] == 'calculation': return 3
			if data[1] == 'maximum' and data[2] == 'semimajor': return 4
			if data[1] == 'ejected': return 5
			t = float(data[0])|units.yr
			if t>t_final: return 0
			for index in range(14):
				if not isfloat(data[index]):
					print("bad file:", input.output_file)
					print("bad line:", line)
					print("bad segment:", data[index])
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
		print('t[yr] R[pc] z phi v_R[km/s] v_z v_phi a[AU] m[MSun] q ecc inc long_asc arg_peri random_number_0 dt[yr] n epsilon_gr t_orbit[s] |de| |di|', file=output_file)	# |de| is the integral of |de/dt|_tidal
		print('perturber: m_per[MSun] Q[AU] eStar iStar OmegaStar omegaStar t_3body[s]', file=output_file)
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
		# integrate the orbit until the next flyby
		time0 = time.time()
		rng = default_rng()
		random_number = rng.exponential()
		random_number_0 = random_number
		r = np.sqrt(k.r()[0]**2+k.r()[1]**2)
		tau_0_value = tau_0 (k.a()|units.AU, k.m()|units.MSun, r|units.pc, Q_max_a, type, m_total, b, m_per=m_per) # doesn't really matter that V=0 here, it's only an estimate for the encounter time
		T = 2*np.pi*(r|units.pc)/sigma(r|units.pc, type, m_total, b)	# approximate outer period
		timeOrbit1 = time.time()
		Q = k._q / (1+k._q)**2
		t_gw = (k.a()|units.AU)/(64/5 * Q * G**3 * (k.m()|units.MSun)**3 / c**5 / (k.a()|units.AU)**3)
		dt = 1.1*min(tau_0_value*random_number, t_gw, (t_final-t))
		n = max(int(dt*input.n/T), 100)
		switch_to_gr = False
		approximation = input.approximation
		de_abs = 0
		di_abs = 0
		while (random_number>0):
			if switch_to_gr: approximation = 2
			ts = np.linspace(0, dt.value_in(units.yr), n+1)#100*n+1) #n is the number of time intervals
			k.integrate(ts, pot=pot, relativity=input.relativity, tidal_effects=input.tidal_effects, gw=input.gw, tau_0=lambda *args: tau_0(args[0]|units.pc, k.m()|units.MSun, args[1]|units.pc, Q_max_a, type, m_total, b, args[2]|units.kms, m_per=m_per).value_in(units.yr), random_number=random_number, rtol=rtol, atol=atol, approximation=approximation, debug_file=input.output_file_2, points_per_period=input.n) #, rtol=1e-3, atol=1e-6)
			t += k.t_fin|units.yr
			if k.merger: break
			random_number = k.probability
			outer_integration_time = k.outer_integration_time
			tidal_time = k.tidal_time
			inner_integration_time = k.inner_integration_time
			epsilon_gr = k.epsilon_gr
			de_abs += k.de_abs
			di_abs += k.di_abs
			if not switch_to_gr: switch_to_gr = k.switch_to_gr
			k = KeplerRing(k.ecc_fin, k.inc_fin, k.long_asc_fin, k.arg_peri_fin, k.r(k.t_fin), k.v(k.t_fin), a=k.a_fin, m=k._m, q=k._q)
			if t>=t_final: break
			R, z, phi = k.r()
			if np.sqrt(R**2+z**2) > 100: break
		timeOrbit2 = time.time()
		timeOrbit += timeOrbit2 - timeOrbit1
		R, z, phi = k.r()
		v_R, v_z, v_phi = k.v()
		if k.merger:
			print(t.value_in(units.yr), "merger", file=output_file)
			return 1
		print(t.value_in(units.yr), R, z, phi, v_R, v_z, v_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), random_number_0, dt.value_in(units.yr), n, epsilon_gr, time.time()-time0, de_abs, di_abs, file=output_file)
		output_file.flush()
		if t>=t_final: return 0
		if np.sqrt(R**2+z**2) > 100:
			print(t.value_in(units.yr), "ejected", file=output_file)
			return 5
		if input.includeEncounters:
			time3body = time.time()
			# sample the perturber parameters
			Q = 0|units.m
			while Q<=Q_min_a*(k.a()|units.AU):
				m_per, aStar, eStar, iStar, OmegaStar, omegaStar = sample_encounter_parameters (k.a()|units.AU, k.m()|units.MSun, np.sqrt(R**2+z**2)|units.pc, phi, Q_max_a, type, m_total, b, [v_R, v_phi, v_z], m_per=m_per)
				Q = aStar*(1-eStar)

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
				result = 0
				third_body_final = 2
				timeDistant1 = time.time()
				dv_binary, a_fin, e_fin, i_fin, Omega_fin, omega_fin = scattering_SA (m1|units.MSun, m2|units.MSun, k.a()|units.AU, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), m_per, aStar, eStar, iStar, OmegaStar, omegaStar)
				timeDistant2 = time.time()
				timeDistant += timeDistant2 - timeDistant1
				 
			print('perturber: ', m_per.value_in(units.MSun), Q.value_in(units.AU), eStar, iStar, OmegaStar, omegaStar, file=output_file)
			output_file.flush()
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
				if input.disableKicks:
					dv_x = 0
					dv_y = 0
					dv_z = 0
				else:
					dv_x = dv_binary[0].value_in(units.kms)	#velocity change during the scattering
					dv_y = dv_binary[1].value_in(units.kms)
					dv_z = dv_binary[2].value_in(units.kms)
				dv_R = (x*dv_x+y*dv_y)/R
				phi_unit_vector = np.cross([0, 0, 1], [x/R, y/R, 0])
				dv_phi = np.dot(phi_unit_vector, [dv_x, dv_y, dv_z])
				k = KeplerRing(e_fin, i_fin.value_in(units.rad), Omega_fin.value_in(units.rad), omega_fin.value_in(units.rad), [R, z, phi], [v_R+dv_R, v_z+dv_z, v_phi+dv_phi], a=a_fin.value_in(units.AU), m=m1+m2, q=min(m1/m2, m2/m1))
				m_bin = m1+m2
				print(t.value_in(units.yr), R, z, phi, v_R+dv_R, v_z+dv_z, v_phi+dv_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), time.time()-time3body, file=output_file)
				output_file.flush()
			else:
				print(t.value_in(units.yr), "something really weird happened: result =", result, file=output_file)
		
		timeTotal = time.time()-timeTotal1
		if timeTotal>input.tmax: 
			print(t.value_in(units.yr), 'maximum calculation time (', str(input.tmax), ' s) exceeded', file=output_file)
			output_file.flush()
			return 3
		if k.a() > input.a_max:
			print(t.value_in(units.yr), 'maximum semimajor axis (', str(input.a_max), ' AU) exceeded', file=output_file)
			output_file.flush()
			return 4

	# timeTotal2 = time.time()
	# print("total time", timeTotal2-timeTotal1, "s", file=output_file)
	# print("close interaction time", timeClose, "s", file=output_file)
	# print("distant interaction time", timeDistant, "s", file=output_file)
	# print("outer orbit integration time", timeOrbit, "s", file=output_file)
	output_file.close()

	return 0

def evolve_binary_encounters_only (input, n_enc, randomize):

	if input.includeWeakEncounters: Q_max_a = input.Q_max_a
	else: Q_max_a = Q_hybrid_a
	
	# Potential
	b = input.b
	m_total = input.m_total
	type = input.potential
	if type=="Plummer": pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 
	elif type=="Hernquist": pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc) 

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
	phi = 0

	# Compute the correct v_phi at pericentre for the selected eccentricity and potential
	v_phi = ecc_to_vel(pot, ecc_out, [R, z, 0])

	# Define the KeplerRing
	k = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=m2/m1)

	output_file = open(input.output_file, 'a')
	print('R[pc] z phi v_R[km/s] v_z v_phi a[AU] m[MSun] q ecc inc long_asc arg_peri t_3body[s]', file=output_file)
	print('perturber: m_per[MSun] Q[AU] eStar iStar OmegaStar omegaStar', file=output_file)
	print(R, z, phi, 0, 0, v_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), file=output_file)
	output_file.flush()

	rtol=input.rtol #1e-11
	atol= rtol*1e-3 #1e-14

	for i in range(n_enc):
		if randomize:
			arg_peri = 2*np.pi*np.random.random_sample()
			long_asc = 2*np.pi*np.random.random_sample()
			inc = np.arccos(2*np.random.random_sample()-1)
			k = KeplerRing(ecc, inc, long_asc, arg_peri, [R, z, 0], [0, 0, v_phi], a=a_in, m=m_bin, q=m2/m1)
			print(R, z, phi, 0, 0, v_phi, k.a(), k.m(), k._q, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), file=output_file)
			output_file.flush()
		time3body = time.time()
		# sample the perturber parameters
		m_per, aStar, eStar, iStar, OmegaStar, omegaStar = sample_encounter_parameters (k.a()|units.AU, k.m()|units.MSun, np.sqrt(R**2+z**2)|units.pc, phi, Q_max_a, type, m_total, b, [0, v_phi, 0], m_per=m_per)
		Q = aStar*(1-eStar)
		print('perturber: ', m_per.value_in(units.MSun), Q.value_in(units.AU), eStar, iStar, OmegaStar, omegaStar, file=output_file)
		output_file.flush()

		# perform the scattering
		q = k._q
		m1 = k.m()/(1+q)
		m2 = k.m()*q/(1+q)
		if Q<=Q_hybrid_a*(k.a()|units.AU):
			result, third_body_final, dv_binary, a_fin, e_fin, i_fin, Omega_fin, omega_fin, n_3body = scattering_hybrid (m1|units.MSun, m2|units.MSun, k.a()|units.AU, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), m_per, aStar, eStar, iStar, OmegaStar, omegaStar)
		else:
			result = 0
			third_body_final = 2
			dv_binary, a_fin, e_fin, i_fin, Omega_fin, omega_fin = scattering_SA (m1|units.MSun, m2|units.MSun, k.a()|units.AU, k.ecc(), k.inc(), k.long_asc(), k.arg_peri(), m_per, aStar, eStar, iStar, OmegaStar, omegaStar)

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
			k1 = KeplerRing(e_fin, i_fin.value_in(units.rad), Omega_fin.value_in(units.rad), omega_fin.value_in(units.rad), [R, z, phi], [v_R+dv_R, v_z+dv_z, v_phi+dv_phi], a=a_fin.value_in(units.AU), m=m1+m2, q=min(m1/m2, m2/m1))
			m_bin = m1+m2
			print(R, z, phi, v_R+dv_R, v_z+dv_z, v_phi+dv_phi, k1.a(), k1.m(), k1._q, k1.ecc(), k1.inc(), k1.long_asc(), k1.arg_peri(), time.time()-time3body, file=output_file)
			output_file.flush()

	output_file.close()

	return 0