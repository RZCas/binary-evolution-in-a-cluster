from binary_evolution_with_flybys import inputParameters, evolve_binary, evolve_binary_noenc, evolve_binary_encounters_only, sample_encounter_parameters, a_h, sigma, G, c, tau_0
import numpy as np
import astropy.units as u
import statistics
import matplotlib.pyplot as plt
from astropy import constants
from amuse.lab import units

_G = constants.G.to(u.pc**3/u.solMass/u.yr**2).value
_c = constants.c.to(u.pc/u.yr).value

t0 = 0

# example 1 - continuation
# t = 15e6
# a_in = 25.0363903724             			# Semi-major axis in AU
# ecc = 0.999997566852            			# Eccentricity
# inc = 1.49994506781			# Inclination with respect to the z-axis
# long_asc = -2.07924540485            	# Longitude of the ascending node
# arg_peri = 0.756529079173		# Arugment of pericentre
# m_tot = 20
# q = 1
# m1 = m_tot / (1+q)
# m2 = m_tot * q / (1+q)
# ecc_out = 0.1/0.3         # Outer orbit eccentricity
# inc_out = 0             # Outer orbit inclination
# a_out = 0.15        # Outer semi-major axis in pc
# potential = "Kepler"
# b = 1
# m_total = 4e6
# t0 = 615037.5659064837

# example 1
# t = 15e6
# a_in = 30             			# Semi-major axis in AU
# ecc = 0.2            			# Eccentricity
# inc = 89.75 * np.pi/180 			# Inclination with respect to the z-axis
# long_asc = 1*np.pi/2            	# Longitude of the ascending node
# arg_peri = 91.0 * np.pi/180		# Arugment of pericentre
# m_tot = 20
# q = 1
# m1 = m_tot / (1+q)
# m2 = m_tot * q / (1+q)
# ecc_out = 0.1/0.3         # Outer orbit eccentricity
# inc_out = 0             # Outer orbit inclination
# a_out = 0.15        # Outer semi-major axis in pc
# potential = "Kepler"
# b = 1
# m_total = 4e6

# example 2
# t = 8e9
# a_in = 250             			# Semi-major axis in AU
# ecc = 0.6            			# Eccentricity
# inc = 89.8 * np.pi/180 			# Inclination with respect to the z-axis
# long_asc = 0*np.pi/2            	# Longitude of the ascending node (Omega)
# arg_peri = 91.0 * np.pi/180		# Arugment of pericentre (omega)
# m_tot = 20
# q = 1
# m1 = m_tot / (1+q)
# m2 = m_tot * q / (1+q)
# ecc_out = 0.2/3.2         # Outer orbit eccentricity
# inc_out = 0             # Outer orbit inclination
# a_out = 1.6        # Outer semi-major axis in pc
# potential = "Hernquist"
# b = 1
# m_total = 1e6

# example 3
# t = 1.1e10
# a_in = 49             			# Semi-major axis in AU
# ecc = 0.5            			# Eccentricity
# inc = 89.9 * np.pi/180 			# Inclination with respect to the z-axis
# long_asc = 0*np.pi/2            	# Longitude of the ascending node
# arg_peri = 91.0 * np.pi/180		# Arugment of pericentre
# m_tot = 2.8
# q = 1
# m1 = m_tot / (1+q)
# m2 = m_tot * q / (1+q)
# ecc_out = 0.2/3.2         # Outer orbit eccentricity
# inc_out = 0             # Outer orbit inclination
# a_out = 1.6        # Outer semi-major axis in pc
# potential = "Hernquist"
# b = 1
# m_total = 1e6


# example 4
t = 9e8
a_in = 30             			# Semi-major axis in AU
ecc = 0.5            			# Eccentricity
inc = 89.3 * np.pi/180 			# Inclination with respect to the z-axis
long_asc = 0*np.pi/2            # Longitude of the ascending node
arg_peri = 45.5 * np.pi/180		# Arugment of pericentre
m_tot = 20
q = 1
m1 = m_tot / (1+q)
m2 = m_tot * q / (1+q)
ecc_out = 0.4/0.6         # Outer orbit eccentricity
inc_out = 0             # Outer orbit inclination
a_out = 0.3        # Outer semi-major axis in pc
potential = "Hernquist"
b = 1
m_total = 10e6

t = 1e1

folder = 'output/test/'

output_file = folder + 'test.txt'
output_file_2 = ''

# forcePrecise = False
# forceApproximate = not forcePrecise
# if forcePrecise:
# 	output_file = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_precise.txt'
# 	output_file_2 = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_precise_evolution.txt'
# if forceApproximate:
# 	output_file = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_approximate.txt'
# 	output_file_2 = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_approximate_evolution.txt'

rtol = 1e-11
tmax = 5e20

# a = 84.0480389205|units.AU
# m_bin = 8.31464608615|units.MSun
# r = 1.5|units.pc
# print(tau_0 (a, m_bin, r, Q_max_a=50, type=potential, m_total=1e6, b=1).value_in(units.yr))

# Q=0.25
# print("t_gw = %.2e" % (((a_in|units.AU)/(64/5 * Q * G**3 * ((m1+m2)|units.MSun)**3 / c**5 / (a_in|units.AU)**3)).value_in(units.yr)/(1+73/24*ecc**2+37/96*ecc**4)*(1-ecc**2)**3.5))
# print(a_h(m1,m2,a_out))
# print(sigma(a_out|units.pc).value_in(units.kms))
# print("t_outer = %.2e" % (np.sqrt(G*((m1+m2)|units.MSun) * (1e6|units.MSun) / (a_out|units.pc)**3).value_in(units.yr)))

input = inputParameters(t=t, a_out=a_out, e_out=ecc_out, inc_out=inc_out, m1=m1, m2=m2, a=a_in, e=ecc, i=inc, Omega=long_asc, omega=arg_peri, output_file=output_file, output_file_2=output_file_2, potential=potential, b=b, m_total = m_total, rtol=rtol, tmax=tmax, t0=t0,
	approximation=0,
	resume=False, 
	includeEncounters=True, 
	includeWeakEncounters=True,
	relativity=True,
	gw=True, 
	n=300, #n=30
	# sameParameters='output/wide_range_1/2.txt',
	# Q_max_a=50,
	disableKicks=False,
	m_per=0.1)
# n_enc=1
# evolve_binary_encounters_only(input, n_enc, randomize=False)
evolve_binary(input)

# result = []
# for n in range (10000):
# 	m_per, aStar, eStar, iStar, OmegaStar, omegaStar = sample_encounter_parameters (a, m_bin, r, phi, Q_max_a=50, type="Plummer", m_total=4e6, b=1, v_bin=[1000,0,0])
# 	result.append(OmegaStar)
# plt.hist(result)
# plt.show()
# print(statistics.mean(result))