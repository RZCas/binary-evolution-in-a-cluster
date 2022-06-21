from binary_evolution_with_flybys import inputParameters, evolve_binary, evolve_binary_noenc, evolve_binary_encounters_only, sample_encounter_parameters, a_h, sigma, G, c, tau_0
import numpy as np
import astropy.units as u
import statistics
import matplotlib.pyplot as plt
from astropy import constants
from amuse.lab import units

_G = constants.G.to(u.pc**3/u.solMass/u.yr**2).value
_c = constants.c.to(u.pc/u.yr).value

t = 3.6e9

# Inner binary parameters
a_in = 10              # Semi-major axis in AU
ecc = 0.01            	# Eccentricity
inc = 0 #89.9 * np.pi/180 #Inclination with respect to the z-axis
long_asc = 1            # Longitude of the ascending node
arg_peri = -0.7 #91.0 * np.pi/180# #    # Arugment of pericentre
m_tot = 30
q = 1
m1 = m_tot / (1+q)
m2 = m_tot * q / (1+q)

# Outer binary parameters
ecc_out = 0.2/3.2         # Outer orbit eccentricity
inc_out = 0             # Outer orbit inclination
a_out = 1.6        # Outer semi-major axis in pc

folder = 'output/perpendicular-hard-plummer/'

# output_file = folder + 'ejected_0.txt'#folder + '0nogw.txt'#folder+'a_in='+str(a_in)+'_e_in='+str(ecc)+'_norelativity.txt'####
output_file = folder + '6.txt'
# output_file_2 = folder + 'ejected_0-evolution.txt'#folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_norelativity_evolution.txt'
output_file_2 = folder + '6-evolution.txt'#folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_norelativity_evolution.txt'

# forcePrecise = False
# forceApproximate = not forcePrecise
# if forcePrecise:
# 	output_file = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_precise.txt'
# 	output_file_2 = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_precise_evolution.txt'
# if forceApproximate:
# 	output_file = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_approximate.txt'
# 	output_file_2 = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_approximate_evolution.txt'

rtol=1e-11
potential = "Plummer"
b = 1
m_total = 1e6
tmax = 5e20

a = 84.0480389205|units.AU
m_bin = 8.31464608615|units.MSun
r = 1.5|units.pc
# print(tau_0 (a, m_bin, r, Q_max_a=50, type=potential, m_total=1e6, b=1).value_in(units.yr))

# Q=0.25
# print("t_gw = %.2e" % (((a_in|units.AU)/(64/5 * Q * G**3 * ((m1+m2)|units.MSun)**3 / c**5 / (a_in|units.AU)**3)).value_in(units.yr)/(1+73/24*ecc**2+37/96*ecc**4)*(1-ecc**2)**3.5))
# print(a_h(m1,m2,a_out))
# print(sigma(a_out|units.pc).value_in(units.kms))
# print("t_outer = %.2e" % (np.sqrt(G*((m1+m2)|units.MSun) * (1e6|units.MSun) / (a_out|units.pc)**3).value_in(units.yr)))

input = inputParameters(t=t, a_out=a_out, e_out=ecc_out, inc_out=inc_out, m1=m1, m2=m2, a=a_in, e=ecc, i=inc, Omega=long_asc, omega=arg_peri, output_file=output_file, output_file_2=output_file_2, potential=potential, b=b, m_total = m_total, rtol=rtol, tmax=tmax, 
	approximation=0,
	resume=True, 
	includeEncounters=True, 
	includeWeakEncounters=True,
	relativity=True,
	gw=True, 
	n=300,
	sameParameters='',
	Q_max_a=50)
n_enc=1
# evolve_binary_encounters_only(input, n_enc, randomize=False)
evolve_binary(input)

# result = []
# for n in range (10000):
# 	m_per, aStar, eStar, iStar, OmegaStar, omegaStar = sample_encounter_parameters (a, m_bin, r, phi, Q_max_a=50, type="Plummer", m_total=4e6, b=1, v_bin=[1000,0,0])
# 	result.append(OmegaStar)
# plt.hist(result)
# plt.show()
# print(statistics.mean(result))