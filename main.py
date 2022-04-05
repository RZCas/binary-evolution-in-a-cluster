from binary_evolution_with_flybys import inputParameters, evolve_binary, evolve_binary_noenc, a_h, sigma_rel, G, c
import numpy as np
import astropy.units as u
from astropy import constants
from amuse.lab import units
_G = constants.G.to(u.pc**3/u.solMass/u.yr**2).value
_c = constants.c.to(u.pc/u.yr).value

t = 1e6

# Inner binary parameters
a_in = 96.165752441              # Semi-major axis in AU
ecc = 0.5            	# Eccentricity
inc = 1.56905099754 #89.9 * np.pi/180 #Inclination with respect to the z-axis
long_asc = 0            # Longitude of the ascending node
arg_peri = -2.30340228521 #91.0 * np.pi/180# #    # Arugment of pericentre
m_tot = 7.08859466824
q = 0.304297545463
m1 = m_tot / (1+q)
m2 = m_tot * q / (1+q)

# print((a_in*(1-ecc)*u.au).to(u.pc), flush=True)

# Outer binary parameters
ecc_out = 0.2/3.2         # Outer orbit eccentricity
inc_out = 0             # Outer orbit inclination
a_out = 1.6        # Outer semi-major axis in pc

folder = 'output/perpendicular-test/'

output_file = folder + '1-noenc.txt'#folder + '0nogw.txt'#folder+'a_in='+str(a_in)+'_e_in='+str(ecc)+'_norelativity.txt'####
output_file_2 = folder + '1-noenc-evolution.txt'#folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_norelativity_evolution.txt'

# forcePrecise = False
# forceApproximate = not forcePrecise
# if forcePrecise:
# 	output_file = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_precise.txt'
# 	output_file_2 = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_precise_evolution.txt'
# if forceApproximate:
# 	output_file = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_approximate.txt'
# 	output_file_2 = folder + 'a_in='+str(a_in)+'_e_in='+str(ecc)+'_approximate_evolution.txt'

rtol=1e-11
potential = "Hernquist"
b = 1
m_total = 1e6
tmax = 5e20

# Q=0.25
# print("t_gw = %.2e" % (((a_in|units.AU)/(64/5 * Q * G**3 * ((m1+m2)|units.MSun)**3 / c**5 / (a_in|units.AU)**3)).value_in(units.yr)/(1+73/24*ecc**2+37/96*ecc**4)*(1-ecc**2)**3.5))
# print(a_h(m1,m2,a_out))
# print(sigma_rel(a_out|units.pc).value_in(units.kms))
# print("t_outer = %.2e" % (np.sqrt(G*((m1+m2)|units.MSun) * (1e6|units.MSun) / (a_out|units.pc)**3).value_in(units.yr)))

input = inputParameters(t=t, a_out=a_out, e_out=ecc_out, inc_out=inc_out, m1=m1, m2=m2, a=a_in, e=ecc, i=inc, Omega=long_asc, omega=arg_peri, output_file=output_file, output_file_2=output_file_2, potential=potential, b=b, m_total = m_total, rtol=rtol, tmax=tmax, 
	approximation=1,
	resume=False, 
	includeEncounters=True, 
	includeWeakEncounters=True,
	relativity=True,
	gw=True, 
	n=300)
# evolve_binary(input)
evolve_binary_noenc(input)