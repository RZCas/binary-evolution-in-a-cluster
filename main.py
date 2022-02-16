from binary_evolution_with_flybys import inputParameters, evolve_binary, evolve_binary_noenc, approximation_test, a_h, sigma_rel, G, c
import numpy as np
import astropy.units as u
from astropy import constants
from amuse.lab import units
_G = constants.G.to(u.pc**3/u.solMass/u.yr**2).value
_c = constants.c.to(u.pc/u.yr).value

t = 1e8

# Inner binary parameters
a_in = 30              # Semi-major axis in AU
ecc = 0.9999            	# Eccentricity
inc = 1           # Inclination with respect to the z-axis
long_asc = 0            # Longitude of the ascending node
arg_peri = 1.5    # Arugment of pericentre
m1 = 5
m2 = 5

# print((a_in*(1-ecc)*u.au).to(u.pc), flush=True)

# Outer binary parameters
ecc_out = 0.1         # Outer orbit eccentricity
inc_out = 0.5             # Outer orbit inclination
a_out = 0.5        # Outer semi-major axis in pc

forcePrecise = False
forceApproximate = not forcePrecise
if forcePrecise:
	output_file = 'output/epsilon_gr_test_2/a_in='+str(a_in)+'_e_in='+str(ecc)+'_precise.txt'
	output_file_2 = 'output/epsilon_gr_test_2/a_in='+str(a_in)+'_e_in='+str(ecc)+'_precise_evolution.txt'
if forceApproximate:
	output_file = 'output/epsilon_gr_test_2/a_in='+str(a_in)+'_e_in='+str(ecc)+'_approximate.txt'
	output_file_2 = 'output/epsilon_gr_test_2/a_in='+str(a_in)+'_e_in='+str(ecc)+'_approximate_evolution.txt'

rtol=1e-11
potential = "Plummer"
tmax = 5e20

Q=0.25
print("t_gw = %.2e" % (((a_in|units.AU)/(64/5 * Q * G**3 * ((m1+m2)|units.MSun)**3 / c**5 / (a_in|units.AU)**3)).value_in(units.yr)/(1+73/24*ecc**2+37/96*ecc**4)*(1-ecc**2)**3.5))
# print(a_h(m1,m2,a_out))
# print(sigma_rel(a_out|units.pc).value_in(units.kms))
# print("t_outer = %.2e" % (np.sqrt(G*((m1+m2)|units.MSun) * (1e6|units.MSun) / (a_out|units.pc)**3).value_in(units.yr)))

input = inputParameters(t=t, a_out=a_out, e_out=ecc_out, inc_out=inc_out, m1=m1, m2=m2, a=a_in, e=ecc, i=inc, Omega=long_asc, omega=arg_peri, output_file=output_file, output_file_2=output_file_2, potential=potential, rtol=rtol, tmax=tmax, 
	forcePrecise=forcePrecise,
	forceApproximate=forceApproximate,
	resume=False, 
	includeEncounters=True, 
	includeWeakEncounters=True, 
	n=10)
# evolve_binary(input)
evolve_binary_noenc(input)
# approximation_test(input)
# print(x)

# Q=0.25
# de = -304/15 * Q * _G**3 * (m1+m2)**3 / _c**5 / (a_in*u.au).to(u.pc).value**4 * (1+121/304*ecc**2)/(1-ecc**2)**2.5 * ecc * t
# da = (64/5 * Q * _G**3 * (m1+m2)**3 / _c**5 / (a_in*u.au).to(u.pc).value**3 * (1+73/24*ecc**2+37/96*ecc**4)/(1-ecc**2)**3.5 * t * u.pc).to(u.au).value
# print(de)