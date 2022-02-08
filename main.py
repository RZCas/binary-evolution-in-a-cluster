from binary_evolution_with_flybys import inputParameters, evolve_binary, evolve_binary_noenc, approximation_test, a_h, sigma_rel
import numpy as np
import astropy.units as u
from astropy import constants
_G = constants.G.to(u.pc**3/u.solMass/u.yr**2).value
_c = constants.c.to(u.pc/u.yr).value

t = 2e5

# Inner binary parameters
a_in = 1              # Semi-major axis in AU
ecc = 0.05            	# Eccentricity
inc = 1           # Inclination with respect to the z-axis
long_asc = 0            # Longitude of the ascending node
arg_peri = 1.5    # Arugment of pericentre
m1 = 5
m2 = 5

# print((a_in*(1-ecc)*u.au).to(u.pc), flush=True)

# Outer binary parameters
ecc_out = 0.1         # Outer orbit eccentricity
inc_out = 0.5             # Outer orbit inclination
a_out = 0.5 #0.5        # Outer semi-major axis in pc

# output_file = 'output/a_dependence_3/0_test.txt'
output_file = 'output/test-noenc.txt'
output_file_2 = 'output/test.pdf'
#'output/gr_ratio_dependence_hernquist_aout=01-2-rtol12.pdf'

rtol=1e-11
potential = "Plummer"
tmax = 5e20

print(a_h(m1,m2,a_out))
print(sigma_rel)
return

input = inputParameters(t=t, a_out=a_out, e_out=ecc_out, inc_out=inc_out, m1=m1, m2=m2, a=a_in, e=ecc, i=inc, Omega=long_asc, omega=arg_peri, output_file=output_file, output_file_2=output_file_2, potential=potential, rtol=rtol, tmax=tmax, 
	forcePrecise=True,
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