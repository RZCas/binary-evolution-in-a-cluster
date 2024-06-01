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

# Example 7 from our paper (Fig. 11)
t = 1.6e6
a_in = 300             			# Semi-major axis in AU
ecc = 0.5            			# Eccentricity
inc = 89.9 * np.pi/180 			# Inclination with respect to the z-axis
long_asc = 0*np.pi/2            # Longitude of the ascending node
arg_peri = 51.9 * np.pi/180		# Arugment of pericentre
m_tot = 20
q = 1
m1 = m_tot / (1+q)
m2 = m_tot * q / (1+q)
ecc_out = 0         # Outer orbit eccentricity
inc_out = 0             # Outer orbit inclination
a_out = 4        # Outer semi-major axis in pc
potential = "Hernquist"
b = 1
m_total = 1e5
output_file = 'example-output.txt'
output_file_2 = ''

rtol = 1e-11
tmax = 5e20

input = inputParameters(t=t, a_out=a_out, e_out=ecc_out, inc_out=inc_out, m1=m1, m2=m2, a=a_in, e=ecc, i=inc, Omega=long_asc, omega=arg_peri, output_file=output_file, output_file_2=output_file_2, potential=potential, b=b, m_total = m_total, rtol=rtol, tmax=tmax, t0=t0,
	approximation=0,
	resume=False, 
	includeEncounters=True, 
	includeWeakEncounters=True,
	relativity=True,
	gw=True, 
	# n=300, #n=30
	# sameParameters='output/wide_range_1/2.txt',
	# Q_max_a=50,
	disableKicks=False,
	m_per=1,
	tidal_effects=True)

# n_enc=1
# evolve_binary_encounters_only(input, n_enc, randomize=False)
evolve_binary(input)
# evolve_binary_noenc(input)