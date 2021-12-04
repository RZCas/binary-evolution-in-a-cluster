from binary_evolution_with_flybys import inputParameters, evolve_binary, approximation_test, detailed_output 
import numpy as np

t = 1e4

# Inner binary parameters
a_in = 0.01              # Semi-major axis in AU
ecc = 0.9            	# Eccentricity
inc = 1           # Inclination with respect to the z-axis
long_asc = 0            # Longitude of the ascending node
arg_peri = 1.5    # Arugment of pericentre
m1 = 5
m2 = 5

# Outer binary parameters
ecc_out = 0.0         # Outer orbit eccentricity
inc_out = 0.5             # Outer orbit inclination
a_out = 0.5        # Outer semi-major axis in pc

output_file='output/detailed-output.pdf'

input = inputParameters(t=t, a_out=a_out, e_out=ecc_out, inc_out=inc_out, m1=m1, m2=m2, a=a_in, e=ecc, i=inc, Omega=long_asc, omega=arg_peri, output_file=output_file, forcePrecise=False)
# evolve_binary(input)
# approximation_test(input)
detailed_output(input)