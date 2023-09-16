import numpy as np
from scipy.integrate import quad
# Calculate A_*

def Phi(r, b, type):
	if type=="Hernquist":
		return -1/(r+b)
	else:
		return -1/np.sqrt(r**2+b**2)

def A(r, b, type):
	if type=="Hernquist":
		return 0.5*b**3*(3*b+r)/r/(r+b)**3
	else:
		return 0.5*b**3*(4*b**2+r**2)/(r**2+b**2)**2.5

def integrand_1 (r, b, type, E, L):
	return A(r, b, type)/np.sqrt(2*(E-Phi(r, b, type))-(L/r)**2)

def integrand_2 (r, b, type, E, L):
	return 1/np.sqrt(2*(E-Phi(r, b, type))-(L/r)**2)

def A_averaged(r_p, r_a, potential_type, b):
	if r_a/r_p-1<0.01:
		return A((r_p+r_a)/2, b, potential_type)
	else:
		E = (r_a**2*Phi(r_a, b, potential_type) - r_p**2*Phi(r_p, b, potential_type)) / (r_a**2 - r_p**2)
		L = np.sqrt(2*(Phi(r_a, b, potential_type)-Phi(r_p, b, potential_type))/(1/r_p**2-1/r_a**2))
		return quad(integrand_1, r_p, r_a, args=(b, potential_type, E, L))[0] / quad(integrand_2, r_p, r_a, args=(b, potential_type, E, L))[0]

print(A(r=4,b=1,type='Hernquist'))