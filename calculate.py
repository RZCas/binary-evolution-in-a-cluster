import numpy as np

gamma = 0.42
e = 0.5
i = 89.9 * np.pi/180
omega = 91 * np.pi/180
eps_gr = 0.0691200754099

def e_max (e, i, omega, eps_gr):

	# returns log_10(1-e_max)

	theta = (1-e**2)*np.cos(i)**2
	D = e**2 * (1 + 10*gamma/(1-5*gamma) * np.sin(i)**2*np.sin(omega)**2) - eps_gr/(3*(1-5*gamma)*np.sqrt(1-e**2))
	sigma = (1+5*gamma)/2 + 5*gamma*theta + (5*gamma-1)/2*D
	jplus2 = (sigma + np.sqrt(sigma**2 - 10*gamma*theta*(1+5*gamma))) / (1+5*gamma)
	jminus2 = (sigma - np.sqrt(sigma**2 - 10*gamma*theta*(1+5*gamma))) / (1+5*gamma)
	eps_weak = 6*(1+5*gamma)*jplus2*np.sqrt(jminus2) #10**-1.6
	eps_strong = 3*(1+5*gamma) #10
	jmin = (eps_gr + np.sqrt(eps_gr**2+eps_weak**2)) / (2*jplus2*eps_strong)

	return np.log10 (1-np.sqrt(1-jmin**2))

print (e_max(e, i, omega, eps_gr))

# cluster4 - ok