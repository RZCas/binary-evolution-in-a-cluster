import numpy as np
from amuse.lab import constants, units
from binary_evolution_with_flybys import a_h
G = constants.G
c = constants.c

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

Q = 14.7374214658|units.AU
e = 2.49391316838
a = Q/(1-e)
m = 15.24176866269581|units.MSun
v = np.sqrt(-G*m/a)
# print(v.value_in(units.kms))

# print(sigma_rel (r=3|units.pc, type="Hernquist", m_total=1e6, b=1).value_in(units.kms))

# GW lifetime for e=0
m1 = 10|units.MSun
m2 = 10|units.MSun
T = 1.4e10|units.yr
a = (256 * G**3 * m1 * m2 * (m1+m2) * T / 5 / c**5)**(1/4)
# print(a.value_in(units.AU))
a = 0.0163828779997|units.AU
T = 5 * a**4 * c**5 / (256*G**3*m1*m2*(m1+m2))
# print(T.value_in(units.Gyr))

# 'destroyed' encounter parameters
a = 192.187989792|units.AU
Q = 2.79899533934|units.AU
eStar = 1.15129217311
m_per = 1|units.MSun
aStar = Q / (1 - eStar)
v = np.sqrt(-G*(m1+m2+m_per)/aStar)
v_b = np.sqrt(G*(m1+m2)/a)
# print(Q/a, v.value_in(units.kms), v_b.value_in(units.kms), v/v_b)

# 'exchange' encounter parameters
a = 8.99502651881|units.AU
Q = 0.788581418745|units.AU
eStar = 1.0009594839
m_per = 1|units.MSun
aStar = Q / (1 - eStar)
v = np.sqrt(-G*(m1+m2+m_per)/aStar)
v_b = np.sqrt(G*(m1+m2)/a)
# print(Q/a, v.value_in(units.kms), v_b.value_in(units.kms), v/v_b)

# paper II, eq. 34
A = 0.007
M = 1e6
b = 1
a = 300
m12 = 20
t1 = 1.7 * (A/0.5)**-1 * (M/1e5)**-1 * b**3 * (m12)*0.5 * (a/10)**-1.5
# print(t1) 

m1 = 10
m2 = 10
r = 2
ah5 = a_h(m1, m2, r, type="Hernquist", m_total=1e5, b=2)
ah6 = a_h(m1, m2, r, type="Hernquist", m_total=1e6, b=2)
print(ah5, ah6)