import os
import glob
import numpy as np
import astropy.units as u
from galpy.potential import evaluaterforces, evaluatePotentials, PotentialError, KeplerPotential, TwoPowerTriaxialPotential, PlummerPotential, HernquistPotential
# from binary_evolution_with_flybys import a_h, sigma
from amuse.lab import units, constants 
from binary_evolution.tools import rarp
from binary_evolution import KeplerRing
_pc = 8000
_kms = 220
G = constants.G
c = constants.c
t_H = 1.4e10|units.yr
H = 15
potential = "Hernquist"

m_per = 1|units.MSun
# def tau_0_factor (a, m_bin, r, Q_max_a=50, type="Plummer", m_total=4e6, b=1, V=0|units.kms):
# 	Q_max = Q_max_a * a
# 	v0 = np.sqrt(G*(m_bin+m_per)/Q_max)
# 	sigma_rel = np.sqrt(sigma(r, type, m_total, b)**2 + V**2)
# 	return (Q_max**2*(1+(v0/sigma_rel)**2))**-1

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
# pyplot.yscale('log')

t_max=1e80
a_out = 2
m_total = 1e6
b = 2
A_ast = 0.3 #A_* for Hernquist potential
pot = HernquistPotential(amp=2*m_total*u.solMass, a=b*u.pc)
# pot = PlummerPotential(amp=m_total*u.solMass, b=b*u.pc) 

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

# old:
# perpendicular-soft-hernquist
# 6 - merged
# 16 - abandoned
# 15 - destroyed
# 11 - merged after an exchange

# new:
# m1=m2=10, mtotal=1e6, hernquist
# 0 - abandoned
# 2 - merged after an exchange
# 7 - merged
# 19 - destroyed
# nokicks, 0 - no kicks
# evolution-hernquist,m_total=1e5,b=1,a_out=4,i=89.9,nokicks,a_in=300-5 - the case where tidal effects dominate
# uniform_mtotal=1e5_hernquist-9 - ejected after an exchange
# ejected/wide_range_mtotal=1e5_plummer_b=1-145 - ejected

root_dir = ["output/m1=m2=10/mtotal=1e6/",
			"output/m1=m2=10/mtotal=1e6/",
			"output/m1=m2=10/mtotal=1e6/",
			"output/m1=m2=10/mtotal=1e6/",
			"output/m1=m2=10/mtotal=1e6_nokicks/",
			'output/hernquist,m_total=1e5,b=1,a_out=4,i=89.9,nokicks,a_in=300/',
			'output/uniform_mtotal=1e5_hernquist/', 
			'output/ejected/wide_range_mtotal=1e5_plummer_b=1-']
pot = [HernquistPotential(amp=2*1e6*u.solMass, a=2*u.pc),
		HernquistPotential(amp=2*1e6*u.solMass, a=2*u.pc),
		HernquistPotential(amp=2*1e6*u.solMass, a=2*u.pc),
		HernquistPotential(amp=2*1e6*u.solMass, a=2*u.pc),
		HernquistPotential(amp=2*1e6*u.solMass, a=2*u.pc),
		HernquistPotential(amp=2*1e5*u.solMass, a=1*u.pc),
		HernquistPotential(amp=2*1e5*u.solMass, a=2*u.pc),
		PlummerPotential(amp=1e5*u.solMass, b=1*u.pc)]
nokicks = [False,False,False,False,True,True,False,False]
indices = [0,2,7,19,0,5,9,145]
fileNames = ['abandoned','exchange','merged','destroyed', 'nokicks', 'tidalDominated', 'ejectedExchange', 'ejected']
# for i in [-2]:#range(len(indices)):
# 	filepath_epsilon = "output/for the paper/"+fileNames[i]+"-epsilon.txt"
# 	with open(filepath_epsilon) as f:
# 		for line in f:
# 			data = line.split()
# 			if data[1]=='nan':
# 				R = float(data[1])
# 				z = float(data[2])
# 				phi = float(data[3])
# 				v_R = float(data[4])
# 				v_z = float(data[5])
# 				v_phi = float(data[6])
# 				a_0 = float(data[7])
# 				m = float(data[8])
# 				q = float(data[9])
# 				e_0 = float(data[10])
# 				i_0 = float(data[11])
# 				Omega_0 = float(data[12])
# 				omega_0 = float(data[13])
# 				k = KeplerRing(e_0, i_0, Omega_0, omega_0, [R, z, phi], [v_R, v_z, v_phi], a=a_0, m=m, q=q)
# 				epsilon_real = np.log10(k.epsilon_gr_real(pot[i], num_periods=200))
# 				print(float(data[0]), epsilon_real)

# t, R, z, phi, v_R, v_z, v_phi, a_0, m, q, e_0, i_0, Omega_0, omega_0 = 998870250.135, 0.0127927546603, -0.032844278183, 0.35972977679, -7.47337299039, -0.43285167261, -0.302352377428, 0.141112877355, 11.3922686477, 0.0962253800297, 0.494480953171, 2.45574246368, 2.49851303152, -1.80625215855
t, R, z, phi, v_R, v_z, v_phi, a_0, m, q, e_0, i_0, Omega_0, omega_0 = 934967023.623, 0.162330925889, -0.0396360845983, 0.656397464457, -6.85600442125, -1.86451112498, 0.1206669983, 0.141122747975, 11.3922686477, 0.0962253800297, 0.496548393201, 2.45885203961, 2.50176580453, 0.556411560659
num_periods = 100
k = KeplerRing(e_0, i_0, Omega_0, omega_0, [R, z, phi], [v_R, v_z, v_phi], a=a_0, m=m, q=q)
epsilon_real = np.log10(k.epsilon_gr_real(pot[-2], num_periods=num_periods))
print('t =', t, "num_periods =", num_periods, ", epsilon_real =", epsilon_real)