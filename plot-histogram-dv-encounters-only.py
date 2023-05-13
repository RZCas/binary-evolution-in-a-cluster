import os
import glob
import numpy as np
import statistics
from binary_evolution_with_flybys import a_h, normalize
from amuse.lab import units, constants 

G = constants.G
c = constants.c
t_H = 1.4e10|units.yr
H = 15
potential = "Hernquist"

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
a_out = 1.6
m_total = 1e6
b = 1
A_ast = 0.3 #A_* for Hernquist potential
# print(a_h(10, 10, a_out, type="Hernquist", m_total=1e6, b=1))

def a_tidal (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	return ((24*G**2*(m|units.MSun)**2/c**2/A)**0.25).value_in(units.AU)
def a_tsec01tH (m, m_cl, b):
	A = A_ast*G*(m_cl|units.MSun)/(b|units.pc)**3
	t1 = 0.1*t_H
	return (((8/(3*A*t1))**2*G*(m|units.MSun))**(1/3)).value_in(units.AU)

filepath = 'output/dynamical_friction_Qmax=200.txt'
color = 'k'
t = []
theta = []
e = []
cosi = []
a = []
de = []
da = []
dl = []
dv_parallel = []
dv_perp_in_orbit = []
dv_perp_to_orbit = []		
t_0 = 0
theta_previous = 0
lineNumber=0
dv_total = 0
with open(filepath) as f:
	for line in f:
		lineNumber+=1
		data = line.split()
		if len(data) > 1:
			if data[0] == 'perturber:' and lineNumber>2:
				Q = float(data[2])
			if lineNumber == 3:
				m = float(data[7])
				q = float(data[8])
				m1 = m/(1+q)
				m2 = m*q/(1+q)
				a_i = float(data[6])
				v_R = float(data[3])
				v_z = float(data[4])
				v_phi = float(data[5])
				v_0 = np.array([v_R, v_phi, v_z])
				R = float(data[0])
				z = float(data[1])
				phi = float(data[2])
				n = [R, 0, z]	#direction from the cluster centre
				l = np.cross(n, v_0)	#outer orbit angular momentum
				k = np.cross(v_0, l)	#vector, perp. to v and lying in the orbital plane (points away from the centre)
			if isfloat(data[0]) and isfloat(data[1]): #non-randomized
			# if lineNumber % 3 == 0 and lineNumber>3: #randomized angles
				v_R = float(data[3])
				v_z = float(data[4])
				v_phi = float(data[5])	
				# x, y, z axes are [R, phi, z]
				v = np.array([v_R, v_phi, v_z])
				if lineNumber > 3 and Q/a_i>80 and Q/a_i<160:
					dv = v - v_0
					dv_parallel.append(np.dot(dv, normalize(v_0)))
					dv_total += dv_parallel[-1] 
					dv_perp_to_orbit.append(np.dot(dv, normalize(l)))
					dv_perp_in_orbit.append(np.dot(dv, normalize(k)))

# figure = pyplot.figure(figsize=(4.5*1.5, 3.5*1.5))
# figure.suptitle(f'$m_1$ = {m1:.1f} $M_\odot$, $m_2$ = {m2:.1f} $M_\odot$, $a$ = {a_i:.1f} AU, $v$ = {np.linalg.norm(v_0):.1f} km/s')

# plot = figure.add_subplot(1,1,1)
# pyplot.xlabel(r'$\delta v$ [km/s]')
# pyplot.yscale('log', nonpositive='clip')
# bins=100
# range=(-2,2)
# plot.hist (dv_parallel, bins, range, histtype='step', color='k', label=r'parallel to $v$')
# plot.hist (dv_perp_to_orbit, bins, range, histtype='step', color='r', label=r'parallel to outer angular momentum')
# plot.hist (dv_perp_in_orbit, bins, range, histtype='step', color='b', label=r'perperdicular to both')
# ylimits = plot.get_ylim()
# plot.plot ((statistics.mean(dv_parallel), statistics.mean(dv_parallel)), ylimits, 'k--')
# plot.plot ((statistics.mean(dv_perp_to_orbit), statistics.mean(dv_perp_to_orbit)), ylimits, 'r--')
# plot.plot ((statistics.mean(dv_perp_in_orbit), statistics.mean(dv_perp_in_orbit)), ylimits, 'b--')
# plot.set_ylim(ylimits)
# plot.legend()

# pyplot.tight_layout()
# pyplot.savefig(filepath[:-4]+'_noweak.pdf')
# pyplot.clf()

print('Total velocity change:', dv_total)
# print('Average velocity change along v:', statistics.mean(dv_parallel))
# print('Average velocity change perpendicular to v:', statistics.mean(dv_perp_in_orbit), ',', statistics.mean(dv_perp_to_orbit))
