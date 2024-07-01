from orbital.utilities import true_anomaly_from_mean
from math import copysign
import numpy as np
import random
import pickle
import sys
import time
import fortran.Mikkola
import astropy.units as u
from astropy import constants
from core import integrate1, orbital_elements_from_nbody, orbital_elements_to_orbital_vectors, orbital_vectors_to_cartesian, compute_eps_SA, compute_eps_oct
G = constants.G

def plot_track(x,y):
	from matplotlib import pyplot
	figure = pyplot.figure(figsize=(10, 10))
	pyplot.rcParams.update({'font.size': 30})
	plot = figure.add_subplot(1,1,1)
	ax = pyplot.gca()
	ax.minorticks_on() 
	ax.locator_params(nbins=3)

	x_label = 'x [au]'
	y_label = 'y [au]'
	pyplot.xlabel(x_label)
	pyplot.ylabel(y_label)

	plot.scatter([0.0], [0.0], color='y', lw=8)
	plot.plot(x.to(u.AU), y.to(u.AU), color = 'b')
	plot.set_xlim(-20, 20)
	plot.set_ylim(-20, 20)

	save_file = 'scattering.png'
	pyplot.savefig(save_file)
	print('\nSaved figure in file', save_file,'\n')
	pyplot.show()

def plot_track_comparison(x,y,maxXY,save_file):
	x_old = []
	y_old = []
	with open('old_code.txt') as f:
		for line in f:
			data = line.split()
			x_old.append(float(data[0]))
			y_old.append(float(data[1]))

	from matplotlib import pyplot
	figure = pyplot.figure(figsize=(12, 12))
	pyplot.rcParams.update({'font.size': 30})
	plot = figure.add_subplot(1,1,1)
	ax = pyplot.gca()
	ax.minorticks_on() 
	ax.locator_params(nbins=3)

	x_label = 'x / a'
	y_label = 'y / a'
	pyplot.xlabel(x_label)
	pyplot.ylabel(y_label)

	plot.scatter([0.0], [0.0], color='y', lw=8)
	new_code = plot.plot(x.to(u.AU), y.to(u.AU), color = 'b', label = 'new code')
	old_code = plot.plot(x_old, y_old, 'r--', label = 'old code')
	plot.legend()
	plot.set_xlim(-maxXY, maxXY)
	plot.set_ylim(-maxXY, maxXY)

	pyplot.savefig(save_file)
	print('\nSaved figure in file', save_file,'\n')
	pyplot.show()

def plot_distances (save_file):

	with open('distances-new', 'rb') as file:
		tr1r2 = pickle.load(file)
	t_new = tr1r2[0]
	r1_new = tr1r2[1]
	r2_new = tr1r2[2]

	t_old = []
	r1_old = []
	r2_old = []
	with open('distances-old.txt') as f:
		for line in f:
			data = line.split()
			t_old.append(float(data[0]))
			r1_old.append(float(data[1]))
			r2_old.append(float(data[2]))	

	from matplotlib import pyplot
	figure = pyplot.figure(figsize=(12, 12))
	pyplot.rcParams.update({'font.size': 30})
	plot = figure.add_subplot(1,1,1,yscale="log")
	ax = pyplot.gca()
	ax.minorticks_on() 
	# ax.locator_params(nbins=3)
	pyplot.xlabel('t')
	pyplot.ylabel('r')

	plot.plot(t_new, r1_new, 'k')
	plot.plot(t_old, r1_old, 'k--')
	plot.plot(t_new, r2_new, 'r')
	plot.plot(t_old, r2_old, 'r--')
	plot.set_xlim(250, 400)
	# plot.set_ylim(-maxXY, maxXY)

	pyplot.savefig(save_file)
	print('\nSaved figure in file', save_file,'\n')
	pyplot.show()

# outliers, all angles
# compare_new_r3maxNbody=100_m1=5_m2=5_a=1_e=0.1_m3=5_v=3.txt
def makePlot(fileName1, e, save_file):
	Qamin = 5
	Qamax = 20
	N = 10
	Qa = [Qamin+n/N*(Qamax-Qamin) for n in range(N+1)]
	# print(Qa)
	totalNumber = np.array([0]*(N+1))
	outlierNumber_hybrid_Omega = np.array([0]*(N+1))
	outlierNumber_hybrid_omega = np.array([0]*(N+1))
	outlierNumber_sa_Omega = np.array([0]*(N+1))
	outlierNumber_sa_omega = np.array([0]*(N+1))
	# outlierFile = open(fileName2, 'w+')

	with open(fileName1) as f:
		lineNumber = 0
		for line in f:
			lineNumber+=1
			if lineNumber>1:
				data = line.split()
				QaValue = float(data[0])
				de_3body = float(data[5])
				de_hybrid = float(data[6])
				de_sa = float(data[7])
				di_3body = float(data[8])
				di_hybrid = float(data[9])
				di_sa = float(data[10])
				da_3body = float(data[11])
				da_hybrid = float(data[12])
				da_sa = float(data[13])
				dOmega_3body = float(data[14])
				dOmega_hybrid = float(data[15])
				dOmega_sa = float(data[16])
				domega_3body = float(data[17])
				domega_hybrid = float(data[18])
				domega_sa = float(data[19])
				n = round(N*(QaValue-Qamin)/(Qamax-Qamin))
				totalNumber[n] += 1
				# if float(data[-1])>0: print('yes!')
				if np.abs(dOmega_hybrid/dOmega_3body-1) > 0.2: outlierNumber_hybrid_Omega[n] += 1
				if np.abs(domega_hybrid/domega_3body-1) > 0.2: 
					outlierNumber_hybrid_omega[n] += 1
					# print(line, file=outlierFile)
				if np.abs(dOmega_sa/dOmega_3body-1) > 0.2: outlierNumber_sa_Omega[n] += 1
				if np.abs(domega_sa/domega_3body-1) > 0.2: outlierNumber_sa_omega[n] += 1

	# outlierFile.close()					
	import matplotlib
	from matplotlib import pyplot
	from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
	matplotlib.rcParams['text.usetex'] = True
	matplotlib.rcParams['mathtext.fontset'] = 'stix'
	matplotlib.rcParams['font.family'] = 'STIXGeneral'
	matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
	figure = pyplot.figure() #(figsize=(12, 12))
	plot = figure.add_subplot(1,1,1)
	ax = pyplot.gca()
	ax.minorticks_on() 
	# ax.xaxis.set_major_locator(MultipleLocator(1))
	# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
	# ax.yaxis.set_major_locator(MultipleLocator(0.1))
	# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
	ax.tick_params(labelsize=14)
	# ax.set_xlabel(r'$r_3$', fontsize=16)
	# ax.set_ylabel(r'$\delta e$', fontsize=16)
	ax.set_xlabel(r'$Q/a$', fontsize=16)
	ax.set_ylabel(r'fraction of $\Delta\Omega$ or $\Delta\omega$ discrepancies $>20\%$', fontsize=16)
	# pyplot.xscale('log')
	# pyplot.yscale('log')
	pyplot.text(0.3, 0.9, '$e='+str(e)+'$', fontsize=16, transform=ax.transAxes)
	# plot.plot(x40, y40, 'k', label=r'3-body, $r_{3,\mathrm{max}}=40$')
	plot.plot(Qa, outlierNumber_hybrid_Omega/totalNumber, 'k', label=r'$\Omega$, hybrid')
	plot.plot(Qa, outlierNumber_sa_Omega/totalNumber, 'r', label=r'$\Omega$, SA')
	plot.plot(Qa, outlierNumber_hybrid_omega/totalNumber, 'k--', label=r'$\omega$, hybrid')
	plot.plot(Qa, outlierNumber_sa_omega/totalNumber, 'r--', label=r'$\omega$, SA')
	# plot.plot(x2, y2, 'r', label='SA')
	# plot.plot(x2neg, y2neg, 'r--')
	# plot.plot(x2, y2, 'k--', label=r'3-body, $\delta t=100$')
	# plot.plot(x_sa, y_sa, 'r', label=r'SA')
	plot.legend(fontsize=16, frameon=False)
	# plot.set_xlim(500, 5000)
	# plot.set_ylim(1e-8, 1e-5)
	pyplot.tight_layout()
	pyplot.savefig(save_file)


def analyze_comparison_results (fileName):
	ejected0_3body = 0
	ejected1_3body = 0
	ejected2_3body = 0
	destroyed_3body = 0
	bound_3body = 0
	timeout_3body = 0
	ejected0_hybrid = 0
	ejected1_hybrid = 0
	ejected2_hybrid = 0
	destroyed_hybrid = 0
	bound_hybrid = 0
	timeout_hybrid = 0

	with open(fileName) as f:
		lineNumber = 0
		for line in f:
			lineNumber+=1
			if lineNumber>1:
				data = line.split()
				result_3body = float(data[4])
				result_hybrid = float(data[5])
				third_body_3body = float(data[6])
				third_body_hybrid = float(data[7])
				if result_3body==0 and third_body_3body==0: ejected0_3body+=1
				if result_3body==0 and third_body_3body==1: ejected1_3body+=1
				if result_3body==0 and third_body_3body==2: ejected2_3body+=1
				if result_3body==1: bound_3body+=1
				if result_3body==2: destroyed_3body+=1
				if result_3body==3: timeout_3body+=1
				if result_hybrid==0 and third_body_hybrid==0: ejected0_hybrid+=1
				if result_hybrid==0 and third_body_hybrid==1: ejected1_hybrid+=1
				if result_hybrid==0 and third_body_hybrid==2: ejected2_hybrid+=1
				if result_hybrid==1: bound_hybrid+=1
				if result_hybrid==2: destroyed_hybrid+=1
				if result_hybrid==3: timeout_hybrid+=1

	print ('ejected_3body: ', ejected0_3body, ejected1_3body, ejected2_3body)
	print ('ejected_hybrid: ', ejected0_hybrid, ejected1_hybrid, ejected2_hybrid)
	print ('destroyed_3body: ', destroyed_3body)
	print ('destroyed_hybrid: ', destroyed_hybrid)
	print ('bound_3body: ', bound_3body)
	print ('bound_hybrid: ', bound_hybrid)
	print ('timeout_3body: ', timeout_3body)
	print ('timeout_hybrid: ', timeout_hybrid)

def scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=50, dt=1000):
	#returns the binary parameters after the interaction

	G_0 = G.to(u.m**3/u.s**2/u.kg).value
	m12_0 = (m1+m2).to(u.kg).value
	m13_0 = (m1+m3).to(u.kg).value
	m23_0 = (m2+m3).to(u.kg).value
	m123_0 = (m1+m2+m3).to(u.kg).value
	a_0 = a.to(u.m).value
	aStar_0 = aStar.to(u.m).value

	#initialize the binary
	if meanAnomaly>6.25: meanAnomaly=0
	trueAnomaly = true_anomaly_from_mean(e, meanAnomaly)
	ex,ey,ez,jx,jy,jz = orbital_elements_to_orbital_vectors(e, i, omega, Omega)
	r_binary, v_binary = orbital_vectors_to_cartesian(
		G_0, m12_0, a_0, trueAnomaly, ex,ey,ez,jx,jy,jz)
	r_binary *= u.m
	v_binary *= u.m/u.s
	x2 = r_binary * m1 / (m1+m2)
	x1 = x2 - r_binary
	v2 = v_binary * m1 / (m1+m2)
	v1 = v2 - v_binary

	velocityUnit = np.sqrt(G*(m1+m2)/a)

	#initialize the 3rd body
	trueAnomalyStar = -np.arccos((aStar/(r3max*a)*(1-eStar**2)-1)/eStar)
	ex,ey,ez,jx,jy,jz = orbital_elements_to_orbital_vectors(
		eStar, iStar, omegaStar, OmegaStar)
	r_starAndTheBinary, v_starAndTheBinary = orbital_vectors_to_cartesian(
		G_0, m123_0, aStar_0, trueAnomalyStar, ex,ey,ez,jx,jy,jz)
	r_starAndTheBinary *= u.m
	v_starAndTheBinary *= u.m/u.s
	x3 = r_starAndTheBinary * (m1+m2) / (m1+m2+m3)
	x_bin = x3 - r_starAndTheBinary
	v3 = v_starAndTheBinary * (m1+m2) / (m1+m2+m3)
	v_bin = v3 - v_starAndTheBinary

	#put the binary into the CM refernce frame
	x1 += x_bin
	x2 += x_bin
	v1 += v_bin
	v2 += v_bin

	# input parameters for chainevolve
	x = np.array([x1[0]/a, x1[1]/a, x1[2]/a, x2[0]/a, x2[1]/a, x2[2]/a, x3[0]/a, x3[1]/a, x3[2]/a])
	v = np.array([v1[0]/velocityUnit, v1[1]/velocityUnit, v1[2]/velocityUnit, v2[0]/velocityUnit, v2[1]/velocityUnit, v2[2]/velocityUnit, v3[0]/velocityUnit, v3[1]/velocityUnit, v3[2]/velocityUnit])
	m = np.array([m1/(m1+m2), m2/(m1+m2), m3/(m1+m2)])
	time = np.array(0)
	N = 3
	dt = 1000
	eps=1.0e-19
	newreg=1
	ksmx=10000
	soft=0
	cmet=np.array([1,0.001,0])
	clight=0
	ixc=0
	NBH=2
	spin=np.array([0,0,0])
	cmxx = np.array([0,0,0])
	cmvx = np.array([0,0,0])
	nmergers=np.array(0)
	mergers=np.array([0,0,0])

	r1 = []
	r2 = []

	result = 0
	# 0 - 3rd body ejected
	# 1 - 3rd body still bound to the binary
	# 2 - all 3 bodies ejected
	# 3 - calculation abandoned
	third_body = 2
	# 0,1,2 - number of the particle which ends up being the 3rd body (initially 2)
	dv_binary = [0,0,0]*u.m/u.s

	while True:
		fortran.Mikkola.chainevolve(N, x, v, m, time, dt, eps, newreg, ksmx, soft, cmet, clight, ixc, NBH, spin,cmxx, cmvx, mergers, nmergers)
		finalTime = time * a / velocityUnit
		x1 = x[0:3] * a
		x2 = x[3:6] * a
		x3 = x[6:9] * a
		v1 = v[0:3] * velocityUnit
		v2 = v[3:6] * velocityUnit
		v3 = v[6:9] * velocityUnit
		x21_0 = (x2-x1).to(u.m).value
		v21_0 = (v2-v1).to(u.m/u.s).value
		a_12, e_12, i_12, Omega_12, omega_12 = orbital_elements_from_nbody(
			G_0, m12_0, x21_0, v21_0)
		a_12 *= u.m
		rBinaryCM = (x1*m1 + x2*m2) / (m1 + m2)
		vBinaryCM = (v1*m1 + v2*m2) / (m1 + m2)  
		r3 = np.linalg.norm(x3 - rBinaryCM)
		if a_12>0*u.m and r3>r3max*a_12: #third body escaped (maybe not forever) and the binary is still bound
			third_body = 2
			starEnergy = m3*np.linalg.norm(v3)**2/2 + (m1+m2)*np.linalg.norm(vBinaryCM)**2/2 - G*m3*(m1+m2)/r3
			if starEnergy<0*u.kg*u.m**2/u.s**2: 
				result = 1
			dv_binary = vBinaryCM - v_bin
			x3bin_0 = (x3-rBinaryCM).to(u.m).value
			v3bin_0 = (v3-vBinaryCM).to(u.m/u.s).value
			aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new = 			orbital_elements_from_nbody(G_0, m123_0, x3bin_0, v3bin_0)
			aStar_new *= u.m
			return result, third_body, finalTime, dv_binary, a_12, e_12, i_12, Omega_12, omega_12, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3/a_12

		x23_0 = (x2-x3).to(u.m).value
		v23_0 = (v2-v3).to(u.m/u.s).value
		a_32, e_32, i_32, Omega_32, omega_32 = orbital_elements_from_nbody(
			G_0, m23_0, x23_0, v23_0)
		a_32 *= u.m
		rBinaryCM = (x3*m3 + x2*m2) / (m3 + m2)
		vBinaryCM = (v3*m3 + v2*m2) / (m3 + m2)  
		r1 = np.linalg.norm(x1 - rBinaryCM)
		if a_32>0*u.m and r1>r3max*a_32:	#first body escaped
			third_body = 0
			starEnergy = m1*np.linalg.norm(v1)**2/2 + (m3+m2)*np.linalg.norm(vBinaryCM)**2/2 - G*m1*(m3+m2)/r1
			if starEnergy<0*u.kg*u.m**2/u.s**2: 
				result = 1
			dv_binary = vBinaryCM - v_bin
			x1bin_0 = (x1-rBinaryCM).to(u.m).value
			v1bin_0 = (v1-vBinaryCM).to(u.m/u.s).value
			aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new = 				orbital_elements_from_nbody(G_0, m123_0, x1bin_0, v1bin_0)
			aStar_new *= u.m
			return result, third_body, finalTime, dv_binary, a_32, e_32, i_32, Omega_32, omega_32, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r1/a_32	

		x31_0 = (x3-x1).to(u.m).value
		v31_0 = (v3-v1).to(u.m/u.s).value
		a_13, e_13, i_13, Omega_13, omega_13 = orbital_elements_from_nbody(
			G_0, m23_0, x31_0, v31_0)
		a_13 *= u.m
		rBinaryCM = (x1*m1 + x3*m3) / (m1 + m3)
		vBinaryCM = (v1*m1 + v3*m3) / (m1 + m3)  
		r2 = np.linalg.norm(x2 - rBinaryCM)
		if a_13>0*u.m and r2>r3max*a_13:	#second body escaped
			third_body = 1
			starEnergy = m2*np.linalg.norm(v2)**2/2 + (m1+m3)*np.linalg.norm(vBinaryCM)**2/2 - G*m2*(m1+m3)/r2
			if starEnergy<0*u.kg*u.m**2/u.s**2: result = 1
			dv_binary = vBinaryCM - v_bin
			x2bin_0 = (x2-rBinaryCM).to(u.m).value
			v2bin_0 = (v2-vBinaryCM).to(u.m/u.s).value
			aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new = 				orbital_elements_from_nbody(G_0, m123_0, x2bin_0, v2bin_0)
			aStar_new *= u.m
			return result, third_body, finalTime, dv_binary, a_13, e_13, i_13, Omega_13, omega_13, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r2/a_13
		
		if a_12<0*u.m and a_13<0*u.m and a_32<0*u.m and r1>r3max*a and r2>r3max*a and r3>r3max*a:	#binary destroyed
			result = 2
			return result, third_body, 0*u.s, 0*u.m/u.s, a, e, i, Omega, omega, aStar, eStar, iStar, OmegaStar, omegaStar, 0

		if finalTime > 1e10*u.yr:	#time limit exceeded
			result = 3
			return result, third_body, 0*u.s, 0*u.m/u.s, a, e, i, Omega, omega, aStar, eStar, iStar, OmegaStar, omegaStar, 0

def normalize (vector):
	return vector / np.linalg.norm(vector)

class Arguments:
	G = 4.0*np.pi**2
	c = 63239.72638679138
	mxstep = 1000000
	verbose = False
	include_quadrupole_terms = True
	include_octupole_terms = True
	include_1PN_terms = False
	N_steps = 1000
	def __init__(self, m1, m2, a, e, i, Omega, omega, M_per, e_per, Q):
		self.m = m1+m2
		self.m1 = m1
		self.m2 = m2
		self.a = a
		self.e = e
		self.i = i
		self.Omega = Omega
		self.omega = omega
		self.M_per = M_per
		self.e_per = e_per
		self.Q = Q

def scattering_SA (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3_i=-1e4, r3_f=1e4):
	#return the binary parameters after the interaction; uses Hamers' SA code

	G_0 = G.to(u.m**3/u.s**2/u.kg).value
	m_0 = (m1+m2).to(u.kg).value
	m_total_0 = (m1+m2+m3).to(u.kg).value
	a_0 = a.to(u.m).value
	aStar_0 = aStar.to(u.m).value

	#initialize the binary
	trueAnomaly = 0	#initial anomaly doesn't matter
	ex,ey,ez,jx,jy,jz = orbital_elements_to_orbital_vectors(e, i, omega, Omega)
	r_binary,v_binary = orbital_vectors_to_cartesian(
		G_0, m_0, a_0, trueAnomaly, ex, ey, ez, jx, jy, jz)

	#initialize the 3rd body
	trueAnomalyStar = 0
	ex,ey,ez,jx,jy,jz = orbital_elements_to_orbital_vectors(eStar, iStar, omegaStar, OmegaStar)
	r_starAndTheBinary, v_starAndTheBinary = orbital_vectors_to_cartesian(
		G_0, m_total_0, aStar_0, trueAnomalyStar, ex, ey, ez, jx, jy, jz) 

	#find the unit vectors for the star reference frame
	h = np.cross(r_starAndTheBinary,v_starAndTheBinary)
	mu = G_0 * m_total_0
	eccVector = np.cross(v_starAndTheBinary,h)/mu - normalize(r_starAndTheBinary)
	e1 = normalize(eccVector)
	e3 = normalize(h)
	e2 = np.cross(e3,e1)
	#coordinate transformation matrix, from cluster to stellar orbit
	A = [e1, e2, e3]
	#binary coordinates in the star reference frame
	r_binary_0 = np.matmul(A, r_binary)
	v_binary_0 = np.matmul(A, v_binary)

	#binary orbital elements in the star reference frame
	a_0, e_0, i_0, Omega_0, omega_0 = orbital_elements_from_nbody(
		G_0, m_0, r_binary_0, v_binary_0)
	Omega_0 -= np.pi
	Q = -aStar*(eStar-1)

	args = Arguments(m1.to(u.solMass).value, m2.to(u.solMass).value, a_0, e_0, i_0, Omega_0, omega_0, m3.to(u.solMass).value, eStar, Q.to(u.m).value)
	args.theta_i = copysign(1,r3_i)*np.arccos((aStar/(np.abs(r3_i)*a)*(1-eStar**2)-1)/eStar)
	args.theta_f = copysign(1,r3_f)*np.arccos((aStar/(np.abs(r3_f)*a)*(1-eStar**2)-1)/eStar)
	if args.theta_f<args.theta_i: 
		args.theta_f+=2*np.pi 
	a_fin_0, e_fin_0, i_fin_0, Omega_fin_0, omega_fin_0 = integrate1(args) 
	Omega_fin_0 += np.pi	

	# calculate the change in the binary CoM velocity
	ex,ey,ez,jx,jy,jz = orbital_elements_to_orbital_vectors(
		eStar, iStar, omegaStar, OmegaStar)
	r_before, v_before = orbital_vectors_to_cartesian(
		G_0, m_total_0, aStar_0, args.theta_i, ex, ey, ez, jx, jy, jz) 
	v_binary_before = -v_before * m3 / (m1+m2+m3)	#assumes r=r_2-r_1, v=v_2-v_1

	ex,ey,ez,jx,jy,jz = orbital_elements_to_orbital_vectors(
		eStar, iStar, omegaStar, OmegaStar)
	r_after, v_after = orbital_vectors_to_cartesian(
		G_0, m_total_0, aStar_0, args.theta_f, ex, ey, ez, jx, jy, jz) 
	v_binary_after = -v_after * m3 / (m1+m2+m3)	#assumes r=r_2-r_1, v=v_2-v_1

	dv_binary = v_binary_after - v_binary_before

	ex,ey,ez,jx,jy,jz = orbital_elements_to_orbital_vectors(e_fin_0,i_fin_0,omega_fin_0,Omega_fin_0)
	theta_bin = 0
	r_fin_0,v_fin_0 = orbital_vectors_to_cartesian(G_0,m_0,a_fin_0,theta_bin,ex,ey,ez,jx,jy,jz)

	A_inv = np.linalg.inv(A)
	r_fin = np.matmul(A_inv, r_fin_0)
	v_fin = np.matmul(A_inv, v_fin_0)

	a_fin, e_fin, i_fin, Omega_fin, omega_fin = orbital_elements_from_nbody(
		G_0, m_0, r_fin, v_fin)

	return dv_binary*u.m/u.s, a_fin*u.m, e_fin, i_fin, Omega_fin, omega_fin

def scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=10000, r3max_nbody=100, meanAnomaly0=-1, n_orbits_max=20):
	# result:
	# 0 - 3rd body ejected
	# 1 - calculation abandoned after more than n_orbits_max bound orbits
	# 2 - all 3 bodies ejected
	# 3 - calculation abandoned after spending too much time in a 3-body phase

	# fileName='/home/alexander/Dropbox/code/flybys-master/debug.txt'
	# file = open(fileName, "a")
	# print('---', file=file)
	# print(aStar, eStar, iStar, OmegaStar, omegaStar, r3max_nbody, file=file)
	binary1 = 0
	binary2 = 1
	dv_binary_final = [0,0,0]*u.m/u.s
	third_body_final = 2
	dv_binary, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3_i=-r3max, r3_f=-r3max_nbody)
	print(dv_binary)
	dv_binary_final += dv_binary
	# print(a_fin_sa/a)
	random.seed()
	if meanAnomaly0<0:
		meanAnomaly = random.random()*2*np.pi
		# print(meanAnomaly, file=file)
		# print(meanAnomaly)
	else:
		meanAnomaly=meanAnomaly0
	result, third_body, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_nbody_new  = scattering (m1, m2, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max_nbody, dt=100)
	dv_binary_final += dv_binary
	if result==2 or result==3:
		return result, third_body_final, dv_binary, a, e, i, Omega, omega, 0
	if third_body==0:
		x = binary1
		binary1 = third_body_final
		third_body_final = x
		m = m1
		m1 = m3
		m3 = m
	if third_body==1:
		x = binary2
		binary2 = third_body_final
		third_body_final = x
		m = m2
		m2 = m3
		m3 = m
	# print(result, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_nbody_new, file=file)
	n_orbits = 0
	while eStar_new<1 and n_orbits<n_orbits_max:
		# print(n_orbits, file=file)
		dv_binary, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a_new, e_new, i_new, Omega_new, omega_new, m3, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3_i=r3max_nbody_new, r3_f=-r3max_nbody)
		dv_binary_final += dv_binary
		# print(a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa, file=file)
		meanAnomaly = random.random()*2*np.pi
		# print(meanAnomaly, file=file)
		# print(meanAnomaly)
		result, third_body, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_nbody_new  = scattering (m1, m2, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa, meanAnomaly, m3, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3max=r3max_nbody, dt=100)
		dv_binary_final += dv_binary
		if result==2 or result==3:
			return result, third_body_final, dv_binary, a, e, i, Omega, omega, 0
		if third_body==0:
			x = binary1
			binary1 = third_body_final
			third_body_final = x
			m = m1
			m1 = m3
			m3 = m
		if third_body==1:
			x = binary2
			binary2 = third_body_final
			third_body_final = x
			m = m2
			m2 = m3
			m3 = m
		# print(result, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_nbody_new, file=file)
		n_orbits += 1
	# if n_orbits<n_orbits_max:
	if eStar_new>=1:
		dv_binary, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a_new, e_new, i_new, Omega_new, omega_new, m3, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3_i=r3max_nbody_new, r3_f=r3max)
		dv_binary_final += dv_binary
	else:
		result = 1
		# print(a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa, file=file)
	# file.close()
	return result, third_body_final, dv_binary_final, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa, n_orbits

def scattering_SA_bound_test (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3_i=-1e4, r3_f=1e4):
	from core import integrate1, orbital_elements_from_nbody, orbital_elements_to_orbital_vectors, orbital_vectors_to_cartesian, compute_eps_SA, compute_eps_oct
	aStar = 5000*a
	eStar = 0.999
	Q = -aStar*(eStar-1)

	args = Arguments(m1.to(u.solMass), m2.to(u.solMass), a.to(u.m), e, i, Omega, omega, m3.to(u.solMass), eStar, Q.to(u.m))
	args.fraction_theta_i = 0.5
	args.fraction_theta_f = 1.5
	a_fin_0, e_fin_0, i_fin_0, Omega_fin_0, omega_fin_0 = integrate1(args) 

	return a_fin_0, e_fin_0, i_fin_0, Omega_fin_0, omega_fin_0

# def scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3max=5000, r3max_nbody=300):
# 	a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3_i=-r3max, r3_f=-r3max_nbody)
# 	random.seed()
# 	meanAnomaly = random.random()*2*np.pi
# 	result, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, v_new, p_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_new  = scattering (m1, m2, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max_nbody, dt=100)
# 	if result!=0:
# 		a_fin_sa = a 
# 		e_fin_sa = e
# 		i_fin_sa = i
# 		Omega_fin_sa = Omega
# 		omega_fin_sa = omega
# 	else:
# 		a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a_new, e_new, i_new, Omega_new, omega_new, m3, v_new, p_new, iStar_new, OmegaStar_new, omegaStar_new, r3_i=r3max_new, r3_f=r3max)
# 	return a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa

def plotHistogram (save_file):
	dE_old = []
	with open('dE-old.txt') as f:
		for line in f:
			data = line.split()
			dE_old.append(float(data[0]))
	dE_new = []
	with open('dE-new.txt') as f:
		for line in f:
			data = line.split()
			dE_new.append(float(data[0]))
	from matplotlib import pyplot
	pyplot.hist(dE_old, 30, histtype = 'step')
	pyplot.hist(dE_new, 30, histtype = 'step')
	pyplot.xlabel('dE')
	pyplot.savefig(save_file)
	print('\nSaved figure in file', save_file,'\n')
	return

# def compare (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, iStar, OmegaStar, omegaStar, Qmin=5, Qmax=20, N=25, N_3body = 300, r3max=400):
# 	# compare 3-body with SA

# 	Q_array = np.logspace(np.log10(Qmin), np.log10(Qmax), N) * a
# 	aStar = -G*(m1+m2+m3)/v**2
# 	fileName = '/home/alexander/flybys-master/comparison_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_i='+str(i)+'_Omega='+str(Omega)+'_omega='+str(omega)+'_meanAnomaly='+str(meanAnomaly)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'_iStar='+str(iStar)+'_OmegaStar='+str(OmegaStar)+'_omegaStar='+str(omegaStar)+'.txt'
# 	file = open(fileName, "w+")
# 	for Q in Q_array:
# 		eStar = 1 - Q/aStar
# 		p = -aStar*np.sqrt(eStar**2-1)
# 		result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
# 		a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar)
# 		file.write(str(Q/a)+" "+str(i_fin_3body-i)+" "+str(i_fin_sa-i)+"\n")
# 	file.close()
# 	return

def destroyedFraction (m1, m2, a, e, m3, v, Qamax=5, N=25, N_3body = 100, r3max=100):
	aStar = -G*(m1+m2+m3)/v**2
	i = 0
	Omega = 0
	omega = 0
	fileName = '/home/alexander/Dropbox/code/flybys-master/destroyed_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'.txt'
	file = open(fileName, "w+")
	file.write('Q/a destroyed_fraction exchange_fraction abandoned_fraction\n')
	random.seed()
	for n in range(1, N):
		Q = n/N*Qamax*a
		N_destroyed = 0
		N_abandoned = 0
		N_exchange = 0
		eStar = 1 - Q/aStar
		p = -aStar*np.sqrt(eStar**2-1)
		for n1 in range(N_3body):
			meanAnomaly = random.random()*2*np.pi
			iStar = np.arccos(random.random()*2-1)
			OmegaStar = random.random()*2*np.pi
			omegaStar = random.random()*2*np.pi
			result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
			if result == 2: N_destroyed+=1
			if result == 1: N_abandoned+=1
			if result == 3: N_exchange+=1
		file.write(str(Q/a)+' '+str(N_destroyed/N_3body)+' '+str(N_exchange/N_3body)+' '+str(N_abandoned/N_3body)+'\n')
	  
def compare (m1, m2, a, e, m3, v, Qamin=5, Qamax=20, N=10, N_3body = 100):
	aStar = -G*(m1+m2+m3)/v**2
	i = 0
	Omega = 0
	omega = 0
	r3maxMin = 1000
	r3maxMax = 10000
	fileName = '/home/alexander/Dropbox/code/flybys-master/compare_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'_N_3body+'+'.txt'
	outliersFileName = '/home/alexander/Dropbox/code/flybys-master/outliers_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'.txt'
	file = open(fileName, "a")
	outliers = open(outliersFileName, "a")
	file.write('Q/a de_3body de_sa\n')
	outliers.write('p, iStar, OmegaStar, omegaStar, meanAnomaly\n')
	random.seed()
	for n in range(0, N):
		Q = (Qamin+n/N*(Qamax-Qamin))*a
		eStar = 1 - Q/aStar
		p = -aStar*np.sqrt(eStar**2-1)
		for n1 in range(N_3body):
			meanAnomaly = random.random()*2*np.pi
			iStar = np.arccos(random.random()*2-1)
			OmegaStar = random.random()*2*np.pi
			omegaStar = random.random()*2*np.pi
			plateau = False
			r3max = r3maxMin
			while (not plateau) and r3max<r3maxMax:
				a_fin_sa, e_fin_sa_1, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
				r3max *= 2
				a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
				if np.abs((e_fin_sa-e)/(e_fin_sa_1-e) - 1) < 0.2: plateau = True
			result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
			if result!=0: print('result = ', result, '\n') 
			file.write(str(Q/a)+' '+str(e_fin_3body-e)+' '+str(e_fin_sa-e)+'\n')
			if (not plateau) or np.abs((e_fin_sa-e)/(e_fin_3body-e) - 1) > 0.2: 
				print(plateau, r3max, p, iStar, OmegaStar, omegaStar, meanAnomaly, file=outliers)

def compareHybrid (m1, m2, a, e, m3, v, Qamin=5, Qamax=20, N=10, N_3body = 1000, r3max=5000, n_orbits_max=20):
	aStar = -G*(m1+m2+m3)/v**2
	i = 0
	Omega = 0
	omega = 0
	r3maxNbodyMin = 100
	r3maxNbodyMax = 1000
	fileName = '/home/alexander/Dropbox/code/flybys-master/compareHybrid_more_r3max='+str(r3max)+'_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'.txt'
	file = open(fileName, "a")
	file.write('Q/a plateau r3maxNbody eStar iStar OmegaStar omegaStar meanAnomaly de_3body de_hybrid de_sa di_3body di_hybrid di_sa da_3body da_hybrid da_sa n_orbits\n')
	file.close()
	random.seed()
	for n in range(0, N):
		Q = (Qamin+n/N*(Qamax-Qamin))*a
		eStar = 1 - Q/aStar
		# p = -aStar*np.sqrt(eStar**2-1)
		for n1 in range(N_3body):
			meanAnomaly = random.random()*2*np.pi
			iStar = np.arccos(random.random()*2-1) 
			OmegaStar = random.random()*2*np.pi
			omegaStar = random.random()*2*np.pi
			file = open(fileName, "a")
			print(Q/a, eStar, iStar, OmegaStar, omegaStar, meanAnomaly, file=file)
			file.close()
			plateau = False 
			r3maxNbody = r3maxNbodyMin
			while (not plateau) and r3maxNbody<r3maxNbodyMax:
				a_fin_hybrid, e_fin_hybrid_1, i_fin_hybrid, Omega_fin_hybrid, omega_fin_hybrid, n_orbits = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=r3maxNbody, n_orbits_max=n_orbits_max)
				r3maxNbody *= 2
				a_fin_hybrid, e_fin_hybrid, i_fin_hybrid, Omega_fin_hybrid, omega_fin_hybrid, n_orbits = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=r3maxNbody, n_orbits_max=n_orbits_max)
				if np.abs((e_fin_hybrid-e)/(e_fin_hybrid_1-e) - 1) < 0.2: plateau = True
				if n_orbits==n_orbits_max:
					file = open(fileName, "a")
					print("more than "+str(n_orbits)+" orbits", r3maxNbody, eStar, iStar, OmegaStar, omegaStar, meanAnomaly, file=file)
					file.close()
					break
			if n_orbits<n_orbits_max:
				result, finalTime, dv_binary, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max)
				a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3_i=-r3max, r3_f=r3max)
				file = open(fileName, "a")
				print(plateau, r3maxNbody, e_fin_3body-e, e_fin_hybrid-e, e_fin_sa-e, i_fin_3body-i, i_fin_hybrid-i, i_fin_sa-i, a_fin_3body/a-1, a_fin_hybrid/a-1, a_fin_sa/a-1, n_orbits, file=file)
				file.close()

def compare_new (m1, m2, a, e, m3, v, r3maxNbody = 50, Qamin=5, Qamax=20, dQa=1.5, N_3body = 1000, r3max=5000, n_orbits_max=20):
	aStar = -G*(m1+m2+m3)/v**2
	i = 0
	Omega = 0
	omega = 0
	fileName = '/Users/axr6631/Dropbox/code/flybys-master/compare_new_r3maxNbody='+str(r3maxNbody)+'_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'.txt'
	file = open(fileName, "a")
	# file.write('Q/a iStar OmegaStar omegaStar de_3body de_hybrid de_sa di_3body di_hybrid di_sa da_3body da_hybrid da_sa dOmega_3body dOmega_hybrid dOmega_sa domega_3body domega_hybrid domega_sa n_orbits\n')
	file.write('Q/a iStar OmegaStar omegaStar result third_body de_3body de_hybrid di_3body di_hybrid da_3body da_hybrid dOmega_3body dOmega_hybrid domega_3body domega_hybrid n_orbits\n')
	file.close()
	random.seed()
	Qa = Qamin
	while Qa <= Qamax:
		Q = Qa*a
		eStar = 1 - Q/aStar
		for n1 in range(N_3body):
			meanAnomaly = random.random()*2*np.pi
			iStar = np.arccos(random.random()*2-1) 
			OmegaStar = random.random()*2*np.pi
			omegaStar = random.random()*2*np.pi
			# file = open(fileName, "a")
			# print(Q/a, eStar, iStar, OmegaStar, omegaStar, meanAnomaly, file=file)
			# file.close()
			result_hybrid, third_body_hybrid, dv_binary_hybrid, a_fin_hybrid, e_fin_hybrid, i_fin_hybrid, Omega_fin_hybrid, omega_fin_hybrid, n_orbits = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=r3maxNbody, n_orbits_max=n_orbits_max)
			# if n_orbits==n_orbits_max:
			# 	file = open(fileName, "a")
			# 	print("more than "+str(n_orbits)+" orbits", eStar, iStar, OmegaStar, omegaStar, meanAnomaly, file=file)
			# 	file.close()
			# else:
			result_3body, third_body_3body, finalTime, dv_binary_3body, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max)
			# a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3_i=-r3max, r3_f=r3max)
			file = open(fileName, "a")
			# print(Q/a, iStar, OmegaStar, omegaStar, result_3body, result_hybrid, third_body_3body, third_body_hybrid, dv_binary_3body, dv_binary_hybrid, e_fin_3body-e, e_fin_hybrid-e, (i_fin_3body-i).to(u.rad), (i_fin_hybrid-i).to(u.rad), a_fin_3body/a-1, a_fin_hybrid/a-1, (Omega_fin_3body-Omega).to(u.rad), (Omega_fin_hybrid-Omega).to(u.rad), (omega_fin_3body-omega).to(u.rad), (omega_fin_hybrid-omega).to(u.rad), n_orbits, file=file)
			print(Q/a, iStar, OmegaStar, omegaStar, result_3body, result_hybrid, third_body_3body, third_body_hybrid, dv_binary_3body, dv_binary_hybrid, e_fin_3body-e, e_fin_hybrid-e, i_fin_3body-i, i_fin_hybrid-i, a_fin_3body/a-1, a_fin_hybrid/a-1, Omega_fin_3body-Omega, Omega_fin_hybrid-Omega, omega_fin_3body-omega, omega_fin_hybrid-omega, n_orbits, file=file)
			file.close()
		Qa += dQa

def compareHybridOnly (m1, m2, a, e, m3, v, r3max=5000, n_orbits_max=20):
	aStar = -G*(m1+m2+m3)/v**2
	i = 0
	Omega = 0
	omega = 0
	# homeDir = '/Users/axr6631'
	homeDir = '/home/alexander'
	nbodyFileName = homeDir+'/Dropbox/code/flybys-master/inaccurate_hybrid_code/compareHybrid_more_r3max='+str(r3max)+'_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'.txt'
	fileName = homeDir+'/Dropbox/code/flybys-master/compareHybridOnly_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'.txt'
	file = open(fileName, "a")
	file.write('Q/a eStar iStar OmegaStar omegaStar de_50 de_100 de_200 di_50 di_100 di_200 da_50 da_100 da_200 n_orbits_50 n_orbits_100 n_orbits_200\n')
	file.close()
	random.seed()
	with open(nbodyFileName) as f:
		lineNumber = 0
		for line in f:
			lineNumber+=1
			if lineNumber % 2 == 0 and lineNumber>3299:
				data = line.split()
				Qa = float(data[0])
				eStar = float(data[1])
				iStar = float(data[2])
				OmegaStar = float(data[3])
				omegaStar = float(data[4])
				file = open(fileName, "a")
				print(Qa, eStar, iStar, OmegaStar, omegaStar, file=file)
				file.close()
				a_fin_hybrid_50, e_fin_hybrid_50, i_fin_hybrid_50, Omega_fin_hybrid_50, omega_fin_hybrid_50, n_orbits_50 = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=50, n_orbits_max=n_orbits_max)
				a_fin_hybrid_100, e_fin_hybrid_100, i_fin_hybrid_100, Omega_fin_hybrid_100, omega_fin_hybrid_100, n_orbits_100 = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=100, n_orbits_max=n_orbits_max)
				a_fin_hybrid_200, e_fin_hybrid_200, i_fin_hybrid_200, Omega_fin_hybrid_200, omega_fin_hybrid_200, n_orbits_200 = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=200, n_orbits_max=n_orbits_max)
				file = open(fileName, "a")
				print(e_fin_hybrid_50-e, e_fin_hybrid_100-e, e_fin_hybrid_200-e, i_fin_hybrid_50-i, i_fin_hybrid_100-i, i_fin_hybrid_200-i, a_fin_hybrid_50/a-1, a_fin_hybrid_100/a-1, a_fin_hybrid_200/a-1, n_orbits_50, n_orbits_100, n_orbits_200, file=file)
				file.close()

def r3maxDependence (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max1=100, r3max2=5000, N=100):
	file = open('/home/alexander/Dropbox/code/flybys-master/r3maxDependence1_dt=1000.txt', "a")
	for r3max in np.logspace(np.log10(r3max1), np.log10(r3max2), N):
		result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max, dt=100)
		print(r3max, e_fin_3body-e, file=file)

def r3maxDependence_sa (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3max1=100, r3max2=5000, N=100):
	file = open('/home/alexander/Dropbox/code/flybys-master/r3maxDependence4_sa.txt', "w+")
	for r3max in np.logspace(np.log10(r3max1), np.log10(r3max2), N):
		a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
		print(r3max, e_fin_sa-e, file=file)

# def meanAnomalyDependence (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, N=100, r3max=100):
# 	file = open('/home/alexander/Dropbox/code/flybys-master/meanAnomalyDependence.txt', "w+")
# 	for meanAnomaly in np.linspace(0, 2*np.pi, N):
# 		result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
# 		print(meanAnomaly, e_fin_3body-e)#, file=file)

def QDependence (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, iStar, OmegaStar, omegaStar, Qamin=5, Qamax=20, N=100, r3max=5000):
	aStar = -G*(m1+m2+m3)/v**2
	file = open('/home/alexander/Dropbox/code/flybys-master/QDependence1_r3max='+str(r3max)+'.txt', "a")
	print('Q/a de_3body de_SA de_hybrid', file=file)
	# print('Q/a de_hybrid', file=file)
	for n in range(0, N):
		Q = (Qamin+n/N*(Qamax-Qamin))*a
		eStar = 1 - Q/aStar
		p = -aStar*np.sqrt(eStar**2-1)
		result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body, v_new, p_new, iStar_new, OmegaStar_new, omegaStar_new, r3_new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
		a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3_i=-r3max, r3_f=r3max)
		a_fin_h, e_fin_h, i_fin_h, Omega_fin_h, omega_fin_h = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
		print(Q/a, e_fin_3body-e, e_fin_sa-e, e_fin_h-e, file=file)

def meanAnomalyDependence (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, N=30, r3max=5000, r3max_nbody=200, dt=100):
	# aStar = -G*(m1+m2+m3)/v**2
	# eStar = np.sqrt(1+p**2/aStar**2)
	# file = open('/home/alexander/Dropbox/code/flybys-master/meanAnomalyDependence1_Q=30_r3max='+str(r3max)+'_r3max_nbody='+str(r3max_nbody)+'.txt', "w+")
	file = open('/home/alexander/Dropbox/code/flybys-master/meanAnomalyDependence1_3body_Q=30_r3max='+str(r3max)+'_dt='+str(dt)+'.txt', "w+")
	# print('meanAnomaly de_3body de_hybrid di_3body di_hybrid da_3body da_hybrid', file=file)
	print('meanAnomaly de_3body', file=file)
	for n in range(0, N):
		meanAnomaly = n/N*2*np.pi
		result, finalTime, dv_binary, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3_new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, dt=dt)
		print(meanAnomaly, finalTime, e_fin_3body-e, file=file)
		# a_fin_h, e_fin_h, i_fin_h, Omega_fin_h, omega_fin_h, n_orbits = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=r3max_nbody, meanAnomaly0=meanAnomaly)
		# print(meanAnomaly, e_fin_3body-e, e_fin_h-e, i_fin_3body-i, i_fin_h-i, a_fin_3body/a-1, a_fin_h/a-1, file=file)
		# print(meanAnomaly, e_fin_h-e, file=file)
		# print(meanAnomaly, e_fin_3body-e, file=file)
	file.close()

def QDependenceHybrid (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, iStar, OmegaStar, omegaStar, Qamin=5, Qamax=20, N=100, r3max=15000, r3max_nbody=1000):
	aStar = -G*(m1+m2+m3)/v**2
	file = open('/home/alexander/Dropbox/code/flybys-master/QDependence1_hybrid_r3max='+str(r3max)+'_r3max_nbody='+str(r3max_nbody)+'.txt', "w+")
	print('Q/a de_hybrid', file=file)
	for n in range(0, N):
		Q = (Qamin+n/N*(Qamax-Qamin))*a
		eStar = 1 - Q/aStar
		p = -aStar*np.sqrt(eStar**2-1)
		a_fin_h, e_fin_h, i_fin_h, Omega_fin_h, omega_fin_h, n_orbits = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=r3max_nbody)
		print(Q/a, e_fin_h-e, file=file)

# def QDependenceHamers (eStar, m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, iStar, OmegaStar, omegaStar, Qamin=5, Qamax=20, N=100, r3max=5000):
# 	file = open('/home/alexander/Dropbox/code/flybys-master/QDependence4_eStar=1.01068346206_r3max='+str(r3max)+'.txt', "w+")
# 	print('Q/a de_3body de_sa', file=file)
# 	for n in range(0, N):
# 		Q = (Qamin+n/N*(Qamax-Qamin))*a
# 		aStar = Q/(1-eStar)
# 		p = -aStar*np.sqrt(eStar**2-1)
# 		v = np.sqrt(-G*(m1+m2+m3)/aStar)
# 		result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
# 		a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3_i=-r3max, r3_f=r3max)
# 		print(Q/a, e_fin_3body-e, e_fin_sa-e, file=file)

def makeDatabase (e):
	file3bodyName = '/Users/axr6631/Dropbox/code/flybys-master/inaccurate_hybrid_code/compareHybrid_more_r3max=5000_m1=5_m2=5_a=1_e='+str(e)+'_m3=5_v=3.txt'
	fileHybridName = '/Users/axr6631/Dropbox/code/flybys-master/compareHybridOnly_m1=5_m2=5_a=1_e='+str(e)+'_m3=5_v=3.txt'
	fileResult = open('/Users/axr6631/Dropbox/code/flybys-master/compare_m1=5_m2=5_a=1_e='+str(e)+'_m3=5_v=3.txt', 'w+')
	print('Q/a eStar iStar OmegaStar omegaStar de_3body de_sa de_50 de_100 de_200 di_3body di_sa di_50 di_100 di_200 da_3body da_sa da_50 da_100 da_200 n_orbits_50 n_orbits_100 n_orbits_200', file=fileResult)
	Qa=[]
	eStar=[]
	iStar=[]
	OmegaStar=[]
	omegaStar=[]
	de_3body=[]
	de_sa=[]
	de_50=[]
	de_100=[]
	de_200=[]
	di_3body=[]
	di_sa=[]
	di_50=[]
	di_100=[]
	di_200=[]
	da_3body=[]
	da_sa=[]
	da_50=[]
	da_100=[]
	da_200=[]
	n_orbits_50=[]
	n_orbits_100=[]
	n_orbits_200=[]	

	lineNumber = 0
	with open(file3bodyName) as f:
		for line in f:
			lineNumber+=1
			data = line.split()
			if lineNumber % 2 == 0:
				Qa.append(float(data[0]))
				eStar.append(float(data[1]))
				iStar.append(float(data[2]))
				OmegaStar.append(float(data[3]))
				omegaStar.append(float(data[4]))
			if lineNumber % 2 == 1 and lineNumber>1:
				de_3body.append(float(data[2]))
				de_sa.append(float(data[4]))
				di_3body.append(float(data[5]))
				di_sa.append(float(data[9]))
				da_3body.append(float(data[11]))
				da_sa.append(float(data[13]))
	lineNumber = 0
	with open(fileHybridName) as f:
		lineNumber = 0
		for line in f:
			lineNumber+=1
			data = line.split()
			if lineNumber % 2 == 1 and lineNumber>1:
				de_50.append(float(data[0]))
				de_100.append(float(data[1]))
				de_200.append(float(data[2]))
				di_50.append(float(data[3]))
				di_100.append(float(data[5]))
				di_200.append(float(data[7]))
				da_50.append(float(data[9]))
				da_100.append(float(data[10]))
				da_200.append(float(data[11]))
				n_orbits_50.append(float(data[12]))
				n_orbits_100.append(float(data[13]))
				n_orbits_200.append(float(data[14]))
	for n in range(len(Qa)):
		print(Qa[n], eStar[n], iStar[n], OmegaStar[n], omegaStar[n], de_3body[n], de_sa[n], de_50[n], de_100[n], de_200[n], di_3body[n], di_sa[n], di_50[n], di_100[n], di_200[n], da_3body[n], da_sa[n], da_50[n], da_100[n], da_200[n], n_orbits_50[n], n_orbits_100[n], n_orbits_200[n], file=fileResult)
	fileResult.close()

if __name__ in ('__main__'):
	# m1 = 5 * u.solMass
	# m2 = 5 * u.solMass
	# a = 1 * u.AU
	# e = 0.1
	# i = 0.
	# Omega = 0
	# omega = 0
	# meanAnomaly = 0

	# m3 = 5 * u.solMass
	# v = 3 * u.km/u.s
	# iStar = 1.17282676772
	# OmegaStar = 1.1054630187304297
	# omegaStar = 6.0707080371921744

	# r3max = 5000
	# r3max_nbody = 100

	# aStar = -G*(m1+m2+m3)/v**2
	# Q = 0.5*a
	# eStar = 1 - Q/aStar

	m1 = 1 * u.solMass
	m2 = 1 * u.solMass
	a = 1 * u.AU
	e = 0.9
	i = 1.
	Omega = 1
	omega = 1
	meanAnomaly = 1

	m3 = 1 * u.solMass
	v = 3 * u.km/u.s
	# p = 1.72123522267e-11*1.98892e+24 * u.m
	# eStar = 1.0114948642451393	
	iStar = 2.
	OmegaStar = 0.1
	omegaStar = np.pi/2

	r3max = 200
	r3max_nbody = 100

	aStar = -G*(m1+m2+m3)/v**2
	Q = 10*a
	eStar = 1 - Q/aStar

	# energy = -G*m1*m2/2/a + m3*v**2/2
	# energy1 = -G*m1*m2/100/a-G*m1*m3/100/a-G*m3*m2/100/a
	# print(energy)
	# print(energy1)

	# compare_new (m1, m2, a, e, m3, v, r3maxNbody = 100, Qamin=0.5, Qamax=3, dQa=0.5, N_3body = 100, r3max=5000, n_orbits_max=20)

	# file = open('/Users/axr6631/Dropbox/code/flybys-master/anomaly.txt', 'w+')
	# for n in range(0, 1000):
	# 	meanAnomaly = n/1000*2*np.pi
	# 	print(meanAnomaly, true_anomaly_from_mean(e, meanAnomaly), file=file)
	# file.close()

	# meanAnomalyDependence (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, N=1000, r3max=r3max, r3max_nbody=r3max_nbody, dt=100)

	result, third_body, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, dt=1000)
	# print(dv_binary)
	# dv_binary, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3_i=-r3max, r3_f=r3max)
	print(result, third_body, finalTime.to(u.s).value, dv_binary.to(u.m/u.s).value, a_new, e_new, i_new, Omega_new, omega_new, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_new, sep='\n')
	# print(dv_binary)
	# for meanAnomaly0 in np.linspace(0, 2*np.pi, 100):
	# result, third_body_final, dv_binary, a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa, n_3body = scattering_hybrid (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max, r3max_nbody=r3max_nbody, meanAnomaly0=-1)
	# print(e_fin_sa-e)
	# print(omega_fin_sa-omega)

	# a_fin_0, e_fin_0, i_fin_0, Omega_fin_0, omega_fin_0 = scattering_SA_bound_test (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3_i=-1e4, r3_f=1e4)
	# print(a_fin_0, e_fin_0, i_fin_0, Omega_fin_0, omega_fin_0)	

	# i = 0
	# Omega = 0
	# omega = 0
	# random.seed()
	# meanAnomaly = random.random()*2*np.pi
	# iStar = np.arccos(random.random()*2-1)
	# OmegaStar = random.random()*2*np.pi
	# omegaStar = random.random()*2*np.pi
	# print(Q/a, iStar, OmegaStar, omegaStar, meanAnomaly)
	# result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body, v_new, p_new, iStar_new, OmegaStar_new, omegaStar_new, r3new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=r3max)
	# print(result)

	# compareHybrid_fixed_r3maxNbody (m1, m2, a, 0.999, m3, v, r3maxNbody=100, Qamin=5, Qamax=20, dQa=1.5, N_3body=100, r3max=5000)
	# compareHybridOnly (m1, m2, a, 0.999, m3, v, r3max=5000, n_orbits_max=20)
	# makeDatabase(0.999)
	# makeDatabase(0.1)

	# QDependenceHybrid (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, iStar, OmegaStar, omegaStar, Qamin=5, Qamax=50, N=100, r3max=r3max, r3max_nbody=r3max_nbody)
	# r3maxDependence (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max1=5000, r3max2=50000, N=5)
	# r3maxDependence_sa (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar, r3max2=50000) 
	
	# aStar = -G*(m1+m2+m3)/v**2
	# eStar = np.sqrt(1+p**2/aStar**2)
	# Q = -aStar*(eStar-1)
	# print(Q/a)

	# Q = 40*a
	# aStar = -G*(m1+m2+m3)/v**2
	# eStar = 1-Q/aStar
	# p = -aStar*np.sqrt(eStar**2-1)

	# print (eStar)
	# meanAnomaly = 0.5

	# meanAnomaly = 0
	# start = time.time()
	# for i in range(300):
	# 	result, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, v_new, p_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, dt=100, r3max=50)
	# end = time.time()
	# # print(finalTime)
	# print(end-start)

	# result, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, v_new, p_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, dt=100, r3max=100)
	# print(finalTime)
	# print(e_new)
	# print(a_new/a-1)
	# meanAnomaly = 3
	# result, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, v_new, p_new, iStar_new, OmegaStar_new, omegaStar_new, r3max_new = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3max=r3max)
	# print(a_new/a-1)
	# dE = G*m1*m2/2*(1/a-1/a_new) + (m1+m2)*dv_binary**2/2
	# EStar = m3*v**2/2
	# print(dv_binary, dE/EStar)
	# print(e_new-e)
	# # print(v_new.to(u.km/u.s))
	# dE = G*m1*m2/2*(1/a-1/a_new) #+ (m1+m2)*dv_binary**2/2
	# EStar = m3*v**2/2
	# print(result)
	# print(dE/EStar)

	# a_fin_sa, e_fin_sa, i_fin_sa, Omega_fin_sa, omega_fin_sa = scattering_SA (m1, m2, a, e, i, Omega, omega, m3, aStar, eStar, iStar, OmegaStar, omegaStar, r3_i=-r3max, r3_f=r3max)
	# print((i_fin_sa-i).to(u.rad))
	
	# print(p/p_new, v/v_new, iStar-iStar_new, OmegaStar-OmegaStar_new, omegaStar-omegaStar_new)

	# homeDir='/Users/axr6631'
	# # homeDir='/home/alexander'
	# fileName1 = homeDir+'/Dropbox/code/flybys-master/compareHybrid_allAngles_r3maxNbody=100_m1=5_m2=5_a=1_e=0.999_m3=5_v=3.txt'
	# fileName2 = homeDir+'/Dropbox/code/flybys-master/de(t)1_Q=30_3body_meanAnomaly=1.5.txt'
	# fileName3 = homeDir+'/Dropbox/code/flybys-master/meanAnomalyDependence1_3body_Q=30_r3max=200_dt=10.txt'
	# fileName4 = homeDir+'/Dropbox/code/flybys-master/meanAnomalyDependence1_3body_Q=30_r3max=200_dt=1.txt'
	# fileName5 = homeDir+'/Dropbox/code/flybys-master/meanAnomalyDependence1_3body_Q=30_r3max=200_dt=0.1.txt'
	# fileName_sa = '/home/alexander/Dropbox/code/flybys-master/r3maxDependence4_sa.txt'
	# makePlot(fileName1, 0.999, homeDir+'/Dropbox/code/flybys-master/outliers_omegas_e=0.999.pdf')
	# analyze_comparison_results (fileName1)
	
	# compare (m1, m2, a, e, m3, v, N_3body=1000)

	# v = 3 * u.km/u.s
	# fileName1 = '/home/alexander/Dropbox/code/flybys-master/destroyed_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'.txt'
	# v = 30 * u.km/u.s
	# fileName2 = '/home/alexander/Dropbox/code/flybys-master/destroyed_m1='+str(m1.to(u.solMass))+'_m2='+str(m2.to(u.solMass))+'_a='+str(a.to(u.AU))+'_e='+str(e)+'_m3='+str(m3.to(u.solMass))+'_v='+str(v.to(u.km/u.s))+'.txt'
	# makePlot(fileName1, fileName2, '/home/alexander/Dropbox/code/flybys-master/destroyed_e=0.1.pdf', e, v.to(u.km/u.s))

	# destroyedFraction (m1, m2, a, e, m3, v, N_3body = 100)

	# Q = 0.2*a
	# aStar = -G*(m1+m2+m3)/v**2
	# eStar = 1 - Q/aStar
	# p = -aStar*np.sqrt(eStar**2-1)
	# result, finalTime, a_fin_3body, e_fin_3body, i_fin_3body, Omega_fin_3body, omega_fin_3body = scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=400)
	# print(result)

	# compare (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, iStar, OmegaStar, omegaStar)

	# print(a.to(u.AU),e,i,Omega,omega)
	# print(scattering_SA (m1, m2, a, e, i, Omega, omega, m3, v, p, iStar, OmegaStar, omegaStar))
	# result, finalTime, a_new, e_new, i_new, Omega_new, omega_new = scattering(m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar, r3max=400)
	# print(a_new, e_new, i_new, Omega_new, omega_new)



	# plotHistogram('histograms.pdf')
	# print ((np.sqrt((G*(m1+m2)/v**2)**2 + p**2) - G*(m1+m2)/v**2) / a)

	# converter = nbody_system.nbody_to_si(m1+m2, a)
	# gravity = Mikkola(converter)
	# help(gravity.parameters)

	
	# plot_distances('distances_2.pdf')

	# rp_max = 0.5*a
	# p_max = rp_max * np.sqrt (1.+2.*G*(m1+m2)/rp_max/v**2)
	# sum = 0*u.m**2/u.s**2

	# p_max = 10*a #corresponds to rp_max~0.5 given v~0.1
	# N = 10000
	# random.seed()
	# for n in range(N):
	# 	meanAnomaly = random.random()*2*np.pi
	# 	iStar = np.arccos(random.random()*2-1)
	# 	OmegaStar = random.random()*2*np.pi
	# 	omegaStar = random.random()*2*np.pi
	# 	p = p_max*np.sqrt(random.random())
	# 	scattering(m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, v, p, iStar, OmegaStar, omegaStar)

	# H1 = 2*np.pi*v**2/(G**2*m1*m2) * p_max**2 * sum/N
	# print(H1)