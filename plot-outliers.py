import numpy as np

def plot_data_omega(fileName1):
	Qamin = 5
	Qamax = 20
	N = 10
	Qa = [Qamin+n/N*(Qamax-Qamin) for n in range(N+1)]
	totalNumber = np.array([0]*(N+1))
	outlierNumber_hybrid_Omega = np.array([0]*(N+1))
	outlierNumber_hybrid_omega = np.array([0]*(N+1))
	outlierNumber_sa_Omega = np.array([0]*(N+1))
	outlierNumber_sa_omega = np.array([0]*(N+1))

	with open(fileName1) as f:
		lineNumber = 0
		for line in f:
			lineNumber+=1
			if lineNumber>1:
				data = line.split()
				QaValue = float(data[0])
				dOmega_3body = float(data[14])
				dOmega_hybrid = float(data[15])
				dOmega_sa = float(data[16])
				domega_3body = float(data[17])
				domega_hybrid = float(data[18])
				domega_sa = float(data[19])
				n = round(N*(QaValue-Qamin)/(Qamax-Qamin))
				totalNumber[n] += 1
				if np.abs(dOmega_hybrid/dOmega_3body-1) > 0.2: outlierNumber_hybrid_Omega[n] += 1
				if np.abs(domega_hybrid/domega_3body-1) > 0.2: outlierNumber_hybrid_omega[n] += 1
				if np.abs(dOmega_sa/dOmega_3body-1) > 0.2: outlierNumber_sa_Omega[n] += 1
				if np.abs(domega_sa/domega_3body-1) > 0.2: outlierNumber_sa_omega[n] += 1
	return Qa, outlierNumber_hybrid_Omega/totalNumber, outlierNumber_sa_Omega/totalNumber, outlierNumber_hybrid_omega/totalNumber, outlierNumber_sa_omega/totalNumber

def plot_data_i(fileName1):
	Qamin = 5
	Qamax = 20
	N = 10
	Qa = [Qamin+n/N*(Qamax-Qamin) for n in range(N+1)]
	totalNumber = np.array([0]*(N+1))
	outlierNumber50 = np.array([0]*(N+1))
	outlierNumber100 = np.array([0]*(N+1))
	outlierNumber200 = np.array([0]*(N+1))
	outlierNumber_sa = np.array([0]*(N+1))

	with open(fileName1) as f:
		lineNumber = 0
		for line in f:
			lineNumber+=1
			if lineNumber>1:
				data = line.split()
				QaValue = float(data[0])
				di_3body = float(data[10])
				di_sa = float(data[11])
				di_50 = float(data[12])
				di_100 = float(data[13])
				di_200 = float(data[14])
				n = round(N*(QaValue-Qamin)/(Qamax-Qamin))
				totalNumber[n] += 1
				if np.abs(di_50/di_3body-1) > 0.2: outlierNumber50[n] += 1
				if np.abs(di_100/di_3body-1) > 0.2: outlierNumber100[n] += 1
				if np.abs(di_200/di_3body-1) > 0.2: outlierNumber200[n] += 1
				if np.abs(di_sa/di_3body-1) > 0.2: outlierNumber_sa[n] += 1
	
	return Qa, outlierNumber50/totalNumber, outlierNumber100/totalNumber, outlierNumber200/totalNumber, outlierNumber_sa/totalNumber		

def plot_data_e(fileName1):
	Qamin = 5
	Qamax = 20
	N = 10
	Qa = [Qamin+n/N*(Qamax-Qamin) for n in range(N+1)]
	totalNumber = np.array([0]*(N+1))
	outlierNumber50 = np.array([0]*(N+1))
	outlierNumber100 = np.array([0]*(N+1))
	outlierNumber200 = np.array([0]*(N+1))
	outlierNumber_sa = np.array([0]*(N+1))

	with open(fileName1) as f:
		lineNumber = 0
		for line in f:
			lineNumber+=1
			if lineNumber>1:
				data = line.split()
				QaValue = float(data[0])
				de_3body = float(data[5])
				de_sa = float(data[6])
				de_50 = float(data[7])
				de_100 = float(data[8])
				de_200 = float(data[9])
				n = round(N*(QaValue-Qamin)/(Qamax-Qamin))
				totalNumber[n] += 1
				if np.abs(de_50/de_3body-1) > 0.2: outlierNumber50[n] += 1
				if np.abs(de_100/de_3body-1) > 0.2: outlierNumber100[n] += 1
				if np.abs(de_200/de_3body-1) > 0.2: outlierNumber200[n] += 1
				if np.abs(de_sa/de_3body-1) > 0.2: outlierNumber_sa[n] += 1
	return Qa, outlierNumber50/totalNumber, outlierNumber100/totalNumber, outlierNumber200/totalNumber, outlierNumber_sa/totalNumber		

Qa, Omega_hybrid_0999, Omega_sa_0999, omega_hybrid_0999, omega_sa_0999 = plot_data_omega('../flybys-master/compareHybrid_allAngles_r3maxNbody=100_m1=5_m2=5_a=1_e=0.999_m3=5_v=3.txt')
Qa, Omega_hybrid_01, Omega_sa_01, omega_hybrid_01, omega_sa_01 = plot_data_omega('../flybys-master/compareHybrid_allAngles_r3maxNbody=100_m1=5_m2=5_a=1_e=0.1_m3=5_v=3.txt')
Qa, i50_0999, i100_0999, i200_0999, i_sa_0999 = plot_data_i('../flybys-master/compare_m1=5_m2=5_a=1_e=0.999_m3=5_v=3.txt')
Qa, i50_01, i100_01, i200_01, i_sa_01 = plot_data_i('../flybys-master/compare_m1=5_m2=5_a=1_e=0.1_m3=5_v=3.txt')
Qa, e50_0999, e100_0999, e200_0999, e_sa_0999 = plot_data_e('../flybys-master/compare_m1=5_m2=5_a=1_e=0.999_m3=5_v=3.txt')
Qa, e50_01, e100_01, e200_01, e_sa_01 = plot_data_e('../flybys-master/compare_m1=5_m2=5_a=1_e=0.1_m3=5_v=3.txt')

import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"

figure = pyplot.figure(figsize=(12, 8)) 
gs = figure.add_gridspec(3, 2, hspace=0)
(e01_plot, e0999_plot), (i01_plot, i0999_plot), (omega01_plot, omega0999_plot) = gs.subplots(sharex=True)

x_label = 0.03
y_label = 0.05
label1 = r'$f(\Delta e>20\%)$'
label2 = r'$f(\Delta i>20\%)$'
label3 = r'$f(\Delta\omega>20\%\;\;\mathrm{or}\;\;\Delta\Omega>20\%)$'

e01_plot.minorticks_on() 
e01_plot.tick_params(labelsize=14)
e01_plot.set_ylabel(label1, fontsize=16)
e01_plot.text(x_label, y_label, '(a)', fontsize=16, transform=e01_plot.transAxes)
e01_plot.text(0.3, 0.9, '$e=0.1$', fontsize=16, transform=e01_plot.transAxes)
e01_plot.plot(Qa, e50_01, 'k', label=r'$r_{\rm 3body}=50$')
e01_plot.plot(Qa, e100_01, 'k--', label=r'$r_{\rm 3body}=100$')
e01_plot.plot(Qa, e200_01, 'k:', label=r'$r_{\rm 3body}=200$')
e01_plot.plot(Qa, e_sa_01, 'r', label=r'orbit-averaged')
e01_plot.legend(fontsize=16, frameon=False)

e0999_plot.minorticks_on() 
e0999_plot.tick_params(labelsize=14)
e0999_plot.set_ylabel(label1, fontsize=16)
e0999_plot.text(x_label, y_label, '(b)', fontsize=16, transform=e0999_plot.transAxes)
e0999_plot.text(0.3, 0.9, '$e=0.999$', fontsize=16, transform=e0999_plot.transAxes)
e0999_plot.plot(Qa, e50_0999, 'k', label=r'$r_{\rm 3body}=50$')
e0999_plot.plot(Qa, e100_0999, 'k--', label=r'$r_{\rm 3body}=100$')
e0999_plot.plot(Qa, e200_0999, 'k:', label=r'$r_{\rm 3body}=200$')
e0999_plot.plot(Qa, e_sa_0999, 'r', label=r'orbit-averaged')
e0999_plot.legend(fontsize=16, frameon=False)

i01_plot.minorticks_on() 
i01_plot.tick_params(labelsize=14)
i01_plot.set_ylabel(label2, fontsize=16)
i01_plot.text(x_label, y_label, '(c)', fontsize=16, transform=i01_plot.transAxes)
i01_plot.text(0.3, 0.9, '$e=0.1$', fontsize=16, transform=i01_plot.transAxes)
i01_plot.plot(Qa, i50_01, 'k', label=r'$r_{\rm 3body}=50$')
i01_plot.plot(Qa, i100_01, 'k--', label=r'$r_{\rm 3body}=100$')
i01_plot.plot(Qa, i200_01, 'k:', label=r'$r_{\rm 3body}=200$')
i01_plot.plot(Qa, i_sa_01, 'r', label=r'orbit-averaged')
i01_plot.legend(fontsize=16, frameon=False)

i0999_plot.minorticks_on() 
i0999_plot.tick_params(labelsize=14)
i0999_plot.set_ylabel(label2, fontsize=16)
i0999_plot.text(x_label, y_label, '(d)', fontsize=16, transform=i0999_plot.transAxes)
i0999_plot.text(0.3, 0.9, '$e=0.999$', fontsize=16, transform=i0999_plot.transAxes)
i0999_plot.plot(Qa, i50_0999, 'k', label=r'$r_{\rm 3body}=50$')
i0999_plot.plot(Qa, i100_0999, 'k--', label=r'$r_{\rm 3body}=100$')
i0999_plot.plot(Qa, i200_0999, 'k:', label=r'$r_{\rm 3body}=200$')
i0999_plot.plot(Qa, i_sa_0999, 'r', label=r'orbit-averaged')
i0999_plot.legend(fontsize=16, frameon=False)

omega01_plot.minorticks_on() 
omega01_plot.tick_params(labelsize=14)
omega01_plot.set_ylabel(label3, fontsize=16)
omega01_plot.text(x_label, y_label, '(e)', fontsize=16, transform=omega01_plot.transAxes)
omega01_plot.text(0.3, 0.9, '$e=0.1$', fontsize=16, transform=omega01_plot.transAxes)
omega01_plot.plot(Qa, Omega_hybrid_01, 'k', label=r'$\Omega$, hybrid')
omega01_plot.plot(Qa, Omega_sa_01, 'r', label=r'$\Omega$, orbit-averaged')
omega01_plot.plot(Qa, omega_hybrid_01, 'k--', label=r'$\omega$, hybrid')
omega01_plot.plot(Qa, omega_sa_01, 'r--', label=r'$\omega$, orbit-averaged')
omega01_plot.legend(fontsize=16, frameon=False)

omega0999_plot.minorticks_on() 
omega0999_plot.tick_params(labelsize=14)
omega0999_plot.set_ylabel(label3, fontsize=16)
omega0999_plot.text(x_label, y_label, '(f)', fontsize=16, transform=omega0999_plot.transAxes)
omega0999_plot.text(0.3, 0.9, '$e=0.999$', fontsize=16, transform=omega0999_plot.transAxes)
omega0999_plot.plot(Qa, Omega_hybrid_0999, 'k', label=r'$\Omega$, hybrid')
omega0999_plot.plot(Qa, Omega_sa_0999, 'r', label=r'$\Omega$, orbit-averaged')
omega0999_plot.plot(Qa, omega_hybrid_0999, 'k--', label=r'$\omega$, hybrid')
omega0999_plot.plot(Qa, omega_sa_0999, 'r--', label=r'$\omega$, orbit-averaged')
omega0999_plot.legend(fontsize=16, frameon=False)

pyplot.tight_layout(rect=[0, 0.03, 1, 0.97])
pyplot.savefig("output/for the paper/outliers.pdf")
pyplot.clf()