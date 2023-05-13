import numpy as np

def make_plot_omega(fileName1, e, save_file):
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
	plot.plot(Qa, outlierNumber_sa_Omega/totalNumber, 'r', label=r'$\Omega$, orbit-averaged')
	plot.plot(Qa, outlierNumber_hybrid_omega/totalNumber, 'k--', label=r'$\omega$, hybrid')
	plot.plot(Qa, outlierNumber_sa_omega/totalNumber, 'r--', label=r'$\omega$, orbit-averaged')
	# plot.plot(x2, y2, 'r', label='SA')
	# plot.plot(x2neg, y2neg, 'r--')
	# plot.plot(x2, y2, 'k--', label=r'3-body, $\delta t=100$')
	# plot.plot(x_sa, y_sa, 'r', label=r'SA')
	plot.legend(fontsize=16, frameon=False)
	# plot.set_xlim(500, 5000)
	# plot.set_ylim(1e-8, 1e-5)
	pyplot.tight_layout()
	pyplot.savefig(save_file)

# compare_m1=5_m2=5_a=1_e=0.999_m3=5_v=3.txt
def make_plot_i(fileName1, e, save_file):
	Qamin = 5
	Qamax = 20
	N = 10
	Qa = [Qamin+n/N*(Qamax-Qamin) for n in range(N)]
	totalNumber = np.array([0]*N)
	outlierNumber50 = np.array([0]*N)
	outlierNumber100 = np.array([0]*N)
	outlierNumber200 = np.array([0]*N)
	outlierNumber_sa = np.array([0]*N)

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
				di_3body = float(data[10])
				di_sa = float(data[11])
				di_50 = float(data[12])
				di_100 = float(data[13])
				di_200 = float(data[14])
				n = round(N*(QaValue-Qamin)/(Qamax-Qamin))
				totalNumber[n] += 1
				# if float(data[-1])>0: print('yes!')
				if np.abs(di_50/di_3body-1) > 0.2: outlierNumber50[n] += 1
				if np.abs(di_100/di_3body-1) > 0.2: outlierNumber100[n] += 1
				if np.abs(di_200/di_3body-1) > 0.2: 
					outlierNumber200[n] += 1
					# print(line, file=outlierFile)	
				if np.abs(di_sa/di_3body-1) > 0.2: 
					outlierNumber_sa[n] += 1
					
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
	ax.set_ylabel(r'fraction of $\Delta i$ discrepancies $>20\%$', fontsize=16)
	# pyplot.xscale('log')
	# pyplot.yscale('log')
	pyplot.text(0.3, 0.9, '$e='+str(e)+'$', fontsize=16, transform=ax.transAxes)
	# plot.plot(x40, y40, 'k', label=r'3-body, $r_{3,\mathrm{max}}=40$')
	plot.plot(Qa, outlierNumber50/totalNumber, 'k', label=r'$r_{\rm 3body}=50$')
	plot.plot(Qa, outlierNumber100/totalNumber, 'k--', label=r'$r_{\rm 3body}=100$')
	plot.plot(Qa, outlierNumber200/totalNumber, 'k:', label=r'$r_{\rm 3body}=200$')
	plot.plot(Qa, outlierNumber_sa/totalNumber, 'r', label=r'orbit-averaged')
	# plot.plot(x2, y2, 'r', label='SA')
	# plot.plot(x2neg, y2neg, 'r--')
	# plot.plot(x2, y2, 'k--', label=r'3-body, $\delta t=100$')
	# plot.plot(x_sa, y_sa, 'r', label=r'SA')
	plot.legend(fontsize=16, frameon=False)
	# plot.set_xlim(500, 5000)
	# plot.set_ylim(1e-8, 1e-5)
	pyplot.tight_layout()
	pyplot.savefig(save_file)

def make_plot_e(fileName1, e, save_file):
	Qamin = 5
	Qamax = 20
	N = 10
	Qa = [Qamin+n/N*(Qamax-Qamin) for n in range(N)]
	totalNumber = np.array([0]*N)
	outlierNumber50 = np.array([0]*N)
	outlierNumber100 = np.array([0]*N)
	outlierNumber200 = np.array([0]*N)
	outlierNumber_sa = np.array([0]*N)

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
				di_3body = float(data[10])
				di_sa = float(data[11])
				di_50 = float(data[12])
				di_100 = float(data[13])
				di_200 = float(data[14])
				n = round(N*(QaValue-Qamin)/(Qamax-Qamin))
				totalNumber[n] += 1
				# if float(data[-1])>0: print('yes!')
				if np.abs(de_50/de_3body-1) > 0.2: outlierNumber50[n] += 1
				if np.abs(de_100/de_3body-1) > 0.2: outlierNumber100[n] += 1
				if np.abs(de_200/de_3body-1) > 0.2: 
					outlierNumber200[n] += 1
					# print(line, file=outlierFile)	
				if np.abs(de_sa/de_3body-1) > 0.2: 
					outlierNumber_sa[n] += 1
					
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
	ax.set_ylabel(r'fraction of $\Delta e$ discrepancies $>20\%$', fontsize=16)
	# pyplot.xscale('log')
	# pyplot.yscale('log')
	pyplot.text(0.3, 0.9, '$e='+str(e)+'$', fontsize=16, transform=ax.transAxes)
	# plot.plot(x40, y40, 'k', label=r'3-body, $r_{3,\mathrm{max}}=40$')
	plot.plot(Qa, outlierNumber50/totalNumber, 'k', label=r'$r_{\rm 3body}=50$')
	plot.plot(Qa, outlierNumber100/totalNumber, 'k--', label=r'$r_{\rm 3body}=100$')
	plot.plot(Qa, outlierNumber200/totalNumber, 'k:', label=r'$r_{\rm 3body}=200$')
	plot.plot(Qa, outlierNumber_sa/totalNumber, 'r', label=r'orbit-averaged')
	# plot.plot(x2, y2, 'r', label='SA')
	# plot.plot(x2neg, y2neg, 'r--')
	# plot.plot(x2, y2, 'k--', label=r'3-body, $\delta t=100$')
	# plot.plot(x_sa, y_sa, 'r', label=r'SA')
	plot.legend(fontsize=16, frameon=False)
	# plot.set_xlim(500, 5000)
	# plot.set_ylim(1e-8, 1e-5)
	pyplot.tight_layout()
	pyplot.savefig(save_file)

make_plot_omega('compareHybrid_allAngles_r3maxNbody=100_m1=5_m2=5_a=1_e=0.999_m3=5_v=3.txt', 0.999, 'outliers_omegas_e=0.999.pdf')
make_plot_omega('compareHybrid_allAngles_r3maxNbody=100_m1=5_m2=5_a=1_e=0.1_m3=5_v=3.txt', 0.1, 'outliers_omegas_e=0.1.pdf')
make_plot_i('compare_m1=5_m2=5_a=1_e=0.999_m3=5_v=3.txt', 0.999, 'outliers_di_e=0.999.pdf')
make_plot_i('compare_m1=5_m2=5_a=1_e=0.1_m3=5_v=3.txt', 0.1, 'outliers_di_e=0.1.pdf')
make_plot_e('compare_m1=5_m2=5_a=1_e=0.999_m3=5_v=3.txt', 0.999, 'outliers_e=0.999.pdf')
make_plot_e('compare_m1=5_m2=5_a=1_e=0.1_m3=5_v=3.txt', 0.1, 'outliers_e=0.1.pdf')
