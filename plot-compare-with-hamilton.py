import os
import glob
import numpy as np

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
figure = pyplot.figure(figsize=(12, 9))

folder = "output/chris_thesis_corrected/"

# t_min = 0.3e7
# t_max = 0.37e7
# t_min = 0.0083e10
# t_max = 0.00843e10
# t_min = 0.0245e10
# t_max = 0.026e10
# t_min = 0.3e7
# t_max = 0.35e7
# t_min = 1.68e7
# t_max = 1.72e7
t_min = 0
t_max = 1e100

example = 2
# example 1 - a=30, e=0.2
# exmaple 2 - a=250, e=0.6
# example 3 - a=49, e=0.5
# example 4 - a=30, e=0.5

n_sa = 2
t_sa = [[] for i in range(n_sa)]
a_sa = [[] for i in range(n_sa)]
e_sa = [[] for i in range(n_sa)]
# Omega = ['0', 'pi4', 'pi2', '3pi4', 'pi', '3pi2']
# label = ['0', r'\pi/4', r'\pi/2', r'3\pi/4', r'\pi', r'3\pi/2']
# if example == 1 or example == 3 or example == 4:
Omega = ['0', 'pi2']
label = ['0', r'\pi/2']
# else:
# 	Omega = ['0', 'pi4']
# 	label = ['0', r'\pi/4']
linestyle = ['k', 'b']
# linestyle = ['k', 'k--', 'b', 'b--', 'g', 'g--']
t_da = []
a_da = []
e_da = []
t_old = []
a_old = []
e_old = []
t_n30 = []
a_n30 = []
e_n30 = []
t_n100 = []
a_n100 = []
e_n100 = []
t_n300 = []
a_n300 = []
e_n300 = []
t_n1000 = []
a_n1000 = []
e_n1000 = []

# My results

# for i in range(n_sa):
# 	lineNumber=0
# 	with open(folder+f"example{example}_Omega={Omega[i]}_n=100.txt") as f:
# 		for line in f:
# 			lineNumber += 1
# 			data = line.split()
# 			if isfloat(data[0]) and float(data[0])>t_min:
# 				if float(data[0])>t_max:
# 					break
# 				t_sa[i].append(float(data[0])/1e9)
# 				a_sa[i].append(float(data[1]))
# 				e_sa[i].append(np.log10(1-float(data[2])))

lineNumber=0
with open(folder+f"example{example}_Omega=0.txt") as f:
	for line in f:
		lineNumber += 1
		data = line.split()
		if isfloat(data[0]) and float(data[0])>t_min:
			if float(data[0])>t_max:
				break
			t_n30.append(float(data[0])/1e9)
			a_n30.append(float(data[1]))
			e_n30.append(np.log10(1-float(data[2])))

lineNumber=0
with open(folder+f"example{example}_Omega=0_n=100.txt") as f:
	for line in f:
		lineNumber += 1
		data = line.split()
		if isfloat(data[0]) and float(data[0])>t_min:
			if float(data[0])>t_max:
				break
			t_n100.append(float(data[0])/1e9)
			a_n100.append(float(data[1]))
			e_n100.append(np.log10(1-float(data[2])))

lineNumber=0
with open(folder+f"example{example}_Omega=0_n=300.txt") as f:
	for line in f:
		lineNumber += 1
		data = line.split()
		if isfloat(data[0]) and float(data[0])>t_min:
			if float(data[0])>t_max:
				break
			t_n300.append(float(data[0])/1e9)
			a_n300.append(float(data[1]))
			e_n300.append(np.log10(1-float(data[2])))

lineNumber=0
with open(folder+f"example{example}_Omega=0_n=1000.txt") as f:
	for line in f:
		lineNumber += 1
		data = line.split()
		if isfloat(data[0]) and float(data[0])>t_min:
			if float(data[0])>t_max:
				break
			t_n1000.append(float(data[0])/1e9)
			a_n1000.append(float(data[1]))
			e_n1000.append(np.log10(1-float(data[2])))

# oldFile = ''
# if example == 2:
# 	oldFile = "a_in=250_e_in=0.6_evolution.txt"
# elif example == 3:
# 	oldFile = "a_in=49_e_in=0.5_evolution.txt"
# if oldFile!='':
# 	lineNumber=0
# 	with open(folder+oldFile) as f:
# 		for line in f:
# 			lineNumber += 1
# 			data = line.split()
# 			if isfloat(data[0]) and float(data[0])>t_min:
# 				if float(data[0])>t_max:
# 					break
# 				t_old.append(float(data[0])/1e9)
# 				a_old.append(float(data[1]))
# 				e_old.append(np.log10(1-float(data[2])))

# Chris' DA result

dt_values = [29.625792, 5635.849503, 16605.296259, 922.723818]
dt = dt_values[example-1]

t = 0
with open(folder+f"da plots/a_Example_{example}") as f:
	for line in f:
		if t>t_max or not isfloat(line):
			break
		if isfloat(line) and t>t_min:
			t_da.append(t/1e9)
			a_da.append(float(line))
		t += dt

t = 0
with open(folder+f"da plots/e_Example_{example}") as f:
	for line in f:
		if t>t_max or not isfloat(line):
			break
		if isfloat(line) and t>t_min:
			e_da.append(np.log10(1-float(line)))
		t += dt

gs = figure.add_gridspec(2, 1)
ax1, ax2 = gs.subplots(sharex=True)

ax1.minorticks_on() 
ax1.tick_params(labelsize=14)
ax1.set_ylabel(r'$a$ [AU]', fontsize=16)
# ax1.plot(t_da, a_da, 'r', label='DA (Hamilton \& Rafikov 2022)')
# for i in range(n_sa):
# 	ax1.plot(t_sa[i], a_sa[i], linestyle[i], label=rf'SA (this paper), $\Omega_0={label[i]}$')
ax1.plot(t_n30, a_n30, 'k', label=r'$n=30$')
ax1.plot(t_n100, a_n100, 'b--', label=r'$n=100$')
ax1.plot(t_n300, a_n300, 'r:', label=r'$n=300$')
ax1.plot(t_n1000, a_n1000, 'g-.', label=r'$n=1000$')
# ax1.plot(t_old, a_old, 'r:', label=r'SA (old), $\Omega_0=0$')
ax1.legend(fontsize=16)

ax2.minorticks_on() 
ax2.tick_params(labelsize=14)
ax2.set_xlabel(r'$t$ [Gyr]', fontsize=16)
ax2.set_ylabel(r'$\log_{10}(1-e)$', fontsize=16)
# ax2.plot(t_da, e_da, 'r')
# for i in range(n_sa):
# 	ax2.plot(t_sa[i], e_sa[i], linestyle[i])
ax2.plot(t_n30, e_n30, 'k')
ax2.plot(t_n100, e_n100, 'b--')
ax2.plot(t_n300, e_n300, 'r:')
ax2.plot(t_n1000, e_n1000, 'g--')
# ax2.plot(t_old, e_old, 'r:')


pyplot.tight_layout()
# pyplot.savefig(folder+f"comparison-with-da-example-{example}-publication.pdf")
pyplot.savefig(folder+f"nDependence.pdf")