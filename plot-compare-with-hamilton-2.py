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
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"

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

# example 1 - a=30, e=0.2
# exmaple 2 - a=250, e=0.6
# example 3 - a=49, e=0.5
# example 4 - a=30, e=0.5

n_sa = 2
t_sa_example2 = [[] for i in range(n_sa)]
a_sa_example2 = [[] for i in range(n_sa)]
e_sa_example2 = [[] for i in range(n_sa)]
t_sa_example3 = [[] for i in range(n_sa)]
a_sa_example3 = [[] for i in range(n_sa)]
e_sa_example3 = [[] for i in range(n_sa)]
t_da_example2 = []
a_da_example2 = []
e_da_example2 = []
t_da_example3 = []
a_da_example3 = []
e_da_example3 = []
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
Omega = ['0', 'pi2']
label = ['0', r'\pi/2']
linestyle = ['k:', 'b--']

# My results

for i in range(n_sa):
	lineNumber=0
	with open(folder+f"example2_Omega={Omega[i]}_n=100.txt") as f:
		for line in f:
			lineNumber += 1
			data = line.split()
			if isfloat(data[0]) and float(data[0])>t_min:
				if float(data[0])>t_max:
					break
				t_sa_example2[i].append(float(data[0])/1e9)
				a_sa_example2[i].append(float(data[1]))
				e_sa_example2[i].append(np.log10(1-float(data[2])))

for i in range(n_sa):
	lineNumber=0
	with open(folder+f"example3_Omega={Omega[i]}_n=100.txt") as f:
		for line in f:
			lineNumber += 1
			data = line.split()
			if isfloat(data[0]) and float(data[0])>t_min:
				if float(data[0])>t_max:
					break
				t_sa_example3[i].append(float(data[0])/1e9)
				a_sa_example3[i].append(float(data[1]))
				e_sa_example3[i].append(np.log10(1-float(data[2])))

lineNumber=0
with open(folder+f"example2_Omega=0.txt") as f:
	for line in f:
		lineNumber += 1
		data = line.split()
		if isfloat(data[0]) and float(data[0])>t_min:# and lineNumber%1000==1:
			if float(data[0])>t_max:
				break
			t_n30.append(float(data[0])/1e9)
			a_n30.append(float(data[1]))
			e_n30.append(np.log10(1-float(data[2])))

lineNumber=0
with open(folder+f"example2_Omega=0_n=100.txt") as f:
	for line in f:
		lineNumber += 1
		data = line.split()
		if isfloat(data[0]) and float(data[0])>t_min:# and lineNumber%1000==1:
			if float(data[0])>t_max:
				break
			t_n100.append(float(data[0])/1e9)
			a_n100.append(float(data[1]))
			e_n100.append(np.log10(1-float(data[2])))

lineNumber=0
with open(folder+f"example2_Omega=0_n=300.txt") as f:
	for line in f:
		lineNumber += 1
		data = line.split()
		if isfloat(data[0]) and float(data[0])>t_min:# and lineNumber%1000==1:
			if float(data[0])>t_max:
				break
			t_n300.append(float(data[0])/1e9)
			a_n300.append(float(data[1]))
			e_n300.append(np.log10(1-float(data[2])))

lineNumber=0
with open(folder+f"example2_Omega=0_n=1000.txt") as f:
	for line in f:
		lineNumber += 1
		data = line.split()
		if isfloat(data[0]) and float(data[0])>t_min:# and lineNumber%1000==1:
			if float(data[0])>t_max:
				break
			t_n1000.append(float(data[0])/1e9)
			a_n1000.append(float(data[1]))
			e_n1000.append(np.log10(1-float(data[2])))

# Chris' DA result

dt_values = [29.625792, 5635.849503, 16605.296259, 922.723818]
dt_example2 = dt_values[1]
dt_example3 = dt_values[2]

t = 0
with open(folder+f"da plots/a_Example_2") as f:
	for line in f:
		if t>t_max or not isfloat(line):
			break
		if isfloat(line) and t>t_min:
			t_da_example2.append(t/1e9)
			a_da_example2.append(float(line))
		t += dt_example2
t = 0
with open(folder+f"da plots/e_Example_2") as f:
	for line in f:
		if t>t_max or not isfloat(line):
			break
		if isfloat(line) and t>t_min:
			e_da_example2.append(np.log10(1-float(line)))
		t += dt_example2

t = 0
with open(folder+f"da plots/a_Example_3") as f:
	for line in f:
		if t>t_max or not isfloat(line):
			break
		if isfloat(line) and t>t_min:
			t_da_example3.append(t/1e9)
			a_da_example3.append(float(line))
		t += dt_example3
t = 0
with open(folder+f"da plots/e_Example_3") as f:
	for line in f:
		if t>t_max or not isfloat(line):
			break
		if isfloat(line) and t>t_min:
			e_da_example3.append(np.log10(1-float(line)))
		t += dt_example3

figure = plt.figure(figsize=(12, 12))

plt.subplot(4, 3, 1)
plt.minorticks_on() 
plt.tick_params(labelsize=14)
plt.xlabel(r'$t$ [Gyr]', fontsize=16)
plt.ylabel(r'$a$ [AU]', fontsize=16)
plt.plot(t_da_example2, a_da_example2, 'r')#, label='DA (Hamilton \& Rafikov 2022)')
for i in range(n_sa):
	plt.plot(t_sa_example2[i], a_sa_example2[i], linestyle[i])#, label=rf'SA (this paper), $\Omega_0={label[i]}$')
plt.xlim(0.01, 10)
plt.xscale('log')
ax = plt.gca()
plt.text(0.025, 0.05 , '(A1)', fontsize=16, transform=ax.transAxes)

plt.subplot(4, 3, 2)
plt.minorticks_on() 
plt.tick_params(labelsize=14)
plt.xlabel(r'$t$ [Gyr]', fontsize=16)
plt.ylabel(r'$a$ [AU]', fontsize=16)
plt.plot(t_da_example3, a_da_example3, 'r')#, label='DA (Hamilton \& Rafikov 2022)')
for i in range(n_sa):
	plt.plot(t_sa_example3[i], a_sa_example3[i], linestyle[i])#, label=rf'SA (this paper), $\Omega_0={label[i]}$')
plt.xlim(0.01, 20)
plt.xscale('log')
ax = plt.gca()
plt.text(0.025, 0.05, '(B1)', fontsize=16, transform=ax.transAxes)

plt.subplot(4, 3, 3)
plt.minorticks_on() 
plt.tick_params(labelsize=14)
plt.xlabel(r'$t$ [Gyr]', fontsize=16)
plt.ylabel(r'$a$ [AU]', fontsize=16)
plt.plot(t_n30, a_n30, 'k')#, label=r'$n=30$')
plt.plot(t_n100, a_n100, 'b--')#, label=r'$n=100$')
plt.plot(t_n300, a_n300, 'r:')#, label=r'$n=300$')
plt.plot(t_n1000, a_n1000, 'g-.')#, label=r'$n=1000$')
plt.xlim(0.01, 2)
plt.xscale('log')
ax = plt.gca()
plt.text(0.025, 0.05, '(C1)', fontsize=16, transform=ax.transAxes)

plt.subplot(4, 1, 2)
plt.minorticks_on() 
plt.tick_params(labelsize=14)
plt.xlabel(r'$t$ [Gyr]', fontsize=16)
plt.ylabel(r'$\log_{10}(1-e)$', fontsize=16)
plt.plot(t_da_example2, e_da_example2, 'r')
for i in range(n_sa):
	plt.plot(t_sa_example2[i], e_sa_example2[i], linestyle[i])
ax = plt.gca()
plt.text(0.007, 0.05, '(A2)', fontsize=16, transform=ax.transAxes)

plt.subplot(4, 1, 3)
plt.minorticks_on() 
plt.tick_params(labelsize=14)
plt.xlabel(r'$t$ [Gyr]', fontsize=16)
plt.ylabel(r'$\log_{10}(1-e)$', fontsize=16)
plt.plot(t_da_example3, e_da_example3, 'r')
for i in range(n_sa):
	plt.plot(t_sa_example3[i], e_sa_example3[i], linestyle[i])
ax = plt.gca()
plt.text(0.007, 0.05, '(B2)', fontsize=16, transform=ax.transAxes)

plt.subplot(4, 1, 4)
plt.minorticks_on() 
plt.tick_params(labelsize=14)
plt.xlabel(r'$t$ [Gyr]', fontsize=16)
plt.ylabel(r'$\log_{10}(1-e)$', fontsize=16)
plt.plot(t_n30, e_n30, 'k')
plt.plot(t_n100, e_n100, 'b--')
plt.plot(t_n300, e_n300, 'r:')
plt.plot(t_n1000, e_n1000, 'g-.')
ax = plt.gca()
plt.text(0.007, 0.05, '(C2)', fontsize=16, transform=ax.transAxes)

# plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig("output/for the paper/comparison-with-da.pdf")