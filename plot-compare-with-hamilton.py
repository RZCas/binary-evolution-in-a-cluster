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
figure = pyplot.figure(figsize=(12, 12))

folder = "output/chris_thesis_corrected/"

# t_min = 0.3e7
# t_max = 0.37e7
# t_min = 0.0083e10
# t_max = 0.00843e10
# t_min = 0.0245e10
# t_max = 0.026e10
t_min = 0.3e7
t_max = 0.35e7

example = 4
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
if example == 1 or example == 4:
	Omega = ['0', 'pi2']
	label = ['0', r'\pi/2']
else:
	Omega = ['0', 'pi4']
	label = ['0', r'\pi/4']
linestyle = ['k', 'b']
# linestyle = ['k', 'k--', 'b', 'b--', 'g', 'g--']
t_da = []
a_da = []
e_da = []
t_old = []
a_old = []
e_old = []

# My results

for i in range(n_sa):
	lineNumber=0
	with open(folder+f"example{example}_Omega={Omega[i]}.txt") as f:
		for line in f:
			lineNumber += 1
			data = line.split()
			if isfloat(data[0]) and float(data[0])>t_min:
				if float(data[0])>t_max:
					break
				t_sa[i].append(float(data[0]))
				a_sa[i].append(float(data[1]))
				e_sa[i].append(np.log10(1-float(data[2])))

# lineNumber=0
# with open(folder+"a_in=250_e_in=0.6_evolution.txt") as f:
# 	for line in f:
# 		lineNumber += 1
# 		data = line.split()
# 		if isfloat(data[0]) and float(data[0])>t_min:
# 			if float(data[0])>t_max:
# 				break
# 			t_old.append(float(data[0]))
# 			a_old.append(float(data[1]))
# 			e_old.append(np.log10(1-float(data[2])))

# Chris' DA results

dt_values = [29.625792, 5635.849503, 16605.296259, 922.723818]
dt = dt_values[example-1]

t = 0
with open(folder+f"da plots/a_Example_{example}") as f:
	for line in f:
		if t>t_max or not isfloat(line):
			break
		if isfloat(line) and t>t_min:
			t_da.append(t)
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

n_plots = 2

plot_a = figure.add_subplot(n_plots,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$a$ [AU]', fontsize=16)
for i in range(n_sa):
	plot_a.plot(t_sa[i], a_sa[i], linestyle[i], label=rf'SA (this paper), $\Omega_0={label[i]}$')
plot_a.plot(t_da, a_da, 'r', label='DA (Hamilton \& Rafikov 2022)')
# plot_a.plot(t_old, a_old, 'r', label=r'SA (old), $\Omega_0=0$')
plot_a.legend()

plot_e = figure.add_subplot(n_plots,1,2)
ax = pyplot.gca()
ax.minorticks_on() 
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$\log_{10}(1-e)$', fontsize=16)
for i in range(n_sa):
	plot_e.plot(t_sa[i], e_sa[i], linestyle[i])
plot_e.plot(t_da, e_da, 'r')
# plot_e.plot(t_old, e_old, 'r')

pyplot.tight_layout()
pyplot.savefig(folder+f"comparison-with-da-example-{example}-emin1.pdf")