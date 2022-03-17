import os
import glob
import numpy as np
gamma = 0.42

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

t_max = 1e60
t_min = 0*1.6e8

folder = "output/noenc_test/cluster/"

n=0

while(n==0):

	t = []
	H = []
	theta = []
	t_0=-1
	lineNumber=0
	epsilon_gr=-1

	with open(folder + str(n) + "test.txt") as f:
		for line in f:
			data = line.split()
			if data[0]=='epsilon_gr': epsilon_gr = float(data[2])

	if epsilon_gr>=0:
		with open(folder + "evolution" + str(n) + ".txt") as f:
			for line in f:
				lineNumber += 1
				data = line.split()
				if isfloat(data[0]) and isfloat(data[1]) and float(data[0])<t_max and float(data[0])>t_min:
					t.append(float(data[0]))
					a = float(data[1])
					e = float(data[2])
					omega = float(data[3])
					i = float(data[4])

					theta.append((1-e**2)*np.cos(i)**2)
					# print(theta[-1])
					H.append((2+3*e**2)*(1-3*gamma*np.cos(i)**2) - 15*gamma*e**2*np.sin(i)**2*np.cos(2*omega) - epsilon_gr/np.sqrt(1-e**2))

		figure = pyplot.figure(figsize=(12, 16))

		plot_H = figure.add_subplot(2,1,1)
		ax = pyplot.gca()
		ax.minorticks_on() 
		ax.tick_params(labelsize=14)
		ax.set_xlabel(r'$t$ [yr]', fontsize=16)
		ax.set_ylabel(r'$H$', fontsize=16)
		plot_H.plot(t, H, 'k')

		plot_theta = figure.add_subplot(2,1,2)
		ax = pyplot.gca()
		ax.minorticks_on() 
		ax.tick_params(labelsize=14)
		ax.set_xlabel(r'$t$ [yr]', fontsize=16)
		ax.set_ylabel(r'$\Theta$', fontsize=16)
		plot_theta.plot(t, theta, 'k')

		pyplot.tight_layout()
		pyplot.savefig(folder + 'conservation' + str(n) + ".pdf")

	n+=1