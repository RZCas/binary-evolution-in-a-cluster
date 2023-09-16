import numpy as np
import glob, os, subprocess
import matplotlib
from pathlib import Path
from matplotlib import pyplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{color}"

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

root_folder = 'storage/m1=m2=10/'
simulation_sets = ['mtotal=1e5/','mtotal=1e6/','mtotal=1e5,enc_only/','mtotal=1e6,enc_only/']
for simulation_set in simulation_sets:
	input_folder = root_folder + simulation_set
	output_folder = root_folder + simulation_set[:-1] + '-e-distribution/'
	Path(output_folder).mkdir(parents=True, exist_ok=True)
	t_0 = 1 #Gyr
	dt = 0.1 #Gyr
	n = 3000 #Number of points uniformly distributed between t_0 and t_0+dt where we calculate the eccentricity
	t = np.linspace(t_0, t_0+dt, n)

	index = 0
	filepath = input_folder + str(index) + '.txt'
	while os.path.isfile(filepath):
		lastline = subprocess.check_output(['tail', '-1', filepath])[:-1]
		if isfloat(lastline.split()[0]) and float(lastline.split()[0])/1e9 < t_0+dt:
			index += 1
			filepath = input_folder + str(index) + '.txt'
			continue
		e = [-1]*n
		i = 0
		t_previous = 0
		e_previous = 0.5
		with open(filepath) as f:
			for line in f:
				data = line.split()
				if len(data)>=11 and isfloat(data[0]) and isfloat(data[1]) and isfloat(data[10]):
					t_current = float(data[0])/1e9
					e_current = float(data[10])
					if t_current < t[i]:
						t_previous = t_current
						continue
					else:
						while t[i] <= t_current:
							e[i] = e_previous + (e_current - e_previous) * (t[i] - t_previous) / (t_current - t_previous)
							i += 1
							if i == n:
								break
						if i == n:
							break

		if i==n:
			output_file = open(output_folder + str(index) + '-e-distribution.txt', 'w+')
			for ecc in e:
				print(ecc, file=output_file)
			fig = pyplot.figure(figsize=(6, 4)) 
			pyplot.hist (e, histtype='step', density=True, bins=np.linspace(0, 1, 31))
			pyplot.xlabel (r'$e$', fontsize=16)
			pyplot.tight_layout()
			pyplot.savefig(output_folder + str(index) + '-e-distribution.pdf')
			pyplot.close()
		index += 1
		filepath = input_folder + str(index) + '.txt'



