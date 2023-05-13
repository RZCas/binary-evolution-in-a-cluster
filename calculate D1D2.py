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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{color}"

# directories = ['storage', 'storage/m1=m2=10']
directories = ['storage/m1=m2=10']
for directory in directories:
	types = os.listdir(directory)
	for type in types:
		if os.path.isfile(directory+'/'+type):	#if it's a file rather than a folder, ignore it
			continue
		root_dir = directory+'/'+type+'/'
	outfile = open (root_dir + type + '-D1D2.txt', 'w+')
	print('index D_1,tidal/(D_1,tidal+D2_flybys)', file=outfile)
	for index in range(0,1000):
		filepath = root_dir + str(index) + '.txt'
		if not os.path.isfile(filepath):
			continue
		print(type, index)

		lineNumber = 0
		shift = 0
		de_abs_tidal = 0
		de2_flybys = 0
		e = 0

		with open(filepath) as f:
			for line in f:
				lineNumber+=1
				data = line.split()
				if len(data) > 1:
					if data[0]=="perturber:" and lineNumber%3==0:
						shift=1
					if isfloat(data[0]) and isfloat(data[1]):
						if len(data)<11:
							break
						t = float(data[0])
						e_prev = e
						e = float(data[10])
						if lineNumber%3==1+shift and lineNumber>1:	#before encounter
							if len(data)<18:
								break
							de_abs_tidal += abs(e - e_prev)
						if lineNumber%3==0+shift and lineNumber>3:	#after encounter
							de2_flybys += (e - e_prev)**2

		D1_tidal = de_abs_tidal / t
		D2 = de2_flybys / t
		print(index, D1_tidal / (D1_tidal + D2), file=outfile)