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

directories = ['storage', 'storage/m1=m2=10']
for directory in directories:
	types = os.listdir(directory)
	for type in types:
		if os.path.isfile(directory+'/'+type):	#if it's a file rather than a folder, ignore it
			continue
		root_dir = directory+'/'+type+'/'
		outfile = open (root_dir + type + '-chi.txt', 'w+')
		print('index chi', file=outfile)
		# print('chi = ', file=outfile)
		for index in range(0,1000):
			filepath = root_dir + str(index) + '.txt'
			if not os.path.isfile(filepath):
				continue
			print(type, index)

			lineNumber = 0
			shift = 0
			de_tidal = 0
			de_flybys = 0
			de_tidal_max = 0
			de_flybys_max = 0
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
							e_prev = e
							e = float(data[10])
							if lineNumber==3:
								e_0 = e
							if lineNumber%3==1+shift and lineNumber>1:	#before encounter
								if len(data)<18:
									break
								de_tidal += e - e_prev
								if abs(de_tidal) > de_tidal_max:
									de_tidal_max = abs(de_tidal)
							if lineNumber%3==0+shift and lineNumber>3:	#after encounter
								de_flybys += e - e_prev
								if abs(de_flybys) > de_flybys_max:
									de_flybys_max = abs(de_flybys)
			print(index, de_tidal_max/(de_tidal_max+de_flybys_max), file=outfile)