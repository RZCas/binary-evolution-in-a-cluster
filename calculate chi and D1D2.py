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

directories = ['storage/m1=m2=10']
for directory in directories:
	# types = ['mtotal=1e7,i=89.9,nokicks', 'mtotal=1e5,ns-ns', 'mtotal=1e6,ns-ns', 'mtotal=1e5,bh-ns', 'mtotal=1e6,bh-ns', 'mtotal=3e4', 'mtotal=2e5', 'mtotal=6e5', 'mtotal=1e5,ain=50', 'mtotal=1e5,ain=200', 'mtotal=1e6,ain=50']#os.listdir(directory)
	types1 = ['1e6,ns-ns,ain=200', '1e6,ns-ns,i=89.9', '1e7,ns-ns', '3e5,ns-ns', '3e6,ns-ns', '3e5_plummer', '3e6_plummer', '1e5,e=0', '1e6,e=0', '1e5,e=0.9', '1e6,e=0.9', '1e5,enc_only', '1e6,enc_only', '1e5,mper=0.5', '1e5,enc_only,mper=0.5']
	types = ['mtotal=' + type for type in types1]
	for type in types:
		if os.path.isfile(directory+'/'+type):	#if it's a file rather than a folder, ignore it
			continue
		root_dir = directory+'/'+type+'/'
		print(root_dir)
		outfile = open (root_dir + type + '-chiD1D2.txt', 'w+')
		print('index chi |D1|/(|D1|+D2) D1_flybys/D2_flybys', file=outfile)
		for index in range(0,1000):
			filepath = root_dir + str(index) + '.txt'
			if not os.path.isfile(filepath):
				continue
			lineNumber = 0
			shift = 0
			de_abs_tidal = 0
			de2_flybys = 0
			de_tidal = 0
			de_flybys = 0
			de_tidal_max = 0
			de_flybys_max = 0
			e = 0
			error = False
			with open(filepath) as f:
				lines = f.read().splitlines()
				if isfloat(lines[-1].split()[0]):
					t_final = 0.9*float(lines[-1].split()[0])
				else: #if the last line is "perturber:..."
					t_final = 0.9*float(lines[-2].split()[0])				
				for line in lines:
					lineNumber+=1 
					data = line.split()
					if len(data) > 1:
						if data[0]=="perturber:" and lineNumber%3==0:
							shift=1
						if isfloat(data[0]) and isfloat(data[1]):
							t = float(data[0])
							if t>t_final or len(data)<11:
								break
							e_prev = e
							if isfloat(data[10]):
								e = float(data[10])
							else:
								print(type, index, 'line', lineNumber, 'critical fucking error')
								error = True
								break
							if lineNumber==3:
								e_0 = e
							if lineNumber%3==1+shift and lineNumber>1:	#before encounter
								if len(data)<18:
									break
								de_tidal += e - e_prev
								de_abs_tidal += abs(e - e_prev)
								if abs(de_tidal) > de_tidal_max:
									de_tidal_max = abs(de_tidal)
							if lineNumber%3==0+shift and lineNumber>3:	#after encounter
								de_flybys += e - e_prev
								de2_flybys += (e - e_prev)**2
								if abs(de_flybys) > de_flybys_max:
									de_flybys_max = abs(de_flybys)
			if not error:
				# de_flybys -= e - e_prev	#to compensate for a possible very strong final encounter
				# de2_flybys -= (e - e_prev)**2
				chi = de_tidal_max / (de_tidal_max + de_flybys_max)
				D1D2 = de_abs_tidal / (de_abs_tidal + de2_flybys)
				print(index, chi, D1D2, de_flybys/de2_flybys, file=outfile)