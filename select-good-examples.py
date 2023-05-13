import os
import glob
import numpy as np
from binary_evolution_with_flybys import a_h
from amuse.lab import units, constants 

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

# types = ['wide_range_1','wide_range_mtotal=1e5','wide_range_m=60','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e6_nokicks_plummer','wide_range_mtotal=1e6_nokicks']
types = ['wide_range_mtotal=1e5', 'wide_range_mtotal=1e5_plummer', 'wide_range_mtotal=1e5_plummer_b=1']
for type in types:
	root_dir = "storage/"+type+"/"
	outfile = open(root_dir+'good_examples.txt', 'w+')
	for index in range(0,100):
		filepath = root_dir + str(index) + '.txt'
		lineNumber=0
		with open(filepath) as f:
			for line in f:
				lineNumber+=1
				data = line.split()
				if lineNumber == 4:
					epsilon = float(data[17])
					if epsilon < 20:
						print(index, file=outfile)
						break