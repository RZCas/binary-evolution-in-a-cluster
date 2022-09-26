# plot the histogram of |de_tidal|/|de_encounters| for a certain cluster parameters
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

types = ['wide_range_1','wide_range_mtotal=1e5','wide_range_m=60','wide_range_mtotal=1e5_nokicks','wide_range_mtotal=1e6_nokicks_plummer','wide_range_mtotal=1e6_nokicks','wide_range_mtotal=1e5_plummer','wide_range_mtotal=1e6_plummer','wide_range_mtotal=1e5_plummer_b=1','wide_range_mtotal=1e7_hernquist']

for type in types[6:]:
	root_dir = "storage/"+type+"/"
	output = open("storage/time-"+type+".txt", "w+")
	print("N model_time[Gyr] time_3body[s] time_tidal[s] 3body_time_fraction", file=output)
	index = 0
	while True:	#for every file in the folder
		shift = 0
		filepath = root_dir + str(index) + '.txt'
		if not os.path.isfile(filepath): 
			break
		lineNumber = 0
		time_3body = 0
		time_tidal = 0
		with open(filepath) as f:
			for line in f:
				lineNumber+=1
				data = line.split()
				if len(data) > 1:
					if isfloat(data[0]):
						model_time = float(data[0])
					if lineNumber%3==1+shift and lineNumber>=4 and isfloat(data[1]):	#before the encounter
						time_tidal += float(data[18])
					if lineNumber%3==0+shift and lineNumber>=6 and isfloat(data[1]):	#after the encounter
						if data[0]=='perturber:':
							shift = 1
						else:
							time_3body += float(data[14])
		print(index, model_time/1e9, time_3body, time_tidal, time_3body/(time_3body+time_tidal), file=output)
		index += 1
	output.close()
