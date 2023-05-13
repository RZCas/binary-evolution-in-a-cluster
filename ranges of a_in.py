import os
import glob
import numpy as np

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

directories = ['storage']
print('ranges of a_in:')
for directory in directories:
	types = os.listdir(directory)
	for type in types:
		if os.path.isfile(directory+'/'+type) or type=='m1=m2=10':	#if it's a file rather than a folder, ignore it
			continue
		a_min = 0
		a_max = 0
		root_dir = directory+'/'+type+'/'
		for index in range(0,1000):
			filepath = root_dir + str(index) + '.txt'
			if not os.path.isfile(filepath):
				continue
			lineNumber=0
			with open(filepath) as f:
				for line in f:
					lineNumber+=1
					data = line.split()
					if lineNumber==3:
						a = float(data[7])
						if a_min == 0 or a<a_min:
							a_min = a
						if a>a_max:
							a_max = a
						break
		print(type+':', a_min, '-', a_max)

