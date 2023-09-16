import os
import glob
import numpy as np

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

root_dir = "storage/m1=m2=10/mtotal=1e6/"
# outfile = open(root_dir+'ejected.txt', 'w+')
for index in range(0,1000):
	filepath = root_dir + str(index) + '.txt'
	with open(filepath) as f:
		for line in f:
			data = line.split()
			if len(data)>8 and isfloat(data[8]):
				m = float(data[8])
			if data[1] == 'ejected':
				if m==20:
					print(index, 'exchange+ejected')
				else:
					print(index, 'ejected')
