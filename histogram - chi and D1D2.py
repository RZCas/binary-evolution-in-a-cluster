import os
import glob
import shutil
import numpy as np
from statistics import stdev
from pathlib import Path

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
# matplotlib.rcParams["figure.titlesize"] = 4

input_files = []
output_files = []

for input_file in glob.glob('output/chi and D1D2/*-chiD1D2.txt'):
	input_files.append(os.path.basename(input_file))
	chi = []
	D1D2 = []
	D1flybysD2 = []
	parameters = ''
	output_file = ''
	subfolder = 'm=1.4-20'
	if 'plummer' in input_file:
		parameters += 'Plummer'
		output_file += 'plummer'
		# subfolder += 'plummer'
	else: 
		parameters += 'Hernquist'
		output_file += 'hernquist'
		# subfolder += 'hernquist'
	if 'total=1e5' in input_file:
		parameters += r', $M_{\rm total}=10^5M_\odot$'
		output_file += ',mtot=1e5'
		# subfolder += ',mtot=1e5'
	elif 'total=3e5' in input_file:
		parameters += r', $M_{\rm total}=3\dot10^5M_\odot$'
		output_file += ',mtot=3e5'
		# subfolder += ',mtot=1e5'
	elif 'total=3e6' in input_file:
		parameters += r', $M_{\rm total}=3\dot10^6M_\odot$'
		output_file += ',mtot=3e6'
	elif 'total=1e7' in input_file:
		parameters += r', $M_{\rm total}=10^7M_\odot$'
		output_file += ',mtot=1e7'
		# subfolder += ',mtot=1e6-1e7'
	elif 'total=3e4' in input_file:
		parameters += r', $M_{\rm total}=3\dot10^4M_\odot$'
		output_file += ',mtot=3e4'
	elif 'total=2e5' in input_file:
		parameters += r', $M_{\rm total}=2\dot10^5M_\odot$'
		output_file += ',mtot=2e5'
	elif 'total=6e5' in input_file:
		parameters += r', $M_{\rm total}=6\dot10^5M_\odot$'
		output_file += ',mtot=6e5'
	else:
		parameters += r', $M_{\rm total}=10^6M_\odot$'
		output_file += ',mtot=1e6'
		# subfolder += ',mtot=1e6-1e7'
	if 'perpendicular' in input_file or 'b=1' in input_file:
		parameters += r', $b=1$ pc'
		output_file += ',b=1'
	else:
		parameters += r', $b=2$ pc'
		output_file += ',b=2'
	if 'perpendicular' in input_file:
		parameters += r', $a_{\rm out}=1.6$ pc'
		output_file += ',aout=1.6'
	elif 'a_out=1' in input_file or 'aout=1' in input_file:
		parameters += r', $a_{\rm out}=1$ pc'
		output_file += ',aout=1'
	elif 'a_out=3' in input_file or 'aout=3' in input_file:
		parameters += r', $a_{\rm out}=3$ pc'
		output_file += ',aout=3'
	elif 'a_out=4' in input_file or 'aout=4' in input_file:
		parameters += r', $a_{\rm out}=4$ pc'
		output_file += ',aout=4'
	else:
		parameters += r', $a_{\rm out}=2$ pc'
		output_file += ',aout=2'
	parameters += ',\n'		
	if 'ns' in input_file and not 'ns-ns' in input_file and not 'bh-ns' in input_file:
		parameters += r'$m_1=1.4M_\odot$, $m_2=1.6M_\odot$'
		output_file += ',m=1.4'
	elif 'ns-ns' in input_file:
		parameters += r'$m_1=m_2=1.4M_\odot$'
		output_file += ',m=1.4'
	elif 'bh-ns' in input_file:
		parameters += r'$m_1=10M_\odot$, $m_2=1.4M_\odot$'
		output_file += ',m1=10,m2=1.4'
	elif 'hernquist,' in input_file or os.path.basename(input_file)[:6] == 'mtotal':
		parameters += r'$m_1=m_2=10M_\odot$'
		output_file += ',m=10'
		subfolder = 'm1=m2=10'
	elif 'm=60' in input_file:
		parameters += r'$m_1=m_2=30M_\odot$'
		output_file += ',m=30'
	elif 'light' in input_file:
		parameters += r'$m_1=m_2=0.1M_\odot$'
		output_file += ',m=0.1'
	else:	
		parameters += r'$m_{1,2}=1.4\dots20M_\odot$'
		output_file += ',m=1.4-20'
	if 'hard' in input_file:
		parameters += r', $a_{\rm in}=3\dots30$ AU'
		output_file += ',ain=3-30'
	elif 'a_in=300' in input_file:
		parameters += r', $a_{\rm in}=300$ AU'
		output_file += ',ain=300'
	elif 'ain=50' in input_file:
		parameters += r', $a_{\rm in}=50$ AU'
		output_file += ',ain=50'
	elif 'ain=200' in input_file:
		parameters += r', $a_{\rm in}=200$ AU'
		output_file += ',ain=200'
	elif 'a_in=300' in input_file:
		parameters += r', $a_{\rm in}=300$ AU'
		output_file += ',ain=300'
	elif 'hernquist,' in input_file or os.path.basename(input_file)[:6] == 'mtotal':
		parameters += r', $a_{\rm in}=100$ AU'
		output_file += ',ain=100'
	elif 'uniform' in input_file or ('wide_range_mtotal' in input_file and not 'wide_range_mtotal=1e6_nokicks' in input_file):
		parameters += r', $a_{\rm in}=50\dots200$ AU'
		output_file += ',ain=50-200'
	else:
		parameters += r', $a_{\rm in}=2\dots200$ AU'
		output_file += ',ain=2-200'
	parameters += r', $e_{\rm in}=0.5$'
	if 'uniform' in input_file or ('hernquist,' in input_file and not 'i=89.9' in input_file) or (os.path.basename(input_file)[:6] == 'mtotal' and not 'i=89.9' in input_file):
		parameters += r', $i=0\dots180^\circ$'
		output_file += ',i=0-180'
	else:
		parameters += r', $i=89.9^\circ$'
		output_file += ',i=89.9'
	if 'nokick' in input_file:
		parameters += '\nno kicks'
		output_file += ',nokicks'
	if 'weakonly' in input_file:
		parameters += ', distant flybys only'
		output_file += ',weakonly'
	if 'noweak' in input_file:
		parameters += ', close flybys only'
		output_file += ',noweak'
	if 'm_per=0.1' in input_file:
		parameters += '\n'+r'$m_{\rm per}=0.1M_\odot$'
		output_file += ',mper=0.1'

	new_txt_file = 'output/chi and D1D2/' + output_file + '.txt'
	shutil.copy(input_file, new_txt_file)

	output_file += '.pdf'
	subfolder += '/'
	Path('output/chi and D1D2/'+subfolder).mkdir(parents=True, exist_ok=True)
	Path(str(Path.home())+'/Dropbox/CLUSTERS/chi and D1D2/'+subfolder).mkdir(parents=True, exist_ok=True)
	output_files.append(output_file[:-4])

	if os.path.basename(input_file)[:6] == 'mtotal':
		evolution_folder = 'output/m1=m2=10/'+os.path.basename(input_file)[:-12]
	else:
		evolution_folder = 'output/'+os.path.basename(input_file)[:-12]
	if os.path.exists(evolution_folder):
		target_folder = 'output/chi and D1D2/'+subfolder+output_file[:-4]
		target_folder_2 = str(Path.home())+'/Dropbox/CLUSTERS/chi and D1D2/'+subfolder+output_file[:-4]
		Path(target_folder).mkdir(parents=True,exist_ok=True)
		Path(target_folder_2).mkdir(parents=True,exist_ok=True)
		for evolution_plot in glob.glob(evolution_folder+'/*evolution*.pdf'):
			shutil.copy(evolution_plot, target_folder)
			shutil.copy(evolution_plot, target_folder_2)
	else:
		print('missing evolution plots:', evolution_folder)

	with open(input_file) as f:
		for line in f:
			data = line.split()
			if isfloat(data[1]):
				chi.append(float(data[1]))
				D1D2.append(float(data[2]))
				D1flybysD2.append(float(data[3]))
				# if abs(D1flybysD2[-1])>100:
				# 	print(input_file, data)

	if len(chi) == 0: continue
	chi = np.array(chi)
	D1D2 = np.array(D1D2)
	D1flybysD2 = np.array(D1flybysD2)

	figure = pyplot.figure(figsize=(9, 7/2))
	figure.suptitle(parameters)
	plot_tidal = figure.add_subplot(1,2,1)
	plot_tidal.set_title(rf'$\chi={chi.mean():.2f}\pm{stdev(chi):.2f}$'+'\n'+rf'$D_1/D_2={D1D2.mean():.2f}\pm{stdev(D1D2):.2f}$', x=0.5, y=0.75)
	plot_tidal.hist (chi, histtype='step', bins=np.linspace(0, 1, 11), label=r'$\chi$')
	plot_tidal.hist (D1D2, histtype='step', bins=np.linspace(0, 1, 11), label=r'$D_1/D_2$')
	plot_tidal.legend()#loc='lower left')

	# figure = pyplot.figure(figsize=(9/2, 7/2))
	plot_flybys = figure.add_subplot(1,2,2)
	plot_flybys.set_title(fr'$\langle D_{{\rm1,flybys}}/D_{{\rm 2,flybys}}\rangle = {D1flybysD2.mean():.3f}$', x=0.5, y=0.75)
	plot_flybys.hist (D1flybysD2, histtype='step', bins=np.linspace(-3, 3, 11), label=r'$D_{\rm1,flybys}/D_{\rm 2,flybys}$')
	plot_flybys.legend()#loc='lower left')
	# pyplot.tight_layout()
	# pyplot.savefig('output/chi and D1D2/'+subfolder+output_file[:-4]+'-D1flybys.pdf')
	# pyplot.savefig(str(Path.home())+'/Dropbox/CLUSTERS/chi and D1D2/'+subfolder+output_file[:-4]+'-D1flybys.pdf')
	# pyplot.clf()

	pyplot.tight_layout()
	pyplot.savefig('output/chi and D1D2/'+subfolder+output_file)
	pyplot.savefig(str(Path.home())+'/Dropbox/CLUSTERS/chi and D1D2/'+subfolder+output_file)
	pyplot.clf()

n = len(input_files)
for i in range(n-1):
	for j in range(i+1,n):
		if output_files[i] == output_files[j]:
			print(input_files[i]+' and '+input_files[j]+' both have '+output_files[i])