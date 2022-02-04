import matplotlib
from matplotlib import pyplot

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{physics}"
figure = pyplot.figure()

t_approx = []
t_precise = []
x_approx = []
x_precise = []
lineNumber = 0
with open('output/outer_orbit_r_ecc_out=0.99_n=100.txt') as f:
	for line in f:
		data = line.split()
		t_precise.append(float(data[0]))
		x_precise.append(float(data[1]))
with open('output/outer_orbit_r_ecc_out=0.99_n=10.txt') as f:
	for line in f:
		data = line.split()
		t_approx.append(float(data[0]))
		x_approx.append(float(data[1]))

plot_a = figure.add_subplot(1,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$r$ [pc]', fontsize=16)
# pyplot.xscale('log')
# pyplot.yscale('log')
# pyplot.text(0.3, 0.9, '$e='+str(e)+'$', fontsize=16, transform=ax.transAxes)
plot_a.plot(t_approx, x_approx, 'r', label=r'interpolation')
plot_a.plot(t_precise, x_precise, 'k', label=r'precise')
plot_a.legend(fontsize=16, frameon=False)
# plot.set_xlim(500, 5000)
# plot.set_ylim(1e-8, 1e-5)

pyplot.tight_layout()
pyplot.savefig('output/compare-orbits-r-ecc_out=0.99.pdf')