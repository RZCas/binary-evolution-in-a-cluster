matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
figure = pyplot.figure()

t_approx = []
t_precise = []
x_approx = []
x_precise = []
lineNumber = 0
with open('output/test.txt') as f:
for line in f:
	
	data = line.split()
	t

plot_a = figure.add_subplot(2,1,1)
ax = pyplot.gca()
ax.minorticks_on() 
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$\Delta a$ [AU]', fontsize=16)
# pyplot.xscale('log')
# pyplot.yscale('log')
# pyplot.text(0.3, 0.9, '$e='+str(e)+'$', fontsize=16, transform=ax.transAxes)
plot_a.plot(ts, k.a_array-a_in, 'k', label=r'GR only')
plot_a.plot(ts, k1.a_array-a_in, 'r', label=r'GR + tidal')
plot_a.legend(fontsize=16, frameon=False)
# plot.set_xlim(500, 5000)
# plot.set_ylim(1e-8, 1e-5)

plot_e = figure.add_subplot(2,1,2)
ax = pyplot.gca()
ax.minorticks_on() 
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$t$ [yr]', fontsize=16)
ax.set_ylabel(r'$\Delta e$', fontsize=16)
# pyplot.xscale('log')
# pyplot.yscale('log')
# pyplot.text(0.3, 0.9, '$e='+str(e)+'$', fontsize=16, transform=ax.transAxes)
plot_e.plot(ts, k.e_array-ecc, 'k', label=r'GR only')
plot_e.plot(ts, k1.e_array-ecc, 'r', label=r'GR + tidal')
plot_e.legend(fontsize=16, frameon=False)
# plot.set_xlim(500, 5000)
# plot.set_ylim(1e-8, 1e-5)

pyplot.tight_layout()
pyplot.savefig(input.output_file)