The function to be used is evolve_binary from binary_evolution_with_flybys.py; an example of its use is given in main.py 

Input parameters file structure:

1) integration time [yr] (shouldn't be much longer than 1e4 if you'd like the integration to be completed in seconds) 
2) outer orbit parameters: semimajor axis [pc], eccentricity (0...1), inclination (0...pi)
3) inner orbit parameters: primary mass [MSun], secondary mass [MSun], semimajor axis [AU] (shouldn't be much lower than 1), eccentricity (0...1), inclination (0...pi), argument of periapsis (0...2pi), longitude of ascending node (0...2pi)

Required python libraries:

amuse (https://amuse.readthedocs.io/en/latest/install/howto-install-AMUSE.html)
amuse-mikkola
orbitalpy
galpy