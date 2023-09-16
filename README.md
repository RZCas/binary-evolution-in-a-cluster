# BESC: Binary Evolution in a Star Cluster
 
The function to be used is evolve\_binary from binary\_evolution\_with\_flybys.py; an example of its use is given in main.py 

Input parameters file structure:

1. Integration time [yr] (shouldn't be much longer than 1e4 if you'd like the integration to be completed in seconds) 
2. Outer orbit parameters: semimajor axis [pc], eccentricity (0...1), inclination (0...pi)
3. Inner orbit parameters: primary mass [MSun], secondary mass [MSun], semimajor axis [AU] (shouldn't be much lower than 1), eccentricity (0...1), inclination (0...pi), argument of periapsis (0...2pi), longitude of ascending node (0...2pi)

Required python libraries:

* amuse (https://amuse.readthedocs.io/en/latest/install/howto-install-AMUSE.html)
* orbitalpy
* galpy