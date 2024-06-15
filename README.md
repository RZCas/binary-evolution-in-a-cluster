# BESC: Binary Evolution in a Star Cluster
 
The function that evolves the binary under the influence of both flybys and tides is `evolve_binary` from `binary_evolution_with_flybys.py`. For the input parameters definition, see the definition of the `inputParameters` class in that file. An example of its use is given in example.py 

## Installation

Required python libraries:

* `amuse` (https://amuse.readthedocs.io/en/latest/install/howto-install-AMUSE.html)
* `orbitalpy`
* `astropy`
* `galpy`

After the installation, run 
```
f2py mikk.pyf -c Mikkola.f -m Mikkola
```
in the `fortran` folder.

## Credits

Paper describing the code: https://arxiv.org/abs/2310.15374

Parts of the code were originally taken from:
* https://github.com/hamers/flybys
* https://github.com/mwbub/binary-evolution
* `AR-CHAIN` (https://ui.adsabs.harvard.edu/abs/2008AJ....135.2398M)