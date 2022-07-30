import Mikkola
# from fortran.flybys3body_fortran import scattering
# import numpy as np
# from amuse.lab import units

N = 3
x = np.array([-1,0,0, 1,0,0, 10,0,0])
v = np.array([0,-1,0, 0,1,0, 0,0,5])
m = np.array([1,1,1])
time = np.array(0)
deltat = 0.001
eps=1.0e-19
newreg=1
ksmx=10000
soft=0
cmet=np.array([1,0.001,0])
clight=0
ixc=0
NBH=2
spin=np.array([0,0,0])
cmxx = np.array([0,0,0])
cmvx = np.array([0,0,0])
nmergers=np.array(0)
mergers=np.array([0,0,0])

fortran.Mikkola.chainevolve(N, x, v, m, time, deltat, eps, newreg, ksmx, soft, cmet, clight, ixc, NBH, spin,cmxx, cmvx, mergers, nmergers)
print(x)
print(v)
print(time)

# m1 = 1|units.MSun
# m2 = 1|units.MSun
# a = 1|units.AU
# e = 0.1
# i=0
# Omega=0
# omega=0
# meanAnomaly=0

# m3 = 1|units.MSun
# aStar = 1|units.AU
# eStar=0
# iStar=0
# OmegaStar=0
# omegaStar=0

# scattering (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar)