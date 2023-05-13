# import fortran.Mikkola
from fortran.flybys3body_fortran import scattering as scattering_new
from flybys3body import scattering as scattering_old
import numpy as np
from amuse.lab import units
import time

# N = 3
# x = np.array([-1,0,0, 1,0,0, 10,0,0])
# v = np.array([0,-1,0, 0,1,0, 0,0,5])
# m = np.array([1,1,1])
# time = np.array(0)
# deltat = 0.001
# eps=1.0e-10
# newreg=1
# ksmx=10000
# soft=0
# cmet=np.array([1,0.001,0])
# clight=0
# ixc=0
# NBH=2
# spin=np.array([0,0,0])
# cmxx = np.array([0,0,0])
# cmvx = np.array([0,0,0])
# nmergers=np.array(0)
# mergers=np.array([0,0,0])

# fortran.Mikkola.chainevolve(N, x, v, m, time, deltat, eps, newreg, ksmx, soft, cmet, clight, ixc, NBH, spin,cmxx, cmvx, mergers, nmergers)
# print(x)
# print(v)
# print(time)

m1 = 10|units.MSun
m2 = 10|units.MSun
a = 1|units.AU
e = 0.8
i=0.1
Omega=0.1
omega=0.8
meanAnomaly=0.1

m3 = 1|units.MSun
aStar = -0.1|units.AU
eStar=7
iStar=0.4
OmegaStar=0.6
omegaStar=0.1

time0 = time.time()
result, third_body, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, shit = scattering_old (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar)
print("old:", result, third_body, a_new.value_in(units.AU), e_new, i_new, Omega_new, omega_new)
print("time:", time.time()-time0)

time0 = time.time()
result, third_body, finalTime, dv_binary, a_new, e_new, i_new, Omega_new, omega_new, aStar_new, eStar_new, iStar_new, OmegaStar_new, omegaStar_new, shit = scattering_new (m1, m2, a, e, i, Omega, omega, meanAnomaly, m3, aStar, eStar, iStar, OmegaStar, omegaStar)
print("new:", result, third_body, a_new.value_in(units.AU), e_new, i_new, Omega_new, omega_new)
print("time:", time.time()-time0)