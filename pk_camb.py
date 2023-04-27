import camb
from camb import model
import numpy as np

H0 = 67.5
ombh2 = 0.022
omch2 = 0.122
ns = 0.965

pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
pars.InitPower.set_params(ns=ns)
kmax = 100  #kmax to use

logmin_z, logmax_z, nz = -3, 3, 10 ** 3
z_array = np.logspace(logmin_z, logmax_z, nz)
logmin_k, logmax_k, nk = -3, 2, 10 ** 4
k_array = np.logspace(logmin_k, logmax_k, nk)

PK = camb.get_matter_power_interpolator(pars, nonlinear=True,
                                        hubble_units=False, k_hunit=False, kmax=kmax,
                                        var1=model.Transfer_Weyl, var2=model.Transfer_Weyl,
                                        zmax=z_array[-1])


Pk_array = np.zeros((nz, nk))
k_array = np.logspace(logmin_k, logmax_k, nk)
for i, k in enumerate(k_array):
    Pk_array[:,i] = PK.P(z_array, k, grid=False) #come from camb convention psi is k^2 psi


np.save("data/P_k_and_z.npy", Pk_array)
np.save("data/k_array.npy", k_array)
np.save("data/z_array.npy", z_array)


