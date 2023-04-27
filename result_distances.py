import numpy as np
import pylab as plt
import distances
import densities
from scipy.interpolate import InterpolatedUnivariateSpline

H0 = 67.5
N_eff = 3.046
T_CMB = 2.7255
ombh2 = 0.022
omch2 = 0.122

rho_c = densities.critical_density(H0, ombh2)
rho_rad = densities.radiation_density(T_CMB, N_eff)
Omega_r = rho_rad/rho_c

ommh2 = ombh2 + omch2
Omega_m = ommh2 / (H0 /  100) ** 2

Omega_Lambda = 1 - Omega_r - Omega_m


z_star = 1100
chi_star = distances.comoving_distance(z_star, H0, Omega_r, Omega_m, Omega_Lambda)
print("----------")
print(f"comoving distance to LSS:  {chi_star} Mpc")
print("----------")

chi_star_giga_light_year = chi_star * 3.262e+6 / 10 ** 9
print("----------")
print(f"comoving distance to LSS:  {chi_star_giga_light_year} Gly")
print("----------")

logmin_z = -4
logmax_z = 4
nz = 10 ** 4

redshift = np.logspace(logmin_z, logmax_z, nz)
chi  = np.zeros(nz)
for i, z in enumerate(redshift):
    chi[i] = distances.comoving_distance(z, H0, Omega_r, Omega_m, Omega_Lambda)
    
plt.figure(figsize=(12,8))
plt.loglog()
plt.plot(redshift, chi)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("redshift", fontsize=18)
plt.ylabel("comoving distance (Mpc)", fontsize=18)
plt.savefig("comoving_distance.png",bbox_inches='tight')
plt.clf()
plt.close()

z_of_chi = InterpolatedUnivariateSpline(chi, redshift)
print("----------")
print(f"redshift for chi = 1 Gpc : {z_of_chi(1000)}")
print("----------")
