import numpy as np
import pylab as plt
import distances
import densities
import lensing

H0 = 67.5
N_eff = 3.046
T_CMB = 2.7255
ombh2 = 0.022
omch2 = 0.122

rho_c = densities.critical_density(H0, ombh2)
rho_rad = densities.radiation_density(T_CMB, N_eff)
Omega_r = rho_rad / rho_c
ommh2 = ombh2 + omch2
Omega_m = ommh2 / (H0 /  100) ** 2
Omega_Lambda = 1 - Omega_r - Omega_m

z_of_chi = distances.interpolate_z_and_chi(H0, Omega_r, Omega_m, Omega_Lambda)

_, k, P_z_and_k_interp = lensing.get_P_z_and_k_interp("lcdm")

plt.figure(figsize=(12,8))
plt.loglog()
zplot = [0, 0.5, 1, 4 ,20]
for zp in zplot:
    plt.plot(k, P_z_and_k_interp(zp, k), label=zp)
plt.xlabel("k Mpc", fontsize=18)
plt.ylabel("$k^{4} P_\Psi\, Mpc^{-3}$", fontsize=18)
plt.legend(fontsize=18)
plt.show()

z_star = 1100

chi_star = distances.comoving_distance(z_star, H0, Omega_r, Omega_m, Omega_Lambda)

chi_array = np.linspace(10, chi_star, 100)
z_array = z_of_chi(chi_array)
dchis = chi_array[1] - chi_array[0]

l_array = np.arange(2, 2000 + 1, dtype=np.float64)
cl_phi = np.zeros(l_array.shape)
for i, l in enumerate(l_array):
    k_array = (l + 0.5) / chi_array
    cl_phi[i] = 0
    for z, chi, k in zip(z_array, chi_array, k_array):
        cl_phi[i] += 4 * dchis * P_z_and_k_interp(z, k) * ((chi_star - chi) / (chi ** 2 * chi_star)) ** 2 / k ** 4
    
fac = np.sqrt(l_array)*(l_array * (l_array + 1)) ** 2 / (2 * np.pi)

compare_with_camb = True
if compare_with_camb == True:
    import camb
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(ns=0.965)

    pars.set_for_lmax(2500,lens_potential_accuracy=2)
    results = camb.get_results(pars)
    cl_camb=results.get_lens_potential_cls(2500)
    ell_camb = np.arange(2, cl_camb[:,0].size)
    #plt.plot(cl_camb[2:,0]*np.sqrt(ell_camb), color='r')


#plt.loglog(l_array, cl_phi * fac, color='b')
plt.figure(figsize=(12,8))
plt.xlim([10,2000])
plt.semilogx()
lb, cb, sigma_b = np.loadtxt("data/result_act.dat", unpack=True)
plt.errorbar(lb, cb, sigma_b, fmt=".",color="red")
plt.plot(l_array, cl_phi * fac, color="black")
plt.xlabel(r"$\ell$", fontsize=22)
plt.ylabel(r"$\ell^{5/2}(\ell+1)^{2}C^{L}_{\ell}/(2\pi)$", fontsize=22)

plt.show()



