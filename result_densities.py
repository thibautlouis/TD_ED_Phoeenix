import densities

H0 = 67.5
N_eff = 3.046
T_CMB = 2.7255
ombh2 = 0.022
omch2 = 0.122

H0_s_minus1 = densities.H0_in_s_minus1(H0)

print("----------")
print(f"H0 in s^{-1} unit: {H0_s_minus1}")
print("----------")

rho_c = densities.critical_density(H0, ombh2)

print("----------")
print(f"rho_c in kg.m^{-3} unit:  {rho_c}")
print("----------")

rho_b = densities.omh2_to_density(H0, ombh2)

print("----------")
print(f"rho_b in kg.m^{-3} unit:  {rho_b}")
print("----------")

rho_cdm = densities.omh2_to_density(H0, omch2)

print("----------")
print(f"rho_cdm in kg.m^{-3} unit:  {rho_cdm}")
print("----------")

mH = 1.6735575 * 10 ** (-27) # hydrogen mass in kg
nH_per_m3 = (rho_b / mH)

print("----------")
print(f"number of hydrogen atom per m^{3}:  {nH_per_m3}")
print("----------")

rho_gamma = densities.radiation_density(T_CMB, 0)
Omega_gamma = rho_gamma/rho_c
print("----------")
print(f"Omega_gamma:  {Omega_gamma}")
print("----------")

rho_rad = densities.radiation_density(T_CMB, N_eff)
Omega_rad= rho_rad/rho_c

print("----------")
print(f"Omega_rad:  {Omega_rad}")
print("----------")
