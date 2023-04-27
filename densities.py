import numpy as np

def H0_in_s_minus1(H0):
    """
    This function compute H0 in s^{-1} unit
    
    Parameters
     ----------
    H0: float
        Hubble constant value in km/s/Mpc

    """
    Mpc_to_km = 3.086 * 10 ** 19
    H0 /= Mpc_to_km
    
    return H0
    
def critical_density(H0, unit="kg/m^3"):
    """
    This function compute the critical density in kg/m^{3}
    optionnaly you can compute it in unit of solarmass/Mpc^3
    
    Parameters
     ----------
    H0: float
        Hubble constant value in km/s/Mpc
    """

    H0 = H0_in_s_minus1(H0)
    G = 6.67430e-11 # m^3/kg/s^2
    rho_c = (3 * H0 ** 2) / (8 * np.pi * G) #kg/m^{3}
    
    if unit == "solarmass/Mpc^3":
        # Convert kg/m^3 to solar mass/Mpc^3
        Msun = 1.989e30 # kg
        Mpc_to_m = 3.086 * 10 ** 22
        rho_c *= Mpc_to_m ** 3 / Msun

    return rho_c

def omh2_to_density(H0, omh2):
    """
    This function convert Omegah2 value into density
    
    Parameters
     ----------
    Omegah2: float
        std convention for matter density
    """
    h = H0 / 100
    rho_c = critical_density(H0)
    rho = omh2 * rho_c / h ** 2
    return rho

def get_Omega_Lambda(rho_r, ombh2, omch2, H0):
    """
    This function compute omega_lambda from other density parameter
    assuming the universe is flat
    
    Parameters
     ----------
    omr: float
        omega radiation
    ombh2: float
        omega baryons x h^2
    omch2: float
        omega cold dark matter x h^2
    H0: float
        hubble constant
    """

    h = H0 / 100
    rho_c = critical_density(H0)
    om_m = (ombh2 + omch2) / h ** 2
    om_lambda = 1 - om_m - (rho_r/ rho_c)
    return om_lambda
    
def radiation_density(T_CMB, N_eff):
    """
    This function compute the radiation density in kg/m^{3}
    it takes as input the CMB temperature in K
    the radiation is the sum of the photons density and neutrino density
    N_eff is the effective number of xtra relativistic degree of freedom
    
    Parameters
     ----------
    H0: float
        Hubble constant value in km/s/Mpc
    Neff: effective number of degree of freedom
    
    """
    rho_gamma_c2 = (np.pi ** 2) / 15 * T_CMB ** 4
    
    kB = 1.380649e-23 # J/K
    hbar = 6.62607015e-34 / (2 * np.pi) # J s
    c = 2.99792458e8 # m/s

    rho_gamma_c2 *= kB ** 4 /(hbar * c) ** 3 # in J / m^{3} (J = kg * m^{2} * s^{-2})
    rho_gamma = rho_gamma_c2 / c ** 2 # radiation density in kg/m^{3}
    
    # add neutrino contribution
    rho_gamma += N_eff * 7/8 * (4/11) ** (4/3) * rho_gamma
        
    return rho_gamma
    

