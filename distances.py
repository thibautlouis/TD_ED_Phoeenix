import numpy as np
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline

def comoving_distance(z, H0, Omega_r, Omega_m, Omega_Lambda):
    """
    This function compute the comoving distance between us at the origin of the
    coordinate system and an object at redshift z, it return its value in Mpc
    
    Parameters
     ----------
    H0: float
        Hubble constant value in km/s/Mpc
    Omega_r: float
        Omega radiation
    Omega_m: float
        Omega matter
    Omega_lambda: float
        Omega dark energy

    """
    c = 2.99792458e8 # m/s
    H0 *= 1000 # m/s/Mpc
    
    def integrand(z):
        H = H0 * np.sqrt(Omega_r * (1 + z) ** 4  + Omega_m * (1 + z) ** 3 + Omega_Lambda)
        return c / H

    chi, _ = quad(integrand, 0, z)
    
    return chi
    
    
def interpolate_z_and_chi(H0, Omega_r, Omega_m, Omega_Lambda, logmin_z=-4, logmax_z=4, nz=10**4):

    """
    This function return the interpolation of z(chi)
    
    Parameters
     ----------
    H0: float
        Hubble constant value in km/s/Mpc
    Omega_r: float
        Omega radiation
    Omega_m: float
        Omega matter
    Omega_lambda: float
        Omega dark energy
    logmin_z: float
        Minimal value of log(z) for the interpolation
    logmax_z: float
        Maximum value of log(z) for the interpolation
    nz: integer
        number of z values
    """

    redshift = np.logspace(logmin_z, logmax_z, nz)
    chi  = np.zeros(nz)
    for i, z in enumerate(redshift):
        chi[i] = comoving_distance(z, H0, Omega_r, Omega_m, Omega_Lambda)
    z_of_chi = InterpolatedUnivariateSpline(chi, redshift)

    return z_of_chi
