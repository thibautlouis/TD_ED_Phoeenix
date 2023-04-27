import numpy as np
from scipy.interpolate import interp2d

def get_P_z_and_k_interp():
    P_z_and_k = np.load("data/P_k_and_z.npy")
    k_array = np.load("data/k_array.npy")
    z_array = np.load("data/z_array.npy")
    P_z_and_k_interp = interp2d(z_array, k_array, P_z_and_k.flatten())
    return z_array, k_array, P_z_and_k_interp
