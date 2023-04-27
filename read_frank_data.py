import numpy as np
import pylab as plt

bin_centers,auto_bandpowers,curl_bandpowers,diagonal_errors,pte,chi,covmat=np.load('./DR6_model_subtract270123.npy',allow_pickle=True)

scaling=(2/np.pi)*np.sqrt(bin_centers)
plt.semilogx()
plt.errorbar(bin_centers,auto_bandpowers*scaling,yerr=diagonal_errors[0]*scaling)
plt.show()

np.savetxt("result_act.dat", np.transpose([bin_centers, auto_bandpowers*scaling, diagonal_errors[0]*scaling]))
