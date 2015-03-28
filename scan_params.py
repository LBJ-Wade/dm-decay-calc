import numpy as np

"""Initialize scanning parameters"""
n_omegam = 25
n_omegar = 25
n_tau	 = 25

omegam_array=np.linspace(0.09, 0.3, n_omegam)
omegar_array=np.linspace(0.001, 0.09, n_omegar)
logtau_array=np.linspace(0, 30, n_tau)
