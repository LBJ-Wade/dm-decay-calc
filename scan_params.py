import numpy as np

"""Initialize scanning parameters"""
n_omega_dm = 5
n_omega_lambda  = 5
n_tau	   = 5

omega_dm_array=np.linspace(1e8, 1e10, n_omega_dm)
omega_lambda_array=np.linspace(0.1, 0.9, n_omega_lambda)
logtau_array=np.linspace(0, 30, n_tau)
