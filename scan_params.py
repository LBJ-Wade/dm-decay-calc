import numpy as np

"""Initialize scanning parameters"""
n_omega_dm = 25
n_omega_lambda  = 25
n_tau	   = 25

# value at z=z_ini
#omega_dm_array=np.linspace(1e8, 1e10, n_omega_dm)
omega_dm_array=np.logspace(5, 8, n_omega_dm)
# value at z=0
omega_lambda_array=np.linspace(0.1, 0.9, n_omega_lambda)
logtau_array=np.linspace(0.03, 30, n_tau)
