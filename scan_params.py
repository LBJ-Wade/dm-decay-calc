import numpy as np

"""Initialize scanning parameters"""
n_omega_dm = 25
n_omega_sm = 25
n_tau	   = 25

# value at z=z_ini
#omega_dm_array=np.linspace(1e8, 1e10, n_omega_dm)
omega_dm_array=np.logspace(7.7, 8.3, n_omega_dm)
# value at z=0
omega_sm_array=np.linspace(0.01, 0.5, n_omega_lambda)
tau_array=np.logspace(-2, 1, n_tau)
