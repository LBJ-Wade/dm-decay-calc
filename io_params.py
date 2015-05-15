import numpy as np

# Remember to update get_chi2 after changing this.
obs={0:'ohd', 
     1:'ohd+H0', 
     2:'ohd+H0+Om0', 
     3:'Om0',
     4:'H0',
     12:'sne', 
     13:'sne+ohd', 
     14:'sne+ohd+H0',
     15:'sne+ohd+H0+Om'}

'''boundary condition'''
omega_drI = 1.  # At large z_ini, no decay happened yet; nonzero as we use log.
omega_b0 = 0 #0.02  # omega_b0 is defined wrt Hn=100 --> Omega_b0_LCDM*h^2=omega_b0 
omega_r0 = 0 #2.47e-5  # Added for formula precision; should be negligible.

'''scanning parameters'''
n_omega_dm = 5 
n_omega_lambda = 5
n_tau	   = 5

# value at z=z_ini
zI = 1000.
omega_dm_array=np.linspace(1e8, 5e8, n_omega_dm)
# value at z=0
omega_lambda_array=np.linspace(0.1, 0.9, n_omega_lambda)
tau_array=np.logspace(-2, 3, n_tau)
