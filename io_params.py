import numpy as np

# Need to be consistent with get_chi2().
obs_dict = {0:'ohd',
     1:'h0', 
     2:'om0', 
     3:'sne'}

# Color for contour level lines.
color_dict = {0:'g',
        1:'m',
        2:'y',
        3:'b',
        10:'r'}

big_font = 16
'''boundary condition'''
omega_drI = 1.  # At large z_ini, no decay happened yet; nonzero as we use log.
omega_b0 = 0 #0.02  # omega_b0 is defined wrt Hn=100 --> Omega_b0_LCDM*h^2=omega_b0 
omega_r0 = 0 #2.47e-5  # Added for formula precision; should be negligible.

'''scanning parameters'''
n_omega_dm = 30
n_omega_lambda = 30
n_tau	   = 30

# value at z=z_ini
zI = 1000.
omega_dm_array=np.linspace(1.0e7, 3.5e8, n_omega_dm)
# value at z=0
omega_lambda_array=np.linspace(0.1, 0.6, n_omega_lambda)
tau_array=np.logspace(-1.7, 3, n_tau)
