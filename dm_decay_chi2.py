# -*- coding: utf-8 -*-
"""
Constraints on the dark matter decay model.

TODO:
    1. Solve the problem of the large discrepancy in the orders

Created on Wed Jan 07 13:47:42 2015

@author: Hao Wang
@date: 2015-05-07
"""

'''Set parameter scope'''
import matplotlib
matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab.
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scan_params import *
import sys

'''
The only argument is choosing OHD('ohd') or SNe('sne')
'''
if len(sys.argv) != 2 or sys.argv[1] == '-h':
    print("[*]Usage: python dm_decay_chi2.py obs")
    print("[*]obs is one of the following numbers:")
    print("[*]0 - ohd")
    print("[*]1 - ohd+H0")
    print("[*]2 - sne")
    print("[*]3 - sne+ohd")
    print("[*]4 - sne+ohd+H0")
    print("[*]5 - sne+ohd+H0+Om")
    raise SystemExit

num_obs = int(sys.argv[1])
obs={0:'ohd', 1:'ohd+H0', 2:'sne', 3:'sne+ohd', 4:'sne+ohd+H0',
        5:'sne+ohd+H0+Om'}

z_ini = 1000.
omega_drI = 1.  # At large z_ini, no decay happened yet; nonzero as we use log.
omega_b0 = 0.02  # omega_b0 is defined wrt Hn=100; Omega_b0_LCDM*h^2=omega_b0 
omega_r0 = 2.47e-5  # The same as above.

'''OHD data'''
ohd_file = np.loadtxt('ohd.txt', usecols=range(0, 3))
# You must sort the data to get reasonable results.
ohd_file = ohd_file[ohd_file[:, 0].argsort()]
# Initial condition is at the largest redshift. So reverse the seq.
ohd_file = ohd_file[::-1]
ohd_data = ohd_file[:, 1]
ohd_error = ohd_file[:, 2]
ohd_z = ohd_file[:, 0]
ohd_x_output = np.log(1./(ohd_z+1.))

'''SNIa data'''
sne_file = np.loadtxt('SCPUnion2.1_mu_vs_z.txt', usecols=range(1, 4), skiprows=5)
sne_cov = np.loadtxt('SCPUnion2.1_covmat_sys.txt')
sne_inv_cov = np.linalg.inv(sne_cov)

# Sorting the data may cause inconsistency between data and invcov.
sne_file = sne_file[sne_file[:, 0].argsort()]
sne_file = sne_file[::-1]
sne_data = sne_file[:, 1]
sne_error = sne_file[:, 2]
sne_z = sne_file[:, 0]

'''constants'''
Hn = 100  # kmsmpc, used to normalize equations

'''
ODE of omega_dm & omega_dr evolution

n.b. The first element of x_output (or z) must correspond to the initial condition
of y_ini in the ode code -- it's imposed by odeint.

x    -> ln(a)
y[0] -> omega_dm
y[1] -> omega_dr

'''
def deriv(y, x, param):
   # w = 1./3
    tau = param[0]
    omega_lambda = param[1]
    omega_b0 = param[2]
    omega_r0 = param[3]
    omega_b = omega_b0*np.exp(-3*x)
    omega_r = omega_r0*np.exp(-4*x)
    if y[0]+y[1]+omega_lambda+omega_b+omega_r < 0:
        print y, x, tau
        print 'Total density is negative now. [Om, Or, Tau] = ', y[0], y[1], tau
        raise SystemExit

    E = np.sqrt(y[0]+y[1]+omega_lambda+omega_b+omega_r)
    du = -y[0]/(tau*E) - 3*y[0]
    dv = +y[0]/(tau*E) - 4*y[1]
    return [du, dv]

'''
ODE of ln(omega_dm) & ln(omega_dr) evolution
To solve the large order-of-magnitude difference accross the x-range.

x    -> ln(a)
y[0] -> ln(omega_dm)
y[1] -> ln(omega_dr)

'''
def deriv2(y, x, param):
    tau = param[0]
    omega_lambda = param[1]
    omega_b0 = param[2]
    omega_r0 = param[3]
    omega_b = omega_b0*np.exp(-3*x)
    omega_r = omega_r0*np.exp(-4*x)
    if np.exp(y[0])+np.exp(y[1])+omega_lambda+omega_b+omega_r < 0:
        print np.exp(y), np.exp(x), tau
        print 'Total density is negative now. [Om, Or, Tau] = ', y[0], y[1], tau
        raise SystemExit
    
    E = np.sqrt(np.exp(y[0])+np.exp(y[1])+omega_lambda+omega_b+omega_r)
    du = -1./(tau*E) - 3.
    dv = np.exp(y[0]-y[1])/(tau*E) - 4.
    return [du, dv]


'''
Solve ODE to get hubble parameter at any x==ln(a); n.b. H0 factor included 
It's a bit redundant now that the deriv function already got E calculated; but
deriv is designed to return only u&v. 
'''
def get_hubble(x, omega_dmI, omega_lambdaI, tau):
    #print 'get_hubble argument x: ', x 
    # x_output's first element must correspond to the initial y[:] value
    x_output = np.insert(x, 0, np.log(1./(1+z_ini)))
    #om_history = odeint(deriv, [omega_dmI, omega_drI], x_output, args=([tau,
    #    omega_lambda, omega_b0, omega_r0],))
    #omega_dmHistory = om_history[:, 0]
    #omega_drHistory = om_history[:, 1]
    om_history = odeint(deriv2, [np.log(omega_dmI), np.log(omega_drI)], 
            x_output, args=([tau, omega_lambdaI, omega_b0, omega_r0],))
    #print om_history
    omega_lambda = omega_lambdaI  # *np.ones(len(x_output))
    omega_b = omega_b0*np.exp(-3*x_output)
    omega_r = omega_r0*np.exp(-4*x_output)
    omega_dmHistory = np.exp(om_history[:, 0])
    omega_drHistory = np.exp(om_history[:, 1])
    return np.sqrt((omega_dmHistory+omega_drHistory+omega_b+omega_r
            )[1:]+omega_lambda)*Hn


'''Luminosity distance is the integration of 1/H over z'''
def integrand(z, omega_dmI, omega_lambdaI, tau):
    x = np.log(1.0/(1.0+z))
    return 1.0/get_hubble(x, omega_dmI, omega_lambdaI, tau)


'''Calculate D_L (in [Mpc]); used in the calculation of mu (for SNIa)'''
def get_dl(z, omega_dmI, omega_lambdaI, tau):
    c = 299792.46  # [km/s]
    return (1.0+z)*c*quad(integrand, 0., z, args=(omega_dmI, omega_lambdaI, tau))[0]


'''
Construct D_L array w.r.t observed redshifts
TODO:
    1. Get 'real' LCDM luminosity distance, instead of an approximation.
    2. Add Omega_dmI as a constraint. (can not be too small)
    3. Vary omega_drI and see its effect.

dl_theory = []
for z in sne_z:
    dl_theory.append(round(get_dl(z, 0.2, 0.01, 1000), 2))
'''


'''Omega_dm0 as an observable'''
def get_chi2_Omegadm0(omega_dmI, omega_lambdaI, tau):
    return 0.


'''Hubble constant as an observable'''
def get_chi2_hubble(omega_dmI, omega_lambdaI, tau):
    H0, H0_error = 73.8, 2.4  # [kmsmpc]
    x_H0 = 0.0  # x==ln(a(H0)=1)
    ohd_theory = get_hubble(x_H0, omega_dmI, omega_lambdaI, tau)
    chi2_hubble = np.power(ohd_theory[0] - H0, 2)/np.power(H0_error, 2)


'''For each parameter set (initial value), get chi2'''   
def get_chi2_ohd(omega_dmI, omega_lambdaI, tau):
    ohd_theory = get_hubble(ohd_x_output, omega_dmI, omega_lambdaI, tau)
    chi2_ohd = np.sum(np.power(ohd_theory-ohd_data, 2)/np.power(ohd_error, 2))
    return chi2_ohd


'''With system error + covariance'''
def get_chi2_sne2(omega_dmI, omega_lambdaI, tau):
    mu = []
    for z in sne_z:
        tmp = 5*np.log10(get_dl(z, omega_dmI, omega_lambdaI, tau))+25.0
        mu.append(tmp)
    chi2_sne=np.dot(mu-sne_data, np.dot(mu-sne_data, sne_inv_cov))
    return chi2_sne


'''No covariance information'''
def get_chi2_sne(omega_dmI, omega_lambdaI, tau):
    chi2_sne=0
    for iz, z in enumerate(sne_z):
        dl = get_dl(z, omega_dmI, omega_lambdaI, tau)
        mu = 5*np.log10(dl)+25.0
        chi2_sne += np.power((mu-sne_data[iz]), 2)/np.power(sne_error[iz], 2)
    return chi2_sne


'''Currently we are not adding up two chi2's, but plot them separately.'''
def get_chi2(omega_dmI, omega_lambdaI, tau):
    if num_obs == 0:
        return get_chi2_ohd(omega_dmI, omega_lambdaI, tau) 
    elif num_obs == 1:
        return get_chi2_sne2(omega_dmI, omega_lambdaI, tau)
    elif num_obs == 2:
        return get_chi2_sne2(omega_dmI, omega_lambdaI, tau) + \
                 get_chi2_ohd(omega_dmI, omega_lambdaI, tau)
    elif num_obs == 3:
        return get_chi2_sne2(omega_dmI, omega_lambdaI, tau) + \
                 get_chi2_ohd(omega_dmI, omega_lambdaI, tau) + \
                 get_chi2_hubble(omega_dmI, omega_lambdaI, tau)
    elif num_obs == 4:
        return get_chi2_sne2(omega_dmI, omega_lambdaI, tau) + \
                 get_chi2_ohd(omega_dmI, omega_lambdaI, tau) + \
                 get_chi2_hubble(omega_dmI, omega_lambdaI, tau) + \
                 get_chi2_Omegadm0(omega_dmI, omega_lambdaI, tau)
    else:
        print("No such observation probe available / Under construction.")
        raise SystemExit


'''Fill in the chi2 matrix'''
chi2 = np.zeros((n_omega_dm, n_omega_lambda, n_tau))
minchi2 = 100000
for iom, omega_dmI in enumerate(omega_dm_array):
    for iol, omega_lambdaI in enumerate(omega_lambda_array):
        for itau, logtau in enumerate(logtau_array):
            tau = np.exp(logtau)
            tmp = get_chi2(omega_dmI, omega_lambdaI, tau)
            print omega_dmI, omega_lambdaI, logtau, tmp
            if minchi2 > tmp:
                minchi2 = tmp
                # peakOm, peakOl, peakTau = omega_dmI, omega_lambdaI, tau

            chi2[iom, iol, itau] = tmp
       
print('minchi2 is ', minchi2)
fid="chi2file-"+obs[num_obs]+'-'+str(n_omega_dm)+'-'+str(n_omega_lambda)+'-'+str(n_tau)+'.bin'
chi2file = chi2.tofile(fid)
#np.savetxt('chi2-hr.txt', chi2, fmt="%.2f")

"""
minchi2 = np.amin(chi2)
dchi2 = chi2-minchi2

'''Integrate over the third (tau) dimension'''
chi2_omol = -2*np.log(np.sum(np.exp(-dchi2/2.), axis=2))
minchi2_omol = np.min(chi2_omol)
[peak_om_index, peak_ol_index] = np.unravel_index(np.argmin(chi2_omol), np.shape(chi2_omol))
peak_om = omega_dm_array[peak_om_index]
peak_ol = omega_lambda_array[peak_ol_index]
dchi2_omol = chi2_omol-minchi2_omol
plt.figure()
contourOmOl = plt.contour(omega_lambda_array, omega_dm_array, dchi2_omol, levels=[2.3, 5.0])
plt.plot(peak_ol, peak_om, 'r*')
# plt.clabel(contourOmOl, inline=1,  extent=(1, 3, 0, 2))
plt.xlabel(r'$\Omega_\Lambda$', fontsize=big)
plt.ylabel(r'$\Omega_m$', fontsize=big)
plt.set_yscale("log")
plt.savefig('ol-om-'+obs[num_obs]+'.eps')

'''Integrate over the second (omega_lambda) dimension'''
chi2_omtau = -2*np.log(np.sum(np.exp(-dchi2/2.), axis=1))
minchi2_omtau = np.min(chi2_omtau)
[peak_om_index, peak_logtau_index] = np.unravel_index(np.argmin(chi2_omtau), np.shape(chi2_omtau))
peak_om = omega_dm_array[peak_om_index]
peak_logtau = logtau_array[peak_logtau_index]
dchi2_omtau = chi2_omtau-minchi2_omtau
plt.figure()
contourOmTau = plt.contour(logtau_array, omega_dm_array, dchi2_omtau, levels=[2.3, 5.0])
plt.plot(peak_logtau, peak_om, 'r*')
# plt.clabel(contourOmTau, inline=1)
plt.xlabel(r'$\tau$', fontsize=big)
plt.ylabel(r'$\Omega_m$', fontsize=big)
plt.set_yscale("log")
plt.savefig('tau-om-'+obs[num_obs]+'.eps')
"""
