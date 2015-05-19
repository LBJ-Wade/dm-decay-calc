# -*- coding: utf-8 -*-
"""
Constraints on the dark matter decay model.
Functionality:
    1. --- 

Usage:
    1. "python dm_decay_chi2.py -h" for usage.
    2. E.g., "python dm_decay_chi2.py 0" calculates chi2 for OHD data.
    3. The resulted chi2 values are stored in .bin files.
    4. These .bin files are processed by contour_plot.py to get contour plots.
    5. E.g., "python contour_plot.py 0,1" plots contours based on OHD&h0 data.

TODO:
    1. Solve the problem of the large discrepancy in the orders

Created on Wed Jan 07 13:47:42 2015

@author: Hao Wang
@date: 2015-05-07
"""

'''Set parameter scope'''
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from io_params import *
import sys

'''
The only argument is choosing OHD('ohd') or SNe('sne')
'''
if len(sys.argv) != 2 or sys.argv[1] == '-h':
    print("[*]Usage: python dm_decay_chi2.py obs")
    print("[*]obs is one of the following numbers:")
    print(obs)
    raise SystemExit

num_obs = int(sys.argv[1])
print('Using', obs[num_obs], 'as observable.')

'''OHD data'''
ohd_file = np.loadtxt('ohd.txt', usecols=range(0, 3))
# Must sort the data w.r.t. redshifts to get reasonable results. 
ohd_file = ohd_file[ohd_file[:, 0].argsort()]
# Initial condition is at the largest redshift. So reverse the seq.
ohd_file = ohd_file[::-1]
ohd_z = ohd_file[:, 0]
ohd_data = ohd_file[:, 1]
ohd_error = ohd_file[:, 2]
ohd_x = np.log(1./(ohd_z+1.))

'''SNIa data'''
sne_file = np.loadtxt('SCPUnion2.1_mu_vs_z.txt', usecols=range(1, 4), skiprows=5)
#sne_cov = np.loadtxt('SCPUnion2.1_covmat_sys.txt')
#sne_inv_cov = np.linalg.inv(sne_cov)

# Sorting the data may cause inconsistency between data and invcov.
sne_file = sne_file[sne_file[:, 0].argsort()]
sne_file = sne_file[::-1]
sne_z = sne_file[:, 0]
sne_data = sne_file[:, 1]
sne_error = sne_file[:, 2]

'''constants'''
Hn = 100  # kmsmpc, used to normalize equations
c = 299792.46  # [km/s]

'''
ODE of ln(omega_dm) & ln(omega_dr) evolution.
Linear variable may cause precision problem as the equation involves up to a z**3
times (10**(9~12) for z=1000) difference.

x    -> ln(a)
y[0] -> ln(omega_dm)
y[1] -> ln(omega_dr)

'''
def deriv(y, x, param):
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
    du = - 3. - 1./(tau*E) 
    dv = - 4. + np.exp(y[0]-y[1])/(tau*E) 
    return [du, dv]


'''
Solve ODE to get density parameter at any x==ln(a). The initial results are
logrithm of the density parameter, so np.exp() them.

n.b. The first element of x_output (or z) must correspond to the initial
condition of y_ini in the ode code -- it's imposed by odeint. So insert a x_ini
to the targeted x array.
'''
def get_omega(x, zI, omega_dmI, omega_lambdaI, tau):
    x=np.array(x)
    x_output = np.insert(x, 0, np.log(1./(1+zI)))
    #print 'get_omega:', z_ini, 'x=', x, 'x_output=', x_output
    #raw_input()
    return np.exp(odeint(deriv, [np.log(omega_dmI), np.log(omega_drI)], 
            x_output, args=([tau, omega_lambdaI, omega_b0, omega_r0],)))


def get_hubble(x, zI, omega_dmI, omega_lambdaI, tau):
    #print 'get_hubble argument x: ', x 
    omega_lambda = omega_lambdaI  
    omega_b = omega_b0*np.exp(-3.0*x)
    omega_r = omega_r0*np.exp(-4.0*x)
    om_history = get_omega(x, zI, omega_dmI, omega_lambdaI, tau)
    omega_dm = om_history[1:, 0]  # Rip off x_ini.
    omega_dr = om_history[1:, 1]
    # print 'get_hubble:', om_history, omega_dm, omega_dr, omega_b, omega_r
    return np.sqrt(omega_dm+omega_dr+omega_b+omega_r+omega_lambda)*Hn


'''Luminosity distance is the integration of 1/H over z'''
def integrand(z, zI, omega_dmI, omega_lambdaI, tau):
    x = np.log(1.0/(1.0+z))
    return 1.0/get_hubble(x, zI, omega_dmI, omega_lambdaI, tau)

'''Calculate D_L (in [Mpc]); used in the calculation of mu (for SNIa)'''
def get_dl(z, zI, omega_dmI, omega_lambdaI, tau):
    return (1.0+z)*c*quad(integrand, 0., z, args=(zI, omega_dmI, omega_lambdaI,
        tau))[0]

def get_dl_union2(z_union2, zI, omega_dmI, omega_lambdaI, tau):
    z_list = np.logspace(np.log10(min(sne_z)), np.log10(max(sne_z)), 20)
    dl_list = []
    for z in z_list:
        dl_list.append(get_dl(z, zI, omega_dmI, omega_lambdaI, tau))

    from scipy import interpolate
    func = interpolate.InterpolatedUnivariateSpline(z_list, dl_list)
    return func(sne_z)


def integrand_lcdm(z, omega_dm0, omega_lambda0):
    H0 =70.
    return 1.0/(np.sqrt(omega_dm0*(1+z)**3 + omega_lambda0 + 
            omega_b0*(1+z)**3 + omega_r0*(1+z)**4) * H0)

def get_dl_lcdm(z, omega_dm0, omega_lambda0):
    return (1.0+z)*c*quad(integrand_lcdm, 0., z, args=(omega_dm0, 
        omega_lambda0))[0]


'''
Construct D_L array w.r.t observed redshifts
TODO:
    1. Vary omega_drI to see its effects.

'''


'''
Omega_dm0 as a single-point observable; James et. al., Deep Lens (1210.2732)
n.b.: omega_dmI*1^2 = Omega_dm0*h^2, due to a different normalization.
om0_error is used directly to get a conservative estimate (James et.al. got 
Omega_dm0 constraints after marginalizing over h; so it's really an approx.
here.
'''
def get_chi2_omegam0(zI, omega_dmI, omega_lambdaI, tau):
    om0, om0_error = 0.262*0.722**2, 0.05
    x_0 = 0.0  # x==ln(a(Om0)=1) = 0
    # get_omega returns [[omega_dm0, omega_dr0], [omega_dm1, omega_dr1]]
    omega_dm0 = get_omega(x_0, zI, omega_dmI, omega_lambdaI, tau)[1, 0]
    # print 'om-params', omega_dmI, omega_lambdaI, tau, omega_dm0
    return np.power(omega_dm0+omega_b0-om0, 2)/np.power(om0_error, 2)


'''Hubble constant as a single-point observable; Riess 3% H0 (1103.2976)'''
def get_chi2_hubble(zI, omega_dmI, omega_lambdaI, tau):
    H0, H0_error = 73.8, 2.4  # [kmsmpc]
    x_0 = 0.0  # x==ln(a(H0)=1) = 0
    ohd_theory = get_hubble(x_0, zI, omega_dmI, omega_lambdaI, tau)
    return np.power(ohd_theory[0] - H0, 2)/np.power(H0_error, 2)


'''For each parameter set (initial value), get chi2'''   
def get_chi2_ohd(zI, omega_dmI, omega_lambdaI, tau):
    ohd_theory = get_hubble(ohd_x, zI, omega_dmI, omega_lambdaI, tau)
    return np.sum(np.power(ohd_theory-ohd_data, 2)/np.power(ohd_error, 2))


'''Interpolation.'''
def get_chi2_sne(zI, omega_dmI, oemga_lambdaI, tau):
    chi2_sne = 0
    dl = get_dl_union2(sne_z, zI, omega_dmI, omega_lambdaI, tau)
    mu = 5*np.log10(dl)+25.0
    chi2_sne = np.sum(np.power(mu-sne_data, 2)/np.power(sne_error, 2))
    return chi2_sne

'''Straightforward calculation. No covariance information'''
def get_chi2_sne1(zI, omega_dmI, omega_lambdaI, tau):
    chi2_sne=0
    for iz, z in enumerate(sne_z):
        dl = get_dl(z, zI, omega_dmI, omega_lambdaI, tau)
        mu = 5*np.log10(dl)+25.0
        chi2_sne += np.power((mu-sne_data[iz]), 2)/np.power(sne_error[iz], 2)
    return chi2_sne

'''
With system error + covariance; 
NOT usable as sne_z is in different order than in sne_inv_cov.
Need to re-order mu before using this method.
'''
def get_chi2_sne2(zI, omega_dmI, omega_lambdaI, tau):
    mu = []
    for z in sne_z:
        tmp = 5*np.log10(get_dl(z, zI, omega_dmI, omega_lambdaI, tau))+25.0
        mu.append(tmp)
    chi2_sne=np.dot(mu-sne_data, np.dot(mu-sne_data, sne_inv_cov))
    return chi2_sne


def get_chi2(num_obs, zI, omega_dmI, omega_lambdaI, tau):

    if num_obs == 0:
        return get_chi2_ohd(zI, omega_dmI, omega_lambdaI, tau)
    elif num_obs == 1:
        return get_chi2_hubble(zI, omega_dmI, omega_lambdaI, tau)
    elif num_obs == 2:
        return get_chi2_omegam0(zI, omega_dmI, omega_lambdaI, tau)
    elif num_obs == 3:
        return get_chi2_sne(zI, omega_dmI, omega_lambdaI, tau)
    else:
        print("No such observation probe available yet.")
        raise SystemExit


'''TEST MODULE'''
import time
class test:
    h = 0.7
    omega_m0 = 0.3
    omega_l0 = 0.7
    omega_mI = omega_m0*h**2 * (1+zI)**3
    omega_lI = omega_l0*h**2 

    def __init__(self, zI, tau):
        self.zI = zI
        self.tau= tau

    def __repr__(self):
        return "Test case:", str(self.__dict__)

test = test(1000., 1000.)
#t1=time.time()
#dl=[]
#
#for z in sne_z:
#    x=np.log(1./(z+1.0))
#    '''
#    omega=get_omega(x, test.zI, test.omega_mI, test.omega_lI, test.tau)
#    print z, omega[1:,:], test.omega_mI*((1+z)/(1+test.zI))**3
#    hubble = get_hubble(x, test.zI, test.omega_mI, test.omega_lI, test.tau)
#    hubble_lcdm = Hn*test.h*np.sqrt(test.omega_m0*(1.+z)**3 + test.omega_l0)
#    print hubble, hubble_lcdm
#    '''
#    dl.append(get_dl(z, test.zI, test.omega_mI, test.omega_lI, test.tau))
#    dl_lcdm = get_dl_lcdm(z, test.omega_m0, test.omega_l0)
#    print dl, dl_lcdm
#t2=time.time()
#print 'time:', t2-t1
t1=time.time()
dl_interp = get_dl_union2(sne_z, test.zI, test.omega_mI, test.omega_lI,
        test.tau)
t2=time.time()
print 'interpolation time', t2-t1
import matplotlib.pyplot as plt
#plt.figure()
#plt.plot(sne_z, dl, 'g-', sne_z, dl_interp, 'r--')
#plt.savefig('dl.eps')
raw_input()

'''Fill in the chi2 matrix'''
chi2 = np.zeros((n_omega_dm, n_omega_lambda, n_tau))
minchi2 = 100000
for iom, omega_dmI in enumerate(omega_dm_array):
    for iol, omega_lambdaI in enumerate(omega_lambda_array):
        for itau, tau in enumerate(tau_array):
            tmp = get_chi2(num_obs, zI, omega_dmI, omega_lambdaI, tau)
            print iom, iol, itau, 'of', n_omega_dm, n_omega_lambda, n_tau,\
            ':', tmp
            if minchi2 > tmp:
                minchi2 = tmp
                # peakOm, peakOl, peakTau = omega_dmI, omega_lambdaI, tau

            chi2[iom, iol, itau] = tmp
       
print('minchi2 is ', minchi2)
fid="chi2file-"+obs[num_obs]+'-'+str(n_omega_dm)+'-'+str(n_omega_lambda)+'-'+str(n_tau)+'.bin'
chi2file = chi2.tofile(fid)
#np.savetxt('chi2-hr.txt', chi2, fmt="%.2f")
