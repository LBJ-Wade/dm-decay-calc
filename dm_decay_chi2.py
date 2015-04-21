# -*- coding: utf-8 -*-
"""
Constraints on dark matter decay model.

TODO:
    1. use classes to reorganize the code

Created on Wed Jan 07 13:47:42 2015

@author: Hao Wang
@date: 2015-01-07
"""

"""Set parameter scope"""
import numpy as np
import math
from scipy.integrate import odeint
from scipy.integrate import quad
from scan_params import *
import sys

obs = sys.argv[1]

"""OHD data"""
ohd_file = np.loadtxt('ohd.txt', usecols=range(0, 3))
ohd_file = ohd_file[ohd_file[:, 0].argsort()]
ohd_data = ohd_file[:, 1]
ohd_error = ohd_file[:, 2]
ohd_z = ohd_file[:, 0]
ohd_x_output = -np.log(ohd_z+1)

"""SNIa data"""
sne_file = np.loadtxt('SCPUnion2.1_mu_vs_z.txt', usecols=range(1, 4), skiprows=5)
sne_cov = np.loadtxt('SCPUnion2.1_covmat_sys.txt')
sne_inv_cov = np.linalg.inv(sne_cov)

sne_file = sne_file[sne_file[:, 0].argsort()]
sne_data = sne_file[:, 1]
sne_error = sne_file[:, 2]
sne_z = sne_file[:, 0]

"""constants"""
mpc = 3.1e22  # m
kmsmpc = 1000/mpc  # s^-1
H0 = 70*kmsmpc  # s^-1
yr = 3.1e7  # s
t0 = 1.4e10*yr  # s
tau0 = H0*t0  # ~1


"""
ODE of omegam & omegar evolution

n.b. The first element of x_output (or z) must correspond to the initial condition
of y_ini in the ode code -- it's imposed by odeint.

"""
def deriv(y, x, param):
    w = 1./3
    tau = param[0]
    omegal = param[1]
    if y[0]+y[1]+omegal < 0:
        print y, x, tau
        print 'Total density is negative now. [Om, Or, Tau] = ', y[0], y[1], tau
        raise SystemExit

    E = math.sqrt(y[0]+y[1]+omegal)
    du = -y[0]/(tau*E) - 3*y[0]
    dv = +y[0]/(tau*E) - 3*(1+w)*y[1]
    return [du, dv]


"""Solve ODE to get hubble parameter at any x==ln(a); n.b. H0 factor included """
def get_hubble(omegam0, omegar0, tau, x):
    omegal = 1-omegam0-omegar0  # omegal = omegal0
    x_tmp = np.array(x)
    x_output = np.insert(x_tmp, 0, 0.0)
    om_history = odeint(deriv, [omegam0, omegar0], x_output, args=([tau, omegal],))
    omegamHistory = om_history[:, 0]
    omegarHistory = om_history[:, 1]
    return np.sqrt(omegamHistory+omegarHistory+omegal)[1:]*(H0/kmsmpc)


"""Luminosity distance is the integration of 1/H over z"""
def integrand(z, omegam0, omegar0, tau):
    x = np.log(1.0/(1.0+z))
    return 1.0/get_hubble(omegam0, omegar0, tau, x)


"""Calculate D_L (dimensional); used in the calculation of mu (for SNIa)"""
def get_dl(omegam0, omegar0, tau, z):
    c = 3*1e5
    return (1.0+z)*c*quad(integrand, 0, z, args=(omegam0, omegar0, tau))[0]

"""Construct D_L array w.r.t observed redshifts"""
dl_theory = []
for z in sne_z:
    dl_theory.append(round(get_dl(0.2, 0.01, 1000, z), 2))


'''For each parameter set (initial value), get chi2'''   
def get_chi2_ohd(omegam0, omegar0, tau):
    ohd_theory = get_hubble(omegam0, omegar0, tau, ohd_x_output)
    chi2_ohd = np.sum(np.power(ohd_theory-ohd_data, 2)/np.power(ohd_error, 2))
    return chi2_ohd


"""With system error + covariance"""
def get_chi2_sne2(omegam0, omegar0, tau):
    mu = []
    for z in sne_z:
        tmp = 5*np.log10(get_dl(omegam0, omegar0, tau, z))+25.0
        mu.append(tmp)
    chi2_sne=np.dot(mu-sne_data, np.dot(mu-sne_data, sne_inv_cov))
    return chi2_sne


"""No covariance information"""
def get_chi2_sne(omegam0, omegar0, tau):
    chi2_sne=0
    for iz, z in enumerate(sne_z):
        dl = get_dl(omegam0, omegar0, tau, z)
        mu = 5*np.log10(dl)+25.0
        chi2_sne += np.power((mu-sne_data[iz]), 2)/np.power(sne_error[iz], 2)
    return chi2_sne


"""Currently we are not adding up two chi2's, but plot them separately."""
def get_chi2(omegam0, omegar0, tau):
    if obs == 'ohd':
        return get_chi2_ohd(omegam0, omegar0, tau)
    elif obs == 'sne':
        return get_chi2_sne2(omegam0, omegar0, tau)
    else:
        print("No such observation probe available.")
        raise SystemExit


'''Fill in the chi2 matrix'''
chi2 = np.zeros((n_omegam, n_omegar, n_tau))
minchi2 = 100000
for iom, omegam0 in enumerate(omegam_array):
    for iol, omegar0 in enumerate(omegar_array):
        for itau, logtau in enumerate(logtau_array):
            tau = np.exp(logtau)
            tmp = get_chi2(omegam0, omegar0, tau)
            print iom, iol, itau, tmp
            if minchi2 > tmp:
                minchi2 = tmp
                peakOm, peakOr, peakTau = omegam0, omegar0, tau
            chi2[iom, iol, itau] = tmp
       
fid="chi2file-"+obs+'-'+str(n_omegam)+'-'+str(n_omegar)+'-'+str(n_tau)+'.bin'
chi2file = chi2.tofile(fid)
