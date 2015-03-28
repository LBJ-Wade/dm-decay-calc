# -*- coding: utf-8 -*-
"""
Constraints on dark matter decay model.

**To do:
    1. break into several parts in modules.
    2. use classes to reorganize the code, then
    3. throw away the current one.

Created on Wed Jan 07 13:47:42 2015

@author: Hao Wang
@date: 2015-01-07
"""

"""Set parameter scope"""
import numpy as np
import math
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import quad
from scan_params import *

"""Set up scan parameters
n_omegam=25
n_omegar=25
n_tau=25
omegam_array=np.linspace(0.09, 0.3, n_omegam)
omegar_array=np.linspace(0.001, 0.09, n_omegar)
logtau_array=np.linspace(0, 30, n_tau)
"""

"""
Read in data.

n.b. The first element of x_output (or z) must correspond to the initial condition
of y_ini in the ode code -- it's imposed by odeint.
"""

"""OHD data"""
ohd_file=np.loadtxt('ohd.txt',usecols=range(0,3))
ohd_file=ohd_file[ohd_file[:,0].argsort()]
ohd_data=ohd_file[:,1]
ohd_error=ohd_file[:,2]
ohd_z=ohd_file[:,0]
ohd_x_output=-np.log(ohd_z+1)

"""SNIa data"""
sne_file=np.loadtxt('SCPUnion2.1_mu_vs_z.txt',usecols=range(1,4),skiprows=5)
sne_cov=np.loadtxt('SCPUnion2.1_covmat_sys.txt')
sne_inv_cov=np.linalg.inv(sne_cov)

sne_file=sne_file[sne_file[:,0].argsort()]
sne_data=sne_file[:,1]
sne_error=sne_file[:,2]
sne_z=sne_file[:,0]

"""constants"""
mpc=3.1e22 #m
kmsmpc=1000/mpc #s^-1
H0=70*kmsmpc #s^-1
yr=3.1e7 #s
t0=1.4e10*yr #s
tau0=H0*t0 #~1

"""ODE of omegam & omegar evolution"""
def deriv(y,x,param):
    w=1./3
    tau=param[0]
    omegal=param[1]
    if y[0]+y[1]+omegal<0:
        print y,x,tau
        print 'Total density is negative now. [Om, Or, Tau] = ', y[0], y[1], tau
        raise SystemExit
#    print 'tau, omegal are ', tau, omegal
#    raw_input('wait:')

    E=math.sqrt(y[0]+y[1]+omegal)
    du = -y[0]/(tau*E) - 3*y[0]
    dv = +y[0]/(tau*E) - 3*(1+w)*y[1]
    return [du,dv]

"""Solve ODE to get hubble parameter at any x==ln(a); n.b. H0 factor included """
def get_hubble(omegam0, omegar0, tau, x):
    omegal=1-omegam0-omegar0 #omegal=omegal0
    x_tmp=np.array(x)
    x_output=np.insert(x_tmp, 0, 0.0)
    #print 'x & x_output', x, x_output
    om_history=odeint(deriv, [omegam0, omegar0], x_output, args=([tau, omegal],)) 
    omegamHistory=om_history[:,0]
    omegarHistory=om_history[:,1]
    #print 'x output again: ', x_output, om_history
    return np.sqrt(omegamHistory+omegarHistory+omegal)[1:]*(H0/kmsmpc)

"""Luminosity distance is the integration of 1/H over z"""
def integrand(z, omegam0, omegar0, tau):
    x=np.log(1.0/(1.0+z))
    return 1.0/get_hubble(omegam0, omegar0, tau, x)

"""Calculate D_L (dimensional); used in the calculation of mu (for SNIa)"""
def get_dl(omegam0, omegar0, tau, z):
    c=3*1e5
    return (1.0+z)*c*quad(integrand, 0, z, args=(omegam0, omegar0, tau))[0]

"""Construct D_L array w.r.t observed redshifts"""
dl_theory=[]
for z in sne_z:
    dl_theory.append(round(get_dl(0.2, 0.01, 1000, z),2))

'''
plt.figure()
plt.plot(sne_z, sne_data)
plt.plot(sne_z, 5*np.log10(dl_theory)+25.0, '+')
#plt.savefig('snePlot.eps')
'''

''' Write in form of ODE
def deriv_dcomov(dcomov, z, param):
    omegam=param[0]
    omegar=param[1]
    tau=param[2]
    x=1./np.log(1.+z)
    
    ddcomov = 1./get_hubble(omegam, omegar, tau, x)
    return [ddcomov]
'''
 
'''For each parameter set (initial value), get chi2'''   
def get_chi2_ohd(omegam0, omegar0, tau):
    ohd_theory=get_hubble(omegam0, omegar0, tau, ohd_x_output)
    chi2_ohd=np.sum(np.power(ohd_theory-ohd_data,2)/np.power(ohd_error,2))       
    return chi2_ohd
    
"""With system error + covariance"""
def get_chi2_sne2(omegam0, omegar0, tau):
    mu=[]
    for z in sne_z:
        tmp=5*np.log10(get_dl(omegam0, omegar0, tau, z))+25.0
        mu.append(tmp)
    chi2_sne=np.dot(mu-sne_data, np.dot(mu-sne_data, sne_inv_cov))
    return chi2_sne
            
"""No covariance information"""
def get_chi2_sne(omegam0, omegar0, tau):
    chi2_sne=0
    for iz, z in enumerate(sne_z):
        dl=get_dl(omegam0, omegar0, tau, z)
        mu=5*np.log10(dl)+25.0
        chi2_sne += np.power((mu-sne_data[iz]),2)/np.power(sne_error[iz], 2)
    return chi2_sne

"""Currently we are not adding up two chi2's, but plot them separately."""
def get_chi2(omegam0, omegar0, tau):
    #return get_chi2_ohd(omegam0, omegar0, tau)
    return get_chi2_sne2(omegam0, omegar0, tau)
        
'''Fill in the chi2 matrix'''
chi2=np.zeros((n_omegam, n_omegar, n_tau))
minchi2=100000
for iom,omegam0 in enumerate(omegam_array):
    for iol, omegar0 in enumerate(omegar_array):
        for itau, logtau in enumerate(logtau_array):
            tau = np.exp(logtau)
            tmp=get_chi2(omegam0, omegar0, tau)
            print iom, iol, itau, tmp
            if minchi2>tmp:
                minchi2=tmp
                peakOm, peakOr, peakTau = omegam0, omegar0, tau
            chi2[iom, iol, itau]=tmp
       
fid="chi2file-"+str(n_omegam)+'-'+str(n_omegar)+'-'+str(n_tau)+'.bin'
chi2file = chi2.tofile(fid)

"""
It's suggested to first obtain the Chi2 result, then plot it out with the
following code.
"""

"""
dchi2=chi2-minchi2
'''Integrate over the third (tau) dimension'''
chi2_omor=-2*np.log(np.sum(np.exp(-dchi2/2),axis=2))
minchi2_omor=np.min(chi2_omor)
[peak_om_index, peak_or_index]=np.unravel_index(np.argmin(chi2_omor),np.shape(chi2_omor))
peak_om=omegam_array[peak_om_index]
peak_or=omegar_array[peak_or_index]
dchi2_omor=chi2_omor-minchi2_omor
contourOmOr=plt.contour(omegar_array, omegam_array, dchi2_omor, levels=[2.3, 5.0])
plt.figure()
plt.plot(peak_or, peak_om, 'r*')
plt.clabel(contourOmOr, inline=1,  extent=(1,3,0,2))
plt.xlabel(r'$\Omega_r$')
plt.ylabel(r'$\Omega_m$')
plt.savefig('or-om-contour.eps')

'''Integrate over the second (omegar) dimension'''
#OM,LOGTAU = np.meshgrid(omegam_array, logtau_array)
chi2_omtau=-2*np.log(np.sum(np.exp(-dchi2/2),axis=1))
minchi2_omtau=np.min(chi2_omtau)
[peak_om_index, peak_logtau_index]=np.unravel_index(np.argmin(chi2_omtau), np.shape(chi2_omtau))
peak_om=omegam_array[peak_om_index]
peak_logtau=logtau_array[peak_logtau_index]
dchi2_omtau=chi2_omtau-minchi2_omtau
contourOmTau=plt.contour(logtau_array, omegam_array, dchi2_omtau, levels=[2.3,5.0])#, levels=[2.3, 5.0]
plt.figure()
plt.plot(peak_logtau, peak_om, 'r*')
plt.clabel(contourOmTau, inline=1)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Omega_m$')
plt.savefig('tau-om-contour.eps')

'''Plot best fitted parameter set
plt.figure()
theory=get_hubble(0.1, 0.06, 35., ohd_x_output)
plt.plot(ohd_z, theory)  
plt.plot(ohd_z, ohd_data, '+')    
'''
"""
