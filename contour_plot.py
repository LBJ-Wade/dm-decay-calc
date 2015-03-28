import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from scan_params import *

"""
n_omegam = 25
n_omegar = 25
n_tau	 = 25
omegam_array=np.linspace(0.09, 0.3, n_omegam)
omegar_array=np.linspace(0.001, 0.09, n_omegar)
logtau_array=np.linspace(0, 30, n_tau)
"""

fid="chi2file-"+str(n_omegam)+'-'+str(n_omegar)+'-'+str(n_tau)+'.bin'

chi = np.fromfile(fid)
chi2= chi.reshape(n_omegam, n_omegar, n_tau)
minchi2=np.amin(chi2)
dchi2=chi2-minchi2

'''Integrate over the third (tau) dimension'''
chi2_omor=-2*np.log(np.sum(np.exp(-dchi2/2),axis=2))
#print chi2_omor
minchi2_omor=np.min(chi2_omor)
[peak_om_index, peak_or_index]=np.unravel_index(np.argmin(chi2_omor),np.shape(chi2_omor))
peak_om=omegam_array[peak_om_index]
peak_or=omegar_array[peak_or_index]
dchi2_omor=chi2_omor-minchi2_omor
plt.figure()
contourOmOr=plt.contour(omegar_array, omegam_array, dchi2_omor, levels=[2.3, 5.0])
plt.plot(peak_or, peak_om, 'r*')
plt.clabel(contourOmOr, inline=1,  extent=(1,3,0,2))
plt.xlabel(r'$\Omega_r$')
plt.ylabel(r'$\Omega_m$')
plt.savefig('or-om.eps')

'''Integrate over the second (omegar) dimension'''
chi2_omtau=-2*np.log(np.sum(np.exp(-dchi2/2),axis=1))
minchi2_omtau=np.min(chi2_omtau)
[peak_om_index, peak_logtau_index]=np.unravel_index(np.argmin(chi2_omtau), np.shape(chi2_omtau))
peak_om=omegam_array[peak_om_index]
peak_logtau=logtau_array[peak_logtau_index]
dchi2_omtau=chi2_omtau-minchi2_omtau
plt.figure()
contourOmTau=plt.contour(logtau_array, omegam_array, dchi2_omtau, levels=[2.3,5.0])
plt.plot(peak_logtau, peak_om, 'r*')
plt.clabel(contourOmTau, inline=1)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Omega_m$')
plt.savefig('tau-om.eps')
