import matplotlib
matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab.
import matplotlib.pyplot as plt
from scan_params import *
import sys

if len(sys.argv) != 2:
    print("[*]Usage: python contour_plot.py obs")
    print("[*]obs='ohd' or 'sne'")

obs = sys.argv[1]
big = 18

fid = "chi2file-"+obs+'-'+str(n_omega_dm)+'-'+str(n_omega_lambda)+'-'+str(n_tau)+'.bin'
#fid = "chi2file-"+str(n_omega_dm)+'-'+str(n_omega_lambda)+'-'+str(n_tau)+'-sum.bin'

chi = np.fromfile(fid)
chi2 = chi.reshape(n_omega_dm, n_omega_lambda, n_tau)
minchi2 = np.amin(chi2)
dchi2 = chi2-minchi2

'''Integrate over the third (tau) dimension'''
chi2_omol = -2*np.log(np.sum(np.exp(-dchi2/2), axis=2))
minchi2_omol = np.min(chi2_omol)
[peak_om_index, peak_ol_index] = np.unravel_index(np.argmin(chi2_omol), np.shape(chi2_omol))
peak_om = omega_dm_array[peak_om_index]
peak_ol = omega_lambda_array[peak_ol_index]
dchi2_omol = chi2_omol-minchi2_omol
plt.figure()
contourOmOl = plt.contour(omega_lambda_array, omega_dm_array, dchi2_omol, levels=[2.3, 5.0])
plt.plot(peak_ol, peak_om, 'r*')
# plt.clabel(contourOmOl, inline=1,  extent=(1, 3, 0, 2))
plt.xlabel(r'$\Omega_r$', fontsize=big)
plt.ylabel(r'$\Omega_m$', fontsize=big)
plt.savefig('or-om-'+obs+'.eps')

'''Integrate over the second (omega_lambda) dimension'''
chi2_omtau = -2*np.log(np.sum(np.exp(-dchi2/2), axis=1))
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
plt.savefig('tau-om-'+obs+'.eps')
