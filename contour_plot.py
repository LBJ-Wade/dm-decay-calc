from io_params import *
import sys

'''
The only argument is choosing OHD('ohd') or SNe('sne')
'''
if len(sys.argv) != 2 or sys.argv[1] == '-h':
    print("[*]Usage: python dm_decay_chi2.py obs")
    print("[*]obs is a sequence formed of the following numbers:")
    print(obs)
    print("[*]Separate your choices with a comma without spaces.")
    print("[*]Duplicate terms will be counted only once.")
    raise SystemExit

#chi2 = np.zeros((n_omega_dm, n_omega_lambda, n_tau))
obs_list = sorted(set(sys.argv[1].split(',')))
chi2list = []

chi2_total = np.zeros((n_omega_dm, n_omega_lambda, n_tau))
name_plot = ''

print('Using {')
num_obs_list = []
for inum, num_obs in enumerate(obs_list):
    num_obs = int(num_obs)
    num_obs_list.append(num_obs)
    print obs_dict[num_obs]
    name_plot += obs_dict[num_obs]
    if inum != len(obs_list)-1:
        name_plot += '-'
    fid = 'chi2file-'+obs_dict[num_obs]+'-'+str(n_omega_dm)+'-'+\
            str(n_omega_lambda)+'-'+str(n_tau)+'.bin'
    chi2list.append(np.fromfile(fid).reshape(n_omega_dm, n_omega_lambda,
        n_tau))
    chi2_total += chi2list[inum]
    
chi2list.append(chi2_total)
num_obs_list.append(10)
print('} as observables.')

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

fig = Figure()
canvas = FigureCanvas(fig)
ax  = fig.add_subplot(121)
ax.set_xlabel(r'$\Omega_\Lambda$', fontsize=big_font)
ax.set_ylabel(r'$\Omega_{dm}(z=1000)$', fontsize=big_font)
bx = fig.add_subplot(122)
bx.set_xlabel(r'$\tau$', fontsize=big_font)
bx.set_xscale("log")

for ichi2, chi2 in enumerate(chi2list):
    minchi2 = np.amin(chi2)
    dchi2 = chi2-minchi2

    '''Integrate over the third (tau) dimension'''
    chi2_omol = -2*np.log(np.sum(np.exp(-dchi2/2.), axis=2))
    minchi2_omol = np.min(chi2_omol)
    dchi2_omol = chi2_omol-minchi2_omol
    color = color_dict.get(num_obs_list[ichi2])
    ax.contour(omega_lambda_array, omega_dm_array, 
            dchi2_omol, levels=[2.3, 5.0], colors=(color,), linestyles=['solid',
            'dashed'])
    # Mark the best-fit point of the comibined chi2.
    if ichi2 == len(chi2list)-1:
        [peak_om_index, peak_ol_index] = np.unravel_index(np.argmin(chi2_omol), 
            np.shape(chi2_omol))
        peak_om = omega_dm_array[peak_om_index]
        peak_ol = omega_lambda_array[peak_ol_index]
        ax.plot(peak_ol, peak_om, 'r*')

    '''Integrate over the second (omega_lambda) dimension'''
    chi2_omtau = -2*np.log(np.sum(np.exp(-dchi2/2.), axis=1))
    minchi2_omtau = np.min(chi2_omtau)
    dchi2_omtau = chi2_omtau-minchi2_omtau
    bx.contour(tau_array, omega_dm_array, dchi2_omtau, 
            levels=[2.3, 5.0], colors=color, linestyles=['solid', 'dashed'])
    if ichi2 == len(chi2list)-1:
        [peak_om_index, peak_tau_index] = np.unravel_index(np.argmin(chi2_omtau), 
                np.shape(chi2_omtau))
        peak_om = omega_dm_array[peak_om_index]
        peak_tau = tau_array[peak_tau_index]
        bx.plot(peak_tau, peak_om, 'r*')



canvas.print_figure(name_plot+'.eps')
