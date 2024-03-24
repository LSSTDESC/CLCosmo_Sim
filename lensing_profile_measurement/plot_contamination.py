import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append('../')
import _redshift_richness_bins as analysis

suff = '_full_coverage'

data_true = np.load(f'../data/stacked_esd_profiles_redmapper_true{suff}.pkl', allow_pickle=True)
data_BPZ = np.load(f'../data/stacked_esd_profiles_redmapper_BPZ{suff}.pkl', allow_pickle=True)
data_flex = np.load(f'../data/stacked_esd_profiles_redmapper_flex{suff}.pkl', allow_pickle=True)

profiles_true = data['stacked profile']
covariances_true = data['stacked covariance']

profiles_true = data['stacked profile']
covariances_true = data['stacked covariance']

profiles_true = data['stacked profile']
covariances_true = data['stacked covariance']
Z_bin = analysis.Z_bin
Obs_bin = analysis.Obs_bin
n_z_bin = len(Z_bin) 
n_m_bin = len(Obs_bin) 
scale = 4
fig, axs = plt.subplots(n_m_bin, n_z_bin, figsize = (10,5))
fig.subplots_adjust(wspace=0, hspace=0)
for i, z_bin in enumerate(Z_bin):
    for j, m_bin in enumerate(Obs_bin):
            label_z =   f'{z_bin[0]:.1f} < z < {z_bin[1]:.1f}'
            label_M = f'{m_bin[0]:.0f} < ' + r'$\lambda$' +f' < {m_bin[1]:.0f}'
            mask_z = (profiles['z_mean'] > z_bin[0])*(profiles['z_mean'] < z_bin[1])
            mask_m = (profiles['obs_mean'] > m_bin[0])*(profiles['obs_mean'] < m_bin[1])
            index = np.arange(len(profiles))
            index_cut = index[mask_m * mask_z]
            f_cut = profiles[index_cut]
            cov = np.array(covariances['cov_t'][index_cut])
            err = cov.T.diagonal()**.5
            axs[j,i].errorbar(f_cut['radius'][0], f_cut['gt'][0] , err[0],
                                c = 'C0',marker = 'o',fmt = ' ', elinewidth = 2, capsize = 5, markersize = 3, markerfacecolor = None)
            axs[j,i].set_ylim(1e12, 3e14)
            axs[j,i].set_xlim(0.4, 12)
            axs[j,i].set_xscale('log')
            axs[j,i].set_yscale('log')
            axs[j,i].tick_params(axis='both', which = 'major', labelsize= 10)
            axs[j,i].legend(frameon = False, loc = 'upper right', fontsize = 10)
            axs[j,i].set_xlabel('R [Mpc]', fontsize = 10)
            axs[j,0].set_ylabel(label_M, fontsize = 10)
            axs[0,i].set_title(label_z, fontsize = 10)
        #except: a=1

for ax in fig.get_axes():
    ax.label_outer()
plt.savefig(f'../fig/stacked_redmapper_profiles_contamination{suff}.png', bbox_inches='tight', dpi=100)
