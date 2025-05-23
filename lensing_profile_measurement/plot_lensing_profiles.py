import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append('../')
import _redshift_richness_bins as analysis

plt.figure(figsize=(7,4))
index = 12
fmt = ['-', '--', '.']
suff = '_full_coverage'
path_to_data = '../../CLCosmo_Sim_database/data/'
data_true = np.load(path_to_data + f'stacked_esd_profiles_redmapper_true{suff}.pkl', allow_pickle=True)
profiles_true = data_true['stacked profile']
z = data_true['stacked profile']['z_mean']
for z_bin in analysis.Z_bin:
    mask = (z > z_bin[0])*(z < z_bin[1])
    r = np.mean(profiles_true['radius'][mask], axis=0)
    ds = np.mean(profiles_true['gt'][mask], axis=0)
    plt.loglog(r, ds, label=f'{z_bin[0]:.1f} < z < {z_bin[1]:.1f}')
plt.xscale('log')
plt.ylabel(r'$\Delta\Sigma$', fontsize=14)
plt.xlabel('R', fontsize=14)
plt.legend()
plt.tick_params(axis='both', which="both", labelsize= 13)
plt.savefig('../fig/stacked_redmapper_profiles_attenuation.png', bbox_inches='tight', dpi=100)


data = np.load(path_to_data + 'ind_profile_redmapper.pkl', allow_pickle=True)
plt.figure(figsize=(7,4))

mask_z = (data['redshift'] < .5)*(data['redshift'] > .2)
mask = mask_z*(data['richness'] < 40)
r = np.average(data['radius_true'][mask],weights = data['W_l_true'][mask],  axis=0)
ds = np.average(data['DSt_true'][mask],weights = data['W_l_true'][mask],  axis=0)
plt.errorbar(r, ds, np.std(data['DSt_true'][mask], axis=0)/np.sqrt(len(data['DSt_true'][mask])), fmt = '-',color = 'C0', marker = 'o', label = '20 < lambda < 40', markerfacecolor = 'w')

mask = mask_z*(data['richness'] > 40)*(data['richness'] < 100)
r = np.average(data['radius_true'][mask], weights = data['W_l_true'][mask], axis=0)
ds = np.average(data['DSt_true'][mask], weights = data['W_l_true'][mask], axis=0)
plt.errorbar(r, ds, np.std(data['DSt_true'][mask], axis=0)/np.sqrt(len(data['DSt_true'][mask])), fmt = '-', color = 'C1',marker = 'o', label = '40 < lambda < 100', markerfacecolor = 'w')

mask = mask_z*(data['richness'] > 100)
r = np.average(data['radius_true'][mask],weights = data['W_l_true'][mask],  axis=0)
ds = np.average(data['DSt_true'][mask],weights = data['W_l_true'][mask],  axis=0)
plt.errorbar(r, ds, np.std(data['DSt_true'][mask], axis=0)/np.sqrt(len(data['DSt_true'][mask])), fmt = '-', color = 'C2',marker = 'o', label = '100 < lambda', markerfacecolor = 'w')

mask_z = (data['redshift'] < 1)*(data['redshift'] > .5)
mask = mask_z*(data['richness'] < 40)
r = np.average(data['radius_true'][mask],weights = data['W_l_true'][mask],  axis=0)
ds = np.average(data['DSt_true'][mask],weights = data['W_l_true'][mask],  axis=0)
plt.errorbar(r, ds, np.std(data['DSt_true'][mask], axis=0)/np.sqrt(len(data['DSt_true'][mask])), fmt = '--', color = 'C0', marker = 'o')

mask = mask_z*(data['richness'] > 40)*(data['richness'] < 100)
r = np.average(data['radius_true'][mask], weights = data['W_l_true'][mask], axis=0)
ds = np.average(data['DSt_true'][mask], weights = data['W_l_true'][mask], axis=0)
plt.errorbar(r, ds, np.std(data['DSt_true'][mask], axis=0)/np.sqrt(len(data['DSt_true'][mask])), fmt = '--', color = 'C1',marker = 'o')

mask = mask_z*(data['richness'] > 100)
r = np.average(data['radius_true'][mask],weights = data['W_l_true'][mask],  axis=0)
ds = np.average(data['DSt_true'][mask],weights = data['W_l_true'][mask],  axis=0)
plt.errorbar(r, ds, np.std(data['DSt_true'][mask], axis=0)/np.sqrt(len(data['DSt_true'][mask])), fmt = '--', color = 'C2',marker = 'o')

plt.plot([], [], '--k', label = '0.5 < z < 1')
plt.plot([], [], '-k', label = '0.2 < z < 0.5')

plt.yscale('log')

#plt.xlim(0.2, 12)
plt.ylabel(r'$\Delta\Sigma$', fontsize=14)
plt.xlabel('R', fontsize=14)
plt.legend()
plt.tick_params(axis='both', which="both", labelsize= 13)
#plt.savefig('../fig/lensing_profiles.png', bbox_inches='tight', dpi=300)

data = np.load(path_to_data + f'stacked_esd_profiles_redmapper_true{suff}.pkl', allow_pickle=True)
profiles = data['stacked profile']
covariances = data['stacked covariance']
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
plt.savefig(f'../fig/stacked_redmapper_profiles_true_source_redshifts{suff}.png', bbox_inches='tight', dpi=100)

plt.figure(figsize=(7,4))
index = 12
fmt = ['-', '--', '.']
label = ['cosmoDC2 - BPZ redshifts', 'cosmoDC2 - FlexZBoost redshifts']    
data_true = np.load(path_to_data + f'stacked_esd_profiles_redmapper_true{suff}.pkl', allow_pickle=True)
profiles_true = data['stacked profile']
z = data_true['stacked profile']['z_mean'][index]
for i, data_name in enumerate([path_to_data + f'stacked_esd_profiles_redmapper_BPZ{suff}.pkl', path_to_data + f'stacked_esd_profiles_redmapper_flex{suff}.pkl']):
    covariances_true = data['stacked covariance']
    data = np.load(data_name, allow_pickle=True)
    profiles = data['stacked profile']
    covariances = data['stacked covariance']
    plt.plot(profiles['radius'][index], (profiles['gt'][index] - profiles_true['gt'][index])/profiles_true['gt'][index],
               fmt[i], label = label[i])
plt.xscale('log')
plt.ylabel(r'$\Delta\Sigma$ (biais wrt to true $z_{\rm source}$, %)', fontsize=14)
plt.xlabel('R', fontsize=14)
plt.legend()
plt.title(f'<z> = {z:.2f}')
plt.tick_params(axis='both', which="both", labelsize= 13)
plt.savefig(f'../fig/stacked_redmapper_profiles_true_BPZ_flex_redshifts{suff}.png', bbox_inches='tight', dpi=100)
