import pickle
import sys
import numpy as np
sys.path.append('../../lensing_profile_measurement/')
import _utils_stacked_lensing_profiles as prf
sys.path.append('../../')
import _redshift_richness_bins as analysis

def load(filename, **kwargs):
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)
def save_pickle(dat, filename, **kwargs):
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()
path_to_data = '../../../CLCosmo_Sim_database/data_vary_fuducial_cosmology/'
ind_profile = np.load(path_to_data + 'ind_profile_redmapper.pkl', allow_pickle = True)
ind_profile['cluster_id'] = ind_profile['id']
suff_coverage = '_full_coverage'
Z_bin = analysis.Z_bin
Obs_bin = analysis.Obs_bin

if True:
    Om_list = np.linspace(0.1, 0.6, 30)
    stacked_profile_vary_Om = {}
    for k, Om_ in enumerate(Om_list):
        #In bins of richness-redshift bins
        print(k)
        profile_true_stack_Om = prf.stacked_profile(profile = ind_profile,
                            r_in = f'Om{k}_radius',
                            gt_in = f'Om{k}_DSt', gx_in = f'Om{k}_DSx',
                            r_out = f'radius',
                            gt_out ='DSt', gx_out = 'DSx',
                            weight = f'Om{k}_W_l',
                            z_name = 'z', obs_name = 'richness',
                            Z_bin = Z_bin, Obs_bin = Obs_bin, add_columns_to_bin = [f'Om{k}_W_l','richness', 'z',])

        profile_true_stack_Om.rename_column('gt', 'DSt')
        profile_true_stack_Om.rename_column('gx', 'DSx')
        
        
        covariance_true_stack_Om = prf.sample_covariance(profile = ind_profile,
                            r_in = f'Om{k}_radius',
                            gt_in = f'Om{k}_DSt', gx_in = f'Om{k}_DSx',
                            r_out = 'radius',
                            gt_out = 'DSt', gx_out = 'DSx',
                            weight = f'Om{k}_W_l',
                            #n_boot = 600,
                            z_name = 'z', obs_name = 'richness',
                            Z_bin = Z_bin, Obs_bin = Obs_bin)
        
        stacked_profile_vary_Om[f'Om{k}_stacked_profile'] = profile_true_stack_Om
        stacked_profile_vary_Om[f'Om{k}_stacked_covariance'] = covariance_true_stack_Om
    
    stacked_profile_vary_Om[f'Omega_m_value'] = Om_list

    save_pickle(stacked_profile_vary_Om, path_to_data + f'stacked_esd_profiles_redmapper_vary_Omega_m.pkl', allow_pickle=True)

if True:
    w_list = np.linspace(-2.5, -0.5, 30)
    stacked_profile_vary_w = {}
    for k, w_ in enumerate(w_list):
        #In bins of richness-redshift bins
        print(k)
        profile_true_stack_w = prf.stacked_profile(profile = ind_profile,
                            r_in = f'w{k}_radius',
                            gt_in = f'w{k}_DSt', gx_in = f'w{k}_DSx',
                            r_out = f'radius',
                            gt_out ='DSt', gx_out = 'DSx',
                            weight = f'w{k}_W_l',
                            z_name = 'z', obs_name = 'richness',
                            Z_bin = Z_bin, Obs_bin = Obs_bin, add_columns_to_bin = [f'w{k}_W_l','richness', 'z',])

        profile_true_stack_w.rename_column('gt', 'DSt')
        profile_true_stack_w.rename_column('gx', 'DSx')
        
        
        covariance_true_stack_w = prf.sample_covariance(profile = ind_profile,
                            r_in = f'w{k}_radius',
                            gt_in = f'w{k}_DSt', gx_in = f'w{k}_DSx',
                            r_out = 'radius',
                            gt_out = 'DSt', gx_out = 'DSx',
                            weight = f'w{k}_W_l',
                            #n_boot = 600,
                            z_name = 'z', obs_name = 'richness',
                            Z_bin = Z_bin, Obs_bin = Obs_bin)
        
        stacked_profile_vary_w[f'w{k}_stacked_profile'] = profile_true_stack_w
        stacked_profile_vary_w[f'w{k}_stacked_covariance'] = covariance_true_stack_w
    
    stacked_profile_vary_w[f'wDE_value'] = w_list

    save_pickle(stacked_profile_vary_w, path_to_data + f'stacked_esd_profiles_redmapper_vary_wDE.pkl', allow_pickle=True)