import pickle
import sys
import numpy as np
import _utils_stacked_lensing_profiles as prf
sys.path.append('../')
import _redshift_richness_bins as analysis

def load(filename, **kwargs):
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)
def save_pickle(dat, filename, **kwargs):
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()
cluster_id_avoid = np.load('../extract_DC2_data/cluster_id_with_partial_coverage.npy', allow_pickle = True)
ind_profile = np.load('../data/ind_profile_redmapper.pkl', allow_pickle = True)
mask_id_avoid = np.invert(np.isin(ind_profile['id'], cluster_id_avoid))
ind_profile = ind_profile[mask_id_avoid]
ind_profile['cluster_id'] = ind_profile['id']
Z_bin = analysis.Z_bin
Obs_bin = analysis.Obs_bin
suff = '_full_coverage'
profile_true_stack = prf.stacked_profile(profile = ind_profile,
                    r_in = 'radius_true',
                    gt_in = 'DSt_true', gx_in = 'DSx_true',
                    r_out = 'radius',
                    gt_out = 'gt', gx_out = 'gx',
                    weight = 'W_l_true',
                    z_name = 'redshift', obs_name = 'richness',
                    Z_bin = Z_bin, Obs_bin = Obs_bin, add_columns_to_bin = ['W_l_true','richness', 'redshift',])


covariance_true_stack = prf.bootstrap_covariance(profile = ind_profile,
                    r_in = 'radius_true',
                    gt_in = 'DSt_true', gx_in = 'DSx_true',
                    r_out = 'radius',
                    gt_out = 'gt', gx_out = 'gx',
                    weight = 'W_l_true',
                    n_boot = 600,
                    z_name = 'redshift', obs_name = 'richness',
                    Z_bin = Z_bin, Obs_bin = Obs_bin)

t = {'stacked profile':profile_true_stack, 'stacked covariance': covariance_true_stack}
save_pickle(t, f'../data/stacked_esd_profiles_redmapper_true{suff}.pkl', allow_pickle=True)


profile_BPZ_stack = prf.stacked_profile(profile = ind_profile,
 r_in = 'radius_bpz',
                    gt_in = 'DSt_bpz', gx_in = 'DSx_bpz',
                    r_out = 'radius',
                    gt_out = 'gt', gx_out = 'gx',
                    weight = 'W_l_bpz',
                    z_name = 'redshift', obs_name = 'richness',
                    Z_bin = Z_bin, Obs_bin = Obs_bin, add_columns_to_bin = ['W_l_bpz','richness', 'redshift'])
covariance_BPZ_stack = prf.bootstrap_covariance(profile = ind_profile,
                    r_in = 'radius_bpz',
                    gt_in = 'DSt_bpz', gx_in = 'DSx_bpz',
                    r_out = 'radius',
                    gt_out = 'gt', gx_out = 'gx',
                    weight = 'W_l_bpz',
                    n_boot = 400,
                    z_name = 'redshift', obs_name = 'richness',
                    Z_bin = Z_bin, Obs_bin = Obs_bin)
t2 = {'stacked profile':profile_BPZ_stack, 'stacked covariance': covariance_BPZ_stack}
save_pickle(t2, f'../data/stacked_esd_profiles_redmapper_BPZ{suff}.pkl', allow_pickle=True)

profile_flex_stack = prf.stacked_profile(profile = ind_profile,
 r_in = 'radius_flex',
                    gt_in = 'DSt_flex', gx_in = 'DSx_flex',
                    r_out = 'radius',
                    gt_out = 'gt', gx_out = 'gx',
                    weight = 'W_l_flex',
                    z_name = 'redshift', obs_name = 'richness',
                    Z_bin = Z_bin, Obs_bin = Obs_bin, add_columns_to_bin = ['W_l_flex','richness', 'redshift'])

covariance_flex_stack = prf.bootstrap_covariance(profile = ind_profile,
                    r_in = 'radius_flex',
                    gt_in = 'DSt_flex', gx_in = 'DSx_flex',
                    r_out = 'radius',
                    gt_out = 'gt', gx_out = 'gx',
                    weight = 'W_l_flex',
                    n_boot = 400,
                    z_name = 'redshift', obs_name = 'richness',
                    Z_bin = Z_bin, Obs_bin = Obs_bin)
t3 = {'stacked profile':profile_flex_stack, 'stacked covariance': covariance_flex_stack}
save_pickle(t3, f'../data/stacked_esd_profiles_redmapper_flex{suff}.pkl', allow_pickle=True)

