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
path_to_data = '../../../CLCosmo_Sim_database/data_vary_fiducial_cosmology/ind_gammat_profile_redmapper.pkl'
ind_profile = np.load(path_to_data, allow_pickle = True)
ind_profile['cluster_id'] = ind_profile['id']
suff_coverage = '_full_coverage'
Z_bin = analysis.Z_bin
Obs_bin = analysis.Obs_bin
stacked_profile_vary = {}

profile_true_stack = prf.stacked_profile(profile = ind_profile,
                    r_in = f'radius',
                    gt_in = f'DSt', gx_in = f'DSx',
                    r_out = f'radius',
                    gt_out ='DSt', gx_out = 'DSx',
                    weight = f'W_l',
                    z_name = 'z', obs_name = 'richness',
                    Z_bin = Z_bin, Obs_bin = Obs_bin, add_columns_to_bin = [f'W_l','richness', 'z',])

profile_true_stack.rename_column('gt', 'DSt')
profile_true_stack.rename_column('gx', 'DSx')

covariance_true_stack = prf.sample_covariance(profile = ind_profile,
                    r_in = f'radius',
                    gt_in = f'DSt', gx_in = f'DSx',
                    r_out = 'radius',
                    gt_out = 'DSt', gx_out = 'DSx',
                    weight = f'W_l',
                    #n_boot = 600,
                    z_name = 'z', obs_name = 'richness',
                    Z_bin = Z_bin, Obs_bin = Obs_bin)

stacked_profile_vary[f'stacked_profile'] = profile_true_stack
stacked_profile_vary[f'stacked_covariance'] = covariance_true_stack

save_pickle(stacked_profile_vary, path_to_data + f'stacked_gammat_profiles_redmapper.pkl', allow_pickle=True)