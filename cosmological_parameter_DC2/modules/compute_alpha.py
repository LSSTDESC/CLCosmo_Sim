import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../')
import _redshift_richness_bins

sys.path.append('../../modeling/')
import CL_MASS_fit_alpha

import pickle

def save_pickle(dat, filename, **kwargs):
    print(f"[DEBUG] Saving pickle to: {filename}")
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()
    print("[DEBUG] Pickle saved successfully.\n")

Om_list = np.linspace(0.1, 0.6, 30)
Om_fid = 0.3
Om_range = [0.1, 0.6]

print(f"[DEBUG] Om_list: {Om_list}")
print(f"[DEBUG] Fiducial Om: {Om_fid}")
print(f"[DEBUG] Om_range: {Om_range}\n")

dir_mass = '../../cluster_mass_measurement/cluster_mass_measurement_vary_cosmology/'
dir_alpha = '../alpha/'

names = ['fid_Om{}cluster-masses_1-halo=nfw+c-M=None_rmin=1.0-rmax=3.5_photoz=Truez.pkl']

for name in names:

    print(f"[DEBUG] Processing file name: {name}")

    subpath = name.split('fid_Om{}')[1].split('.pkl')[0]
    where = dir_mass + subpath + '/' + name

    print(f"[DEBUG] Full file path constructed:\n    {where}\n")
    
    print("[DEBUG] Running store_log10M...")
    LogM10_cMNone_, LogM10_cMNone_err_ = CL_MASS_fit_alpha.store_log10M(
        _redshift_richness_bins.Z_bin, 
        _redshift_richness_bins.Obs_bin, 
        Om_list, 
        where
    )
    print("[DEBUG] store_log10M completed.")
    print(f"[DEBUG] LogM10_cMNone_ shape: {np.shape(LogM10_cMNone_)}")
    print(f"[DEBUG] LogM10_cMNone_err_ shape: {np.shape(LogM10_cMNone_err_)}\n")
    
    print("[DEBUG] Running mass–richness model fit...")
    results = CL_MASS_fit_alpha.fit_Mass_Model(
        LogM10_cMNone_, LogM10_cMNone_err_,
        _redshift_richness_bins.Z_bin,
        _redshift_richness_bins.Obs_bin, 
        Om_list, Om_fid,
        Om_range = Om_range,
    )
    print("[DEBUG] fit_Mass_Model completed.\n")
    
    (log10M0, Alpha, log10M0_err, Alpha_err,
     sigma_log10M0, sigma_Alpha, sigma_log10M0_err, sigma_Alpha_err) = results

    params = dict()

    params['Z_bin'] = _redshift_richness_bins.Z_bin
    params['Richness_bin'] = _redshift_richness_bins.Obs_bin
    
    params['log10M0'] = log10M0
    params['log10M0_err'] = log10M0_err
    params['Alpha'] = Alpha
    params['Alpha_err'] = Alpha_err

    params['sigmaM_log10M0'] = sigma_log10M0
    params['sigmaM_log10M0_err'] = sigma_log10M0_err
    params['sigmaM_Alpha'] = sigma_Alpha
    params['sigmaM_Alpha_err'] = sigma_Alpha_err

    name_save = (
        dir_alpha +
        f'alpha_Omfid={Om_fid:.2f}_fit_from_{Om_range[0]:.2f}_to_{Om_range[1]:.2f}_for_' +
        name.split('cluster-masses_')[1]
    )

    print(f"[DEBUG] Output file will be:\n    {name_save}\n")

    # Uncomment when ready to save
    #save_pickle(params, name_save,)
