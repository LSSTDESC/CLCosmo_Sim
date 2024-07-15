import numpy as np
import pickle
def load(filename, **kwargs):
    """Loads GalaxyCluster object to filename using Pickle"""
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)

def load_data(analysis_metadata, Z_bin, Richness_bin, z_corner, rich_corner):

    # Count
    print('[load data]: counts')
    table_redmapper = load('../../data/lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl')
    N_obs, proxy_edges, z_edges = np.histogram2d(table_redmapper['redshift'], 
                                                                table_redmapper['richness'],
                                                           bins=[z_corner, rich_corner])
    if analysis_metadata['type'] == 'N': return N_obs
    
    #Masses
    if analysis_metadata['type'] == 'M' or analysis_metadata['type'] == 'MxN':
        print('[load data]: stacked masses and errors')
        mass_data = load(analysis_metadata['mass_file'])['masses']
        log10Mass = np.zeros(N_obs.T.shape)
        log10Mass_err = np.zeros(N_obs.T.shape)
        for i, richness_bin in enumerate(Richness_bin):
            for j, redshift_bin in enumerate(Z_bin):
                maskz = (mass_data['z_mean'] > redshift_bin[0])*(mass_data['z_mean'] < redshift_bin[1])
                maskr = (mass_data['obs_mean'] > richness_bin[0])*(mass_data['obs_mean'] < richness_bin[1])
                mask = maskz * maskr
                mass_in_bin = mass_data[mask]
                log10Mass[i,j] = np.array(mass_in_bin['log10M200c_WL'])
                log10Mass_err[i,j] = np.array(mass_in_bin['err_log10M200c_WL'])
        return N_obs, log10Mass, log10Mass_err
            
    # Lensing
    if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN':
        print('[load data]: stacked lensing profiles and errors')
        data = np.load(analysis_metadata['lensing_data'], allow_pickle=True)
        profiles = data['stacked profile']
        covariances = data['stacked covariance']
        r = profiles['radius'][0]
        DS_obs = np.zeros([len(r), len(Richness_bin), len(Z_bin)])
        Err_obs = np.zeros([len(r), len(Richness_bin), len(Z_bin)])

        for i, z_bin in enumerate(Z_bin):
            mask_z = (profiles['z_mean'] > z_bin[0])*(profiles['z_mean'] < z_bin[1])
            for j, richness_bin in enumerate(Richness_bin):
                mask_richness = (profiles['obs_mean'] > richness_bin[0])*(profiles['obs_mean'] < richness_bin[1])
                mask_tot = mask_z * mask_richness
                DS_obs[:,j,i] = profiles['gt'][mask_tot][0]
                Err_obs[:,j,i] = covariances['cov_t'][mask_tot][0].diagonal()**.5

        r_grid = np.zeros(DS_obs.shape)
        for i, z_bin in enumerate(Z_bin):
            for j, richness_bin in enumerate(Richness_bin):
                r_grid[:,j,i] = r

        rup = analysis_metadata['radius_max']
        rlow = analysis_metadata['radius_min']
        mask = (r_grid > rlow)*(r_grid < rup)

        DS_obs_mask = DS_obs[mask]
        DS_obs = DS_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])

        Err_obs_mask = Err_obs[mask]
        Err_obs = Err_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
        
        return N_obs, DS_obs, Err_obs