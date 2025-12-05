import numpy as np
import pickle
from scipy.interpolate import interp1d
def load(filename, **kwargs):
    """Loads GalaxyCluster object to filename using Pickle"""
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)

def load_data_vary_cosmology(analysis_metadata):

    Z_bin, Richness_bin = analysis_metadata['redshift_bins'], analysis_metadata['richness_bins']
    z_corner, rich_corner = analysis_metadata['redshift_corner'], analysis_metadata['richness_corner']

    # Count
    print('[load data]: counts')
    path_to_data = '../../../CLCosmo_Sim_database/data_vary_fuducial_cosmology/'
    table_redmapper = load(path_to_data + 'lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl')
    N_obs, proxy_edges, z_edges = np.histogram2d(table_redmapper['redshift'], 
                                                                table_redmapper['richness'],
                                                           bins=[z_corner, rich_corner])
    if analysis_metadata['type'] == 'N': 
        return N_obs
            
    # Lensing
    if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN':
        print('[load data]: stacked lensing profiles and errors')
        Om_list = np.linspace(0.1, 0.6, 30)
        w_list = np.linspace(-2.5, -0.5, 30)
        path = '../../../CLCosmo_Sim_database/data_vary_fuducial_cosmology/'
        data = np.load(path+'stacked_esd_profiles_redmapper_vary_Omega_m.pkl', allow_pickle=True)

        r = analysis_metadata['radius_centers']
        DS_obs = np.zeros([len(r), len(Richness_bin), len(Z_bin)])
        Err_obs = np.zeros([len(r), len(Richness_bin), len(Z_bin)])

        r_grid = np.zeros(DS_obs.shape)
        for i, z_bin in enumerate(Z_bin):
            for j, richness_bin in enumerate(Richness_bin):
                r_grid[:,j,i] = r

        rup = analysis_metadata['radius_max']
        rlow = analysis_metadata['radius_min']
        mask_is_in_fit_range = (r_grid > rlow)*(r_grid < rup)
        mask_is_in_fit_range_1darray = (r > rlow)*(r < rup)

        DS_interpolator = {}
        DSerr_interpolator = {}

        for i, z_bin in enumerate(Z_bin):
            
            for j, richness_bin in enumerate(Richness_bin):
                
                data_interp_binOm_ds = []
                data_interp_binOm_ds_err = []
                
                for k, Om_ in enumerate(Om_list):
                    
                    ds_datastackOm_k = data['Om'+str(k)+'_stacked_profile']
                    cov_ds_datastackOm_k = data['Om'+str(k)+'_stacked_covariance']
                    maskz = (ds_datastackOm_k['z_mean'] > z_bin[0])*(ds_datastackOm_k['z_mean'] < z_bin[1])
                    maskrichness = (ds_datastackOm_k['obs_mean'] > richness_bin[0])*(ds_datastackOm_k['obs_mean'] < richness_bin[1])
                    maskbin = maskz * maskrichness
                    ds_datastackOm_in_redshift_richness_bin = ds_datastackOm_k[maskbin]
                    cov_ds_datastackOm_in_redshift_richness_bin = cov_ds_datastackOm_k[maskbin]
                    x = np.array(ds_datastackOm_in_redshift_richness_bin['radius'][0])
                    ds = np.array(ds_datastackOm_in_redshift_richness_bin['DSt'][0])
                    ds_err = np.array(cov_ds_datastackOm_in_redshift_richness_bin['cov_t'][0].diagonal()**.5)
                    data_interp_binOm_ds.append(ds)
                    data_interp_binOm_ds_err.append(ds_err)
                
                lin_interp_DS = interp1d(Om_list, data_interp_binOm_ds, axis=0, kind='linear',
                                        bounds_error=False, fill_value='extrapolate')

                lin_interp_DS_err = interp1d(Om_list, data_interp_binOm_ds_err, axis=0, kind='linear',
                                        bounds_error=False, fill_value='extrapolate')

                DSerr_interpolator['z'+str(i) +'_' + 'lambda'+str(j)] = lin_interp_DS_err
                DS_interpolator['z'+str(i) +'_' + 'lambda'+str(j)] = lin_interp_DS

        def DS_obs_fct(Om):

            DS_obs = np.zeros([len(r), len(Richness_bin), len(Z_bin)])
            for i, z_bin in enumerate(Z_bin):
                for j, richness_bin in enumerate(Richness_bin):
                    DS_obs[:,j,i] = DS_interpolator['z'+str(i) +'_' + 'lambda'+str(j)]([Om])
            return DS_obs[mask_is_in_fit_range].reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])

        def Err_obs_fct(Om):

            DSerr_obs = np.zeros([len(r), len(Richness_bin), len(Z_bin)])
            for i, z_bin in enumerate(Z_bin):
                for j, richness_bin in enumerate(Richness_bin):
                    DSerr_obs[:,j,i] = DSerr_interpolator['z'+str(i) +'_' + 'lambda'+str(j)]([Om])
            return DSerr_obs[mask_is_in_fit_range].reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
            

        r_grid = np.zeros(DS_obs.shape)
        for i, z_bin in enumerate(Z_bin):
            for j, richness_bin in enumerate(Richness_bin):
                r_grid[:,j,i] = r
        
        return r, N_obs, DS_obs_fct, Err_obs_fct, mask_is_in_fit_range

    if analysis_metadata['type'] == 'M' or analysis_metadata['type'] == 'MxN':
        print('[load data]: Alpha parameters and log10M0')
        alpha_log10m0 = load(analysis_metadata['cosmo_alpha_log10m0'])
        alpha = alpha_log10m0['alpha']
        log10m0 = alpha_log10m0['log10m0']
        return N_obs, alpha, log10m0, alpha_err, log10m0_err


def load_data(analysis_metadata):

    Z_bin, Richness_bin = analysis_metadata['redshift_bins'], analysis_metadata['richness_bins']
    z_corner, rich_corner = analysis_metadata['redshift_corner'], analysis_metadata['redshift_corner']

    # Count
    print('[load data]: counts')
    path_to_data = '../../../CLCosmo_Sim_database/data/'
    table_redmapper = load(path_to_data + 'lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl')
    N_obs, proxy_edges, z_edges = np.histogram2d(table_redmapper['redshift'], 
                                                                table_redmapper['richness'],
                                                           bins=[z_corner, rich_corner])
    if analysis_metadata['type'] == 'N': 
        return N_obs
    
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
        mask_is_in_fit_range = (r_grid > rlow)*(r_grid < rup)

        DS_obs_mask = DS_obs[mask_is_in_fit_range]
        DS_obs = DS_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])

        Err_obs_mask = Err_obs[mask_is_in_fit_range]
        Err_obs = Err_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
        
        return N_obs, DS_obs, Err_obs, mask_is_in_fit_range
