import matplotlib.pyplot as plt
import pickle
from astropy.coordinates import SkyCoord, match_coordinates_3d, match_coordinates_sky
import sys
import emcee
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.table import Table, QTable, hstack, vstack
from astropy import units as u
import corner
from astropy.coordinates import SkyCoord, match_coordinates_3d
cosmo_astropy = FlatLambdaCDM(H0=71.0, Om0=0.265, Ob0 = 0.0448)
import iminuit
from iminuit import Minuit
cosmo_astropy.critical_density(0.4).to(u.Msun / u.Mpc**3).value

import _utils_cluster_mass_measurement_cosmology_dependent as fit_mass
import _analysis_cluster_mass_measurement as analysis_WL_mean_mass

def save_pickle(dat, filename, **kwargs):
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()

import clmm
import clmm.dataops
from clmm.dataops import compute_tangential_and_cross_components, make_radial_profile, make_bins
from clmm.galaxycluster import GalaxyCluster
import clmm.utils as u
import astropy.units as un
from clmm import Cosmology
from clmm.support import mock_data as mock
import pyccl as ccl

code, analysisname, index_analysis = sys.argv
analysis_WL_metadata = analysis_WL_mean_mass.analysis[str(analysisname)][int(index_analysis)]

#stacked_profiles
#data = np.load(analysis_WL_metadata['data_path'], allow_pickle=True)
#profiles = data['stacked profile']
#covariances = data['stacked covariance']

fix_c = False if analysis_WL_metadata['cM_relation'] == 'None' else True
two_halo_bool = True if analysis_WL_metadata['two_halo']=='True' else False

Om_list = np.linspace(0.1, 0.6, 30)

H0_true = 71
h = H0_true/100
Omega_b_true = 0.02258 / (h**2)
Omega_c_true = 0.1109 / (h**2)
Omega_m_true = Omega_b_true + Omega_c_true
sigma8_true = 0.8
ns_true = 0.963

path = '../../../CLCosmo_Sim_database/data_vary_fiducial_cosmology/'
data = np.load(path+'stacked_esd_profiles_redmapper_vary_Omega_m.pkl', allow_pickle=True)

for k, Om_ in enumerate(Om_list):

    data_profiles = data[f'Om{k}_stacked_profile']
    print(data_profiles.colnames)
    data_covariances = data[f'Om{k}_stacked_covariance']
    
    cosmo_ccl = ccl.Cosmology(Omega_c = Om_  - Omega_b_true, Omega_b = Omega_b_true,
                            h = H0_true/100, sigma8 = sigma8_true, n_s=ns_true, w0=-1, wa=0)

    mass_fit =  fit_mass.fit_WL_cluster_mass(cosmology_ccl = cosmo_ccl, profile = data_profiles, covariance = data_covariances,
                                             colname_radius='radius',
                                                colname_cluster_z='z_mean',
                                                colname_covariance='cov_t',
                                                colname_gt='DSt',                     
                                             a = 0, b =  analysis_WL_metadata['radius_min'], 
                                             rmax = analysis_WL_metadata['radius_max'], 
                                             two_halo_term =two_halo_bool, 
                                             
                                             fix_c = fix_c,
                                             halo_model=analysis_WL_metadata['density_profile'],
                                             mc_relation=analysis_WL_metadata['cM_relation'], 
                                             method='minuit')
    
    res = {'masses':mass_fit, 'analysis':analysis_WL_metadata}
    path = '../../../CLCosmo_Sim/cluster_mass_measurement/cluster_mass_measurement_vary_cosmology/'
    save_pickle(res, path+ 'fid_Om'+str(k) + analysis_WL_metadata['name_save'])
