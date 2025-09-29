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

import _utils_cluster_mass_measurement as fit_mass
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

cosmo = Cosmology(H0 = 71.0, Omega_dm0 = 0.265 - 0.0448, Omega_b0 = 0.0448, Omega_k0 = 0.0)
cosmo_astropy = FlatLambdaCDM(H0=71.0, Om0=0.265, Ob0 = 0.0448)
cosmo_clmm = Cosmology(H0 = 71.0, Omega_dm0 = 0.265 - 0.0448, Omega_b0 = 0.0448, Omega_k0 = 0.0)
cosmo_ccl  = ccl.Cosmology(Omega_c=0.265-0.0448, Omega_b=0.0448, h=0.71, A_s=2.1e-9, n_s=0.96, Neff=0, Omega_g=0)

code, analysisname, index_analysis = sys.argv
analysis_WL_metadata = analysis_WL_mean_mass.analysis[str(analysisname)][int(index_analysis)]

#stacked_profiles
data = np.load(analysis_WL_metadata['data_path'], allow_pickle=True)
profiles = data['stacked profile']
covariances = data['stacked covariance']

fix_c = False if analysis_WL_metadata['cM_relation'] == 'None' else True
two_halo_bool = True if analysis_WL_metadata['two_halo']=='True' else False
mass_fit =  fit_mass.fit_WL_cluster_mass(profile = profiles, covariance = covariances,
                                         a = 0, b =  analysis_WL_metadata['radius_min'], 
                                         rmax = analysis_WL_metadata['radius_max'], 
                                         two_halo_term =two_halo_bool, 
                                         fix_c = fix_c,
                                         halo_model=analysis_WL_metadata['density_profile'],
                                         mc_relation=analysis_WL_metadata['cM_relation'], 
                                         method='minuit')

res = {'masses':mass_fit, 'analysis':analysis_WL_metadata}
save_pickle(res, analysis_WL_metadata['name_save'])
