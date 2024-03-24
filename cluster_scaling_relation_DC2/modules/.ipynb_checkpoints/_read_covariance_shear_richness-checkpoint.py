import sys
import numpy as np
import pickle
sys.path.append('../../')
import _redshift_richness_bins as analysis

Z_bin = analysis.Z_bin
Richness_bin = analysis.Obs_bin

def load(filename, **kwargs):
    """Loads GalaxyCluster object to filename using Pickle"""
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs) 
    
def read_covariance(photoz = 'True', radius = None):
    
    path_cov='/pbs/throng/lsst/users/cpayerne/cov_CLMassDC2/data/covariance/'
    
    cov_DS_selection_bias_obs = np.zeros([len(radius), len(Richness_bin), len(Z_bin)])
    cov_err_DS_selection_bias_obs = np.zeros([len(radius), len(Richness_bin), len(Z_bin)])
    
    for i, z_bin in enumerate(Z_bin):
        data_cov =  load(path_cov+photoz+f'_z{z_bin[0]:.2f}-{z_bin[1]:.2f}_cov.pkl')
        for j, richness_bin in enumerate(Richness_bin):
            cov_DS_selection_bias_obs[:,j,i] = 1e12*data_cov['cov']
            cov_err_DS_selection_bias_obs[:,j,i] = 1e12*data_cov['cov_err']
    
    return cov_DS_selection_bias_obs, cov_err_DS_selection_bias_obs

def read_gamma():
    
    Gamma1 = np.zeros([len(Richness_bin), len(Z_bin)])
    Gamma1_dictionnary = load('HMF_Tinker08_beta1.pkl')
    for i, richness_bin in enumerate(Richness_bin):
        for j, redshift_bin in enumerate(Z_bin):
            name = f'richness={richness_bin[0]}-{richness_bin[1]}_z={redshift_bin[0]:.2f}-{redshift_bin[1]:.2f}'
            Gamma1[i,j] = Gamma1_dictionnary[name]
    
    return Gamma1