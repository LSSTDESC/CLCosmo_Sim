import sys
import numpy as np
import pickle
sys.path.append('../../')
import _redshift_richness_bins as analysis

def load(filename, **kwargs):
    """Loads GalaxyCluster object to filename using Pickle"""
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs) 
    
def read_covariance(Richness_bin, Z_bin, photoz = 'TrueZ', radius = None):
    
    if photoz == 'Truez': photoz_load = 'TrueZ'
    elif photoz == 'BPZ': photoz_load = 'BPZ'
    elif photoz == 'flex': photoz_load = 'FlexZ'
    
    path_cov = '../../../CLCosmo_Sim_database/shear-richness_covariance_data/'
    
    cov_DS_selection_bias_obs = np.zeros([len(radius), len(Richness_bin), len(Z_bin)])
    cov_err_DS_selection_bias_obs = np.zeros([len(radius), len(Richness_bin), len(Z_bin)])
    
    for i, z_bin in enumerate(Z_bin):
        data_cov =  load(path_cov + photoz_load + '_lambda20-70'+f'z{z_bin[0]:.2f}-{z_bin[1]:.2f}_cov.pkl')
        for j, richness_bin in enumerate(Richness_bin):
            cov_DS_selection_bias_obs[:,j,i] = data_cov['cov']
            cov_err_DS_selection_bias_obs[:,j,i] = data_cov['cov_err']
    
    return cov_DS_selection_bias_obs, cov_err_DS_selection_bias_obs

def read_gamma(Richness_bin, Z_bin, Gamma1_name = 'HMF_Tinker08_beta1.pkl'):
    
    Gamma1 = np.zeros([len(Richness_bin), len(Z_bin)])
    Gamma1_dictionnary = load(Gamma1_name)
    for i, richness_bin in enumerate(Richness_bin):
        for j, redshift_bin in enumerate(Z_bin):
            name = f'richness={richness_bin[0]}-{richness_bin[1]}_z={redshift_bin[0]:.2f}-{redshift_bin[1]:.2f}'
            Gamma1[i,j] = Gamma1_dictionnary[name]
    
    return Gamma1