import sys
sys.path.append('../../../modeling/')
import CL_COUNT_modeling_richness_mass_relation as scalingrel
import numpy as np

def mean_scaling_relation(logm_array, z, chains, log10m0 = 14, z0 = 0.5):
    
    mean_rel = []
    for i in range(len(chains)):
        theta_rm = [log10m0, z0] + list(chains[i])
        mean_rel.append(scalingrel.proxy_mu_f(logm_array, z, theta_rm))
    return np.mean(mean_rel, axis=0), np.std(mean_rel, axis=0)
