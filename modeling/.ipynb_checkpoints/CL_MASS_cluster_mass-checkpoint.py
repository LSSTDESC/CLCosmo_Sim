import numpy as np
import pyccl as ccl
import scipy
import clmm
from scipy.integrate import quad,simps, dblquad
from scipy import interpolate
import CL_COUNT_cluster_abundance as cl_count

def binning(edges): return [[edges[i],edges[i+1]] for i in range(len(edges)-1)]


def Cluster_dNd0mega_Mass_ProxyZ(bins, integrand_count = None, grids = None): 
    r"""
    """  
    
    richness_grid, logm_grid, z_grid = grids['richness_grid'], grids['logm_grid'], grids['z_grid']
    z_bins, richness_bins = bins['redshift_bins'], bins['richness_bins']
            
    Nstacked_masses = np.zeros([len(richness_bins), len(z_bins)])

    for i, richness_bin in enumerate(richness_bins):
        #resize richness-axis
        index_richness_grid_cut, richness_grid_cut = cl_count.reshape_axis(richness_grid, richness_bin)
        for j, z_bin in enumerate(z_bins):
            #resize redshift-axis
            index_z_grid_cut, z_grid_cut = cl_count.reshape_axis(z_grid, z_bin)
            m_grid_mat = 10 ** np.tile(logm_grid, (len(z_grid_cut), 1)).T
            
            integrand_cut = cl_count.reshape_integrand(integrand_count, index_richness_grid_cut, index_z_grid_cut)
            
            Nstacked_masses[i,j] = simps(simps(simps(integrand_cut, 
                                         richness_grid_cut, axis=0) * m_grid_mat,
                                         logm_grid, axis=0), 
                                         z_grid_cut, axis=0)
                
    
    return Nstacked_masses
