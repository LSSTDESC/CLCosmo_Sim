import numpy as np
import pyccl as ccl
import scipy
import clmm
from scipy.integrate import quad,simps, dblquad
from scipy import interpolate
import CL_COUNT_cluster_abundance as cl_count
import pyccl as ccl

def binning(edges): return [[edges[i],edges[i+1]] for i in range(len(edges)-1)]

def compute_cluster_lensing(R, cM, logm_grid, z_grid, cosmo_ccl, cosmo_clmm, two_halo = False):
    r"""
    Compute the lensing signals on a logm; z grid
    Attributes:
    -----------
    R: array
        radius for the lensing profiles
    cM: str
        name of a mass-richness relation
    logm_grid: array
        grid for mass axis
    z_grid: array
        grid for redshift axis
    cosmo_ccl: object
        CCL cosmology object
    cosmo_clmm: object
        CLMM cosmology object
    two_halo: Bool
        use two halo term or not
    Returns:
    --------
    cluster_lensing_grid: array
        tabulated lensing signal
    """
    
    #concentration
    concentration_grid = np.zeros([len(logm_grid), len(z_grid)])
    deff = ccl.halos.massdef.MassDef(200, 'critical')
    if cM == 'Diemer15':
        conc = ccl.halos.concentration.ConcentrationDiemer15(mass_def=deff)
    if cM == 'Duffy08':
        conc = ccl.halos.concentration.ConcentrationDuffy08(mass_def=deff)
    if cM == 'Prada12':
        conc = ccl.halos.concentration.ConcentrationPrada12(mass_def=deff)
    if cM == 'Bhattacharya13':
        conc = ccl.halos.concentration.ConcentrationBhattacharya13(mass_def=deff)
    for i, z in enumerate(z_grid):
        lnc = np.log(conc._concentration(cosmo_ccl, 10**logm_grid, 1./(1. + z))) 
        concentration_grid[:,i] = np.exp(lnc)  
    #excess surface density
    
    cluster_lensing_grid = np.zeros([ len(R), len(logm_grid), len(z_grid)])
    
    for i in range(len(logm_grid)):
        for j in range(len(z_grid)):
            cluster_lensing_grid[:,i,j] = clmm.compute_excess_surface_density(R, 10**logm_grid[i], 
                                                                            concentration_grid[i,j], z_grid[j], 
                                                                            cosmo_clmm, delta_mdef=200,
                                                                            halo_profile_model='nfw',
                                                                            massdef='critical')
    if two_halo == True:
        
        definition_200m = ccl.halos.massdef.MassDef(200, 'matter')
        halobias_200m_fct = ccl.halos.hbias.tinker10.HaloBiasTinker10(mass_def='200m', mass_def_strict=True)
        m200c_to_m200m = ccl.halos.massdef.mass_translator(mass_in=deff, mass_out=definition_200m, concentration=conc)
        logm_grid_200m = np.zeros([len(logm_grid), len(z_grid)])
        halo_bias_200m = np.zeros([len(logm_grid), len(z_grid)])

        for i, z in enumerate(z_grid):
            logm_grid_200m[:,i] = np.log10(m200c_to_m200m(cosmo_ccl, 10**logm_grid, 1/(1+z)))
            halo_bias_200m[:,i] = halobias_200m_fct(cosmo_ccl, 10**logm_grid_200m[:,i], 1/(1+z))

        cluster_lensing_grid_2h = np.zeros([ len(R), len(logm_grid), len(z_grid)])
        for j in range(len(z_grid)):
            ds_2h = clmm.compute_excess_surface_density_2h(R, z_grid[j], cosmo=cosmo_clmm)
            for i in range(len(logm_grid)):
                cluster_lensing_grid_2h[:,i,j] = ds_2h * halo_bias_200m[i,j]
        cluster_lensing_grid = cluster_lensing_grid + cluster_lensing_grid_2h
        
    return cluster_lensing_grid

def Cluster_dNd0mega_Lensing_ProxyZ(bins, integrand_count = None, cluster_lensing = None, lensing_radius = None, grids = None): 
    r"""
    """  
    
    richness_grid, logm_grid, z_grid = grids['richness_grid'], grids['logm_grid'], grids['z_grid']
    z_bins, richness_bins = bins['redshift_bins'], bins['richness_bins']
            
    stacked_lensing_profile = np.zeros([len(lensing_radius), len(richness_bins), len(z_bins)])
    for i, richness_bin in enumerate(richness_bins):
        #resize richness-axis
        index_richness_grid_cut, richness_grid_cut = cl_count.reshape_axis(richness_grid, richness_bin)
        for j, z_bin in enumerate(z_bins):
            #resize redshift-axis
            index_z_grid_cut, z_grid_cut = cl_count.reshape_axis(z_grid, z_bin)
            integrand_cut = cl_count.reshape_integrand(integrand_count, index_richness_grid_cut, index_z_grid_cut)
            #pdf
            cluster_lensing_cut = cluster_lensing[:,:,index_z_grid_cut]
            for k, r in enumerate(lensing_radius):
                stacked_lensing_profile[k,i,j] = simps(simps(simps(integrand_cut, 
                                         richness_grid_cut, axis=0) * cluster_lensing_cut[k,:,:], 
                                         logm_grid, axis=0), 
                                         z_grid_cut, axis=0)
                
    
    return stacked_lensing_profile
