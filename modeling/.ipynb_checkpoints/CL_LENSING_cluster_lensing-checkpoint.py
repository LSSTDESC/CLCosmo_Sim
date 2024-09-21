import numpy as np
import pyccl as ccl
import scipy
import clmm
import clmm.utils as u
from scipy.integrate import quad,simps, dblquad
from scipy import interpolate
import CL_COUNT_cluster_abundance as cl_count
import pyccl as ccl
import cluster_toolkit as ct
from tqdm.auto import tqdm, trange
import multiprocessing
import time
global fct_n

def binning(edges): return [[edges[i],edges[i+1]] for i in range(len(edges)-1)]

def compute_cluster_lensing(R, cM, logm_grid, z_grid, cosmo_ccl, cosmo_clmm, 
                            two_halo = False, 
                            boost_factor=False,
                            return_two_halo_term=False, miscentering=False, 
                            fmis=0.3, Rmis=0.2, multiprocessing=False):
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
    return_two_halo_term: Bool
        return the two halo term separate to the one-halo term
    boost_factor: Bool
        use the boost factor or not
     miscentering: Bool
         use miscentering or not
    fmis: float
        fraction of miscentered object
    Rmis: float
        Rmis of the Gamma distribution
    multiprocessing: Bool
        use multiprocessing or not to compute the miscentered stacked signal
         
    Returns:
    --------
    cluster_lensing_grid: array
        tabulated lensing signal
    """
    #concentration
    concentration_grid = np.zeros([len(logm_grid), len(z_grid)])
    deff = ccl.halos.massdef.MassDef(200, 'critical')
    if cM == 'Diemer15': conc = ccl.halos.concentration.ConcentrationDiemer15(mass_def=deff)
    if cM == 'Duffy08': conc = ccl.halos.concentration.ConcentrationDuffy08(mass_def=deff)
    if cM == 'Prada12': conc = ccl.halos.concentration.ConcentrationPrada12(mass_def=deff)
    if cM == 'Bhattacharya13': conc = ccl.halos.concentration.ConcentrationBhattacharya13(mass_def=deff)
    for i, z in enumerate(z_grid):
        concentration_grid[:,i] = conc._concentration(cosmo_ccl, 10**logm_grid, 1./(1. + z))
    
    #excess surface density
    if boost_factor:
        boost_factor_grid = np.zeros([ len(R), len(logm_grid), len(z_grid)])
        rs_grid = np.zeros([len(logm_grid), len(z_grid)])
        for i, z in enumerate(z_grid):
            r200c = clmm.compute_rdelta(mdelta=10**logm_grid, redshift=z, 
                                        cosmo=cosmo_clmm, massdef='critical', delta_mdef=200)
            rs_grid[:,i] = r200c / concentration_grid[:,i]
            for j, logm in enumerate(logm_grid):
                boost_factor_grid[:,j,i] = u.compute_nfw_boost(R, rscale=rs_grid[j,i], boost0=.1)
    
    cluster_lensing_grid = np.zeros([ len(R), len(logm_grid), len(z_grid)])
    
    if not miscentering:
        for i in range(len(logm_grid)):
            for j in range(len(z_grid)):
                cluster_lensing_grid[:,i,j] = excess_surface_density_fct(R, logm_grid[i], concentration_grid[i,j], z_grid[j], cosmo_clmm)
            
    elif miscentering:
        if multiprocessing:
            cluster_lensing_grid_mis = np.zeros([ len(R), len(logm_grid), len(z_grid)])
            Logm_grid = np.tile(logm_grid, (len(z_grid),1)).T
            Z_grid = np.tile(z_grid, (len(logm_grid),1))
            Logm_grid_flatten, Z_grid_flatten, Concentration_grid_flatten = Logm_grid.flatten(), Z_grid.flatten(), concentration_grid.flatten()
            global excess_surface_density_mis_fct_n
            def excess_surface_density_mis_fct_n(n):
                return excess_surface_density_mis_fct(R, Logm_grid_flatten[n], Concentration_grid_flatten[n], Z_grid_flatten[n], cosmo_clmm)[0]
            n_points = len(Logm_grid_flatten)
            res_mp=map(excess_surface_density_mis_fct_n, np.arange(n_points), ncores=8, ordered=True,)
            excess_surface_density_mis = np.array(res_mp).T.reshape((len(R), len(logm_grid), len(z_grid)))
            for i in range(len(logm_grid)):
                for j in range(len(z_grid)):
                    excess_surface_density = excess_surface_density_fct(R, logm_grid[i], concentration_grid[i,j], z_grid[j], cosmo_clmm)
                    cluster_lensing_grid[:,i,j] = excess_surface_density_mis[:,i,j] * fmis + (1 - fmis) * excess_surface_density

        else:
            for i in range(len(logm_grid)):
                for j in range(len(z_grid)):
                    excess_surface_density = excess_surface_density_fct(R, logm_grid[i], concentration_grid[i,j], z_grid[j], cosmo_clmm)
                    ti = time.time()
                    excess_surface_density_mis, surface_density_mis = excess_surface_density_mis_fct(R, logm_grid[i], 
                                                                                                     concentration_grid[i,j], z_grid[j], 
                                                                                                     cosmo_clmm, Rmis=Rmis)
                    tf = time.time()
                    cluster_lensing_grid[:,i,j] = excess_surface_density_mis * fmis + (1 - fmis) * excess_surface_density
    if two_halo:
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
    
    if boost_factor: cluster_lensing_grid=cluster_lensing_grid/boost_factor_grid                                          
        
    if return_two_halo_term: return cluster_lensing_grid_2h, cluster_lensing_grid
    else: return cluster_lensing_grid

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


def miscentering_distribution(R_off, sigma_off, which='gamma'):
    r"""
    Attributes:
    ----------
    R_off : array, float
        Miscentering radius in Mpc
    sigma_off : float
        miscentering scattering in Mpc
    Returns:
    -------
    P_off : the probability density function of the miscentering radius 
    at R_off for the rayleigh distribution
    """
    if which=='raleigh': return ( R_off/(sigma_off**2) ) * np.exp(-0.5*(R_off/sigma_off)**2)
    elif which=='gamma': return ( R_off/(sigma_off**2) ) * np.exp(-R_off/sigma_off)

def excess_surface_density_mis_fct(R, logm, concentration, z, cosmo_clmm, Rmis=0.11):
    "predict the off-centered stacked excess surface density profile"
    R_Sigma = np.logspace(np.log10(0.0001), np.log10(32), 500)
    Sigma_mis = clmm.compute_surface_density(R_Sigma, 10**logm, 
                                             concentration, z, 
                                             cosmo_clmm, delta_mdef=200, 
                                             halo_profile_model='nfw', 
                                             massdef='critical')
    excess_surface_density_mis, surface_density_mis = predict_sigma_excess_miscentering_single(Sigma_mis, R, R_Sigma, Rmis, z, "gamma", cosmo_clmm)
    return excess_surface_density_mis, surface_density_mis

def excess_surface_density_fct(R, logm, concentration, z, cosmo_clmm):
    "predict the standard well-centered excess surface density profile"
    excess_surface_density = clmm.compute_excess_surface_density(R, 10**logm, 
                                                                concentration, z, 
                                                                cosmo_clmm, delta_mdef=200,
                                                                halo_profile_model='nfw',
                                                                massdef='critical')
    return excess_surface_density


def predict_sigma_excess_miscentering_single(Sigma, R, R_Sigma, Rmis, cluster_z, kernel, cosmo):
    r"""
    Attributes:
    ----------
    Sigma : array
        the surface mass density computed for R_Sigma Msun/Mpc^2
    R : array, float
        Radius for DeltaSigma with miscentering to be computed in Mpc
    R_Sigma: array
        radial axis for Sigma in Mpc
    Rmis : float
        the miscentering scattering in Mpc
    cluster_z : float
        the cluster redshift
    kernel : string
        'rayleigh' for the rayleigh distribution
        'gamma' for the gamma distribution
    moo : modelling object of clmm wiht 'cosmo' attribute
    Returns:
    -------
    The miscentered excess surface density relative to Sigma in Msun/Mpc^2
    """
    
    Omega_m, h = cosmo.get_Omega_m(cluster_z), cosmo['h']
    Sigma_ct_unit = Sigma/((10**12)*h*(1 + cluster_z)**2) #comoving hMsun/Mpc^2
    Rcomoving = h*(1 + cluster_z)*R #comoving Mpc/h
    Rsigma = h*(1 + cluster_z)*R_Sigma #comoving Mpc/h
    Rmis = h*(1 + cluster_z)*Rmis #comoving Mpc/h
    cluster_m_mis, concentration_mis = 3e14, 5 #Msun, no_units
    Sigma_mis = ct.miscentering.Sigma_mis_at_R(Rsigma, Rsigma, Sigma_ct_unit, cluster_m_mis*h, concentration_mis, Omega_m, Rmis, kernel= kernel) #comoving hMsun/Mpc^2
    DeltaSigma_mis = ct.miscentering.DeltaSigma_mis_at_R(Rcomoving, Rsigma, Sigma_mis) #comoving hMsun/Mpc^2
    return DeltaSigma_mis*((10**12)*h*(1 + cluster_z)**2), Sigma_mis*((10**12)*h*(1 + cluster_z)**2) #comoving hMsun/Mpc^2


def _map_f(args):
    f, i = args
    return f(i)

def map(func, iter, ncores=3, ordered=True):
    
    ncpu = multiprocessing.cpu_count()
    print('You have {0:1d} CPUs'.format(ncpu))
    pool = multiprocessing.Pool(processes=ncpu) 
    inputs = ((func,i) for i in iter)
    res_list = []
    if ordered: pool_map = pool.imap
    else: pool_map = pool.imap_unordered
    with tqdm(total=len(iter), desc='# progress ...') as pbar:
        for res in pool_map(_map_f, inputs):
            try :
                pbar.update()
                res_list.append(res)
            except KeyboardInterrupt:
                pool.terminate()
    pool.close()
    pool.join()
    return res_list