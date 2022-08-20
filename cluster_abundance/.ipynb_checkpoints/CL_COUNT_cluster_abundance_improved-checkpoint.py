import numpy as np
import pyccl as ccl
import scipy
from scipy.integrate import quad,simps, dblquad
from scipy import interpolate

def binning(edges): return [[edges[i],edges[i+1]] for i in range(len(edges)-1)]

def reshape_axis(axis, bounds):
    
    index = np.arange(len(axis))
    mask = (axis > bounds[0])*(axis < bounds[1])
    index_cut = index[mask]
    axis_cut = axis[mask]
    axis_cut[0], axis_cut[-1] = bounds
    return index_cut, axis_cut

def dndlog10M(logm, z, cosmo, hmd):
    r"""
    Attributes:
    -----------
    log10M : array
        \log_{10}(M), M dark matter halo mass
    z : float
        halo redshift
    cosmo: CCL cosmology object
        cosmological parameters
    hmd: CCL hmd object
        halo definition
    Returns:
    --------
    hmf : array
        halo mass function for the corresponding masses and redshift
    """
    hmf = hmd.get_mass_function(cosmo, 10**np.array(logm), 1./(1. + z))
    return hmf

def dVdzdOmega(z, cosmo):
    r"""
    Attributes:
    ----------
    z : float
        redshift
    cosmo: CCL cosmology
        cosmological parameters
    Returns:
    -------
    dVdzdOmega_value : float
        differential comoving volume 
    """
    a = 1./(1. + z)
    da = ccl.background.angular_diameter_distance(cosmo, a)
    ez = ccl.background.h_over_h0(cosmo, a) 
    dh = ccl.physical_constants.CLIGHT_HMPC / cosmo['h']
    dVdzdOmega_value = dh * da * da/( ez * a ** 2)
    return dVdzdOmega_value

def completeness(logm, z, theta_completeness):
    r"""
    Attributes:
    -----------
    log10M : array
        \log_{10}(M), M dark matter halo mass
    z : float
        halo redshift
    theta_completeness: array
        parameters of completeness
    Returns:
    --------
    completeness : array
        completeness of cluster detection
    """
    logm0, c0, c1, nc = theta_completeness
    logm_scale = logm0 + c0 + c1*(1 + z)
    m_rescale = 10**logm/(10**logm_scale)
    return m_rescale**nc/(m_rescale**nc+1)

def purity(richness, z, theta_purity):
    r"""
    Attributes:
    -----------
    richness : array
        cluster richness
    z : float
        cluster redshift
    theta_purity: array
        parameters of purity
    Returns:
    --------
    purity : array
        purity of cluster detection
    """
    richness0, p0, p1, np_ = theta_purity
    richness_scale = np.log(richness0) + p0 + p1*(1 + z)
    richness_rescale = richness/richness_scale
    return richness_rescale**np_/(richness_rescale**np_+1)

def richness_mass_relation(richness, logm, z, theta_rm):
    r"""
    Attributes:
    -----------
    richness : array
        cluster richness
    logm: array
        logm of halo mass
    z : float
        cluster redshift
    theta_rm: array
        parameters of purity
    Returns:
    --------
    rm : array
        richness-mass relation
    """
    log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = theta_rm
    proxy_mu = proxy_mu0 + proxy_muz * np.log((1+z)/(1 + z0)) + proxy_mulog10m * (logm-log10m0)
    proxy_sigma = proxy_sigma0 + proxy_sigmaz * np.log((1+z)/(1 + z0)) + proxy_sigmalog10m * (logm-log10m0)
    return (1/richness)*np.exp(-.5*(np.log(richness)-proxy_mu)**2/proxy_sigma**2)/np.sqrt(2*np.pi*proxy_sigma**2)

def compute_dNdzdlogMdOmega_grid(logm_grid, z_grid, cosmo, hmd):
    r"""
    Attributes:
    -----------
    logm_grid : array
        log10M tabulated axis
    z_grid : array
        redshift tabulated axis
    cosmo : CCL cosmology
        cosmology
    hmd; CCL hmd object
        halo definition
    Returns:
    --------
    dN_dzdlogMdOmega : array
        tabulated dzdlogMdOmega
    """
    dNdzdlogMdOmega_grid = np.zeros([len(logm_grid), len(z_grid)])
    for i, z in enumerate(z_grid):
        dNdzdlogMdOmega_grid[:,i] = dndlog10M(logm_grid ,z, cosmo, hmd) * dVdzdOmega(z, cosmo)
    return dNdzdlogMdOmega_grid
    
def compute_purity_grid(richness_grid, z_grid, theta_purity):
    r"""
    based on https://arxiv.org/pdf/1611.05468.pdf
    Attributes:
    -----------
    richness_grid : array
        richness tabulated axis
    z_grid : array
        redshift tabulated axis
    theta_purity : array
        parameters of purity
    Returns:
    --------
    purity_grid : array
        tabulated purity
    """
    R, Z = np.meshgrid(richness_grid, z_grid)
    R_flat, Z_flat = R.T.flatten(), Z.T.flatten()
    purity_flat = purity(R_flat, Z_flat, theta_purity)
    purity_grid = np.reshape(purity_flat, [len(richness_grid), len(z_grid)])
    return purity_grid

def compute_completeness_grid(logm_grid, z_grid, theta_completeness):
    r"""
    based on https://arxiv.org/pdf/1611.05468.pdf
    Attributes:
    -----------
    logm_grid : array
        log10M tabulated axis
    z_grid : array
        redshift tabulated axis
    theta_completeness : array
        parameters of purity
    Returns:
    --------
    completeness_grid : array
        tabulated completeness
    """
    Logm, Z = np.meshgrid(logm_grid, z_grid)
    Logm_flat, Z_flat = Logm.T.flatten(), Z.T.flatten()
    completeness_flat = completeness(Logm_flat, Z_flat, theta_completeness)
    completeness_grid = np.reshape(completeness_flat, [len(logm_grid), len(z_grid)])
    return completeness_grid

def compute_richness_mass_relation_grid(richness_grid, logm_grid, z_grid,  theta_rm):
    r"""
    based on https://arxiv.org/pdf/1904.07524.pdf
    Attributes:
    -----------
    richness_grid : array
        richness tabulated axis
    logm_grid : array
        log10M tabulated axis
    z_grid : array
        redshift tabulated axis
    theta_rm : array
        parameters of purity
    Returns:
    --------
    rm_grid : array
        tabulated richness-mass relation
    """
    rm_relation_grid = np.zeros([len(richness_grid), len(logm_grid), len(z_grid)])
    for i in range(len(logm_grid)):
        for j in range(len(z_grid)):
            pdf = richness_mass_relation(richness_grid, logm_grid[i], z_grid[j], theta_rm)
            rm_relation_grid[:,i,j] = pdf
            #add normalization
    return rm_relation_grid

def Cluster_Abundance_ProxyZ(z_bins, richness_bins, 
                             dNdzdlogMdOmega, purity, completeness, richness_mass_relation,
                             richness_grid, logm_grid, z_grid, 
                             add_purity, add_completeness): 
    r"""
    Attributes:
    -----------
    z_bins : array
        list of redshift bins
    richness_bins : array
        list of richness bins
    dNdzdlogMdOmega : array
        dNdzdlogMdOmega grid
    purity : array
        purity grid
    completeness : array
        completeness grid
    richness_mass_relation : array
        richness_mass_relation grid
    richness_grid : array
        richness grid
    logm_grid : array
        logm grid
    z_grid : array
        redshift grid
    add_purity : Bool
        add purity ?
    add_completeness : Bool
        add completeness ?
    Returns:
    --------
    dNdOmega : ndarray
        Cluster abundance prediction in redshift and proxy bins
    """  
    dNdOmega = np.zeros([len(richness_bins), len(z_bins)])
    if add_purity:
        if add_completeness:
            integrand = dNdzdlogMdOmega * richness_mass_relation * completeness / purity
        else:
            integrand = dNdzdlogMdOmega * richness_mass_relation / purity
    else:
        if add_completeness:
            integrand = dNdzdlogMdOmega * richness_mass_relation * completeness
        else:
            integrand = dNdzdlogMdOmega * richness_mass_relation
            
    for i, richness_bin in enumerate(richness_bins):
        #resize richness-axis
        index_richness_grid_cut, richness_grid_cut = reshape_axis(richness_grid, richness_bin)
        for j, z_bin in enumerate(z_bins):
            #resize redshift-axis
            index_z_grid_cut, z_grid_cut = reshape_axis(z_grid, z_bin)
            integrand_cut = integrand[index_richness_grid_cut,:,:]
            integrand_cut = integrand_cut[:,:,index_z_grid_cut]
            integral = simps(simps(simps(integrand_cut, 
                                         richness_grid_cut, axis=0), 
                                         logm_grid, axis=0), 
                                         z_grid_cut, axis=0)
            dNdOmega[i,j] = integral
    return dNdOmega, integrand

def abundance_proxy(z_bins, richness_bins, cosmo, hmd, theta_purity, theta_completeness, theta_rm,
             richness_grid, logm_grid, z_grid,
             dNdzdlogMdOmega=None, purity=None, completeness=None, richness_mass_relation=None, 
             add_purity=True, add_completeness=True,
             compute_dNdzdlogMdOmega=True, compute_purity=True, compute_completeness=True, compute_richness_mass_relation=True):
    r"""
    Attributes:
    -----------
    z_bins : array
        list of redshift bins
    richness_bins : array
        list of richness bins
    cosmo : CCL Cosmology
        cosmological aprameters
    hmd : CCL hmd object
        halo definition
    theta_purity : array
        purity parameters
    theta_completeness : array
        completeness parameters
    theta_rm : array
        richness-mass relation parameters
    richness_grid : array
        richness grid
    logm_grid : array
        logm grid
    z_grid : array
        redshift grid
    dNdzdlogMdOmega : array
        dNdzdlogMdOmega grid
    purity : array
        purity grid
    completeness : array
        completeness grid
    richness_mass_relation : array
        richness_mass_relation grid
    compute_dNdzdlogMdOmega : array
        compute dNdzdlogMdOmega grid ?
    compute_purity : array
        compute_purity grid ?
    compute_completeness : array
        compute_completeness grid ?
    compute_richness_mass_relation : array
        compute_richness_mass_relation grid ?
    Returns:
    --------
    dNdOmega : ndarray
        Cluster abundance prediction in redshift and proxy bins
    integrand : ndarray
        3d integrand
    """  
    shape_integrand = [len(richness_grid), len(logm_grid), len(z_grid)]
            
    if compute_purity:
        purity_ = compute_purity_grid(richness_grid, z_grid, theta_purity)
        purity = np.zeros(shape_integrand)
        for i in range(shape_integrand[1]):
            purity[:,i,:] = purity_
            
    if compute_completeness:
        completeness_ = compute_completeness_grid(logm_grid, z_grid, theta_completeness)
        completeness = np.zeros(shape_integrand)
        for i in range(shape_integrand[0]):
            completeness[i,:,:] = completeness_
     
    if compute_richness_mass_relation:
        richness_mass_relation = compute_richness_mass_relation_grid(richness_grid, logm_grid, z_grid,  theta_rm)
        
    if compute_dNdzdlogMdOmega:
        dNdzdlogMdOmega_ = compute_dNdzdlogMdOmega_grid(logm_grid, z_grid, cosmo, hmd)
        dNdzdlogMdOmega = np.zeros(shape_integrand)
        for i in range(shape_integrand[0]):
            dNdzdlogMdOmega[i,:,:] = dNdzdlogMdOmega_
        
    dNdOmega, integrand = Cluster_Abundance_ProxyZ(z_bins, richness_bins,
                             dNdzdlogMdOmega, purity, completeness, richness_mass_relation,
                             richness_grid, logm_grid, z_grid, add_purity, add_completeness)
    
    return dNdOmega, integrand

def Cluster_Abundance_MZ(z_bins, logm_bins, dNdzdlogMdOmega, completeness,
                             logm_grid, z_grid, add_completeness): 
    r"""
    Attributes:
    -----------
    z_bins : array
        list of redshift bins
    richness_bins : array
        list of richness bins
    dNdzdlogMdOmega : array
        dNdzdlogMdOmega grid
    purity : array
        purity grid
    completeness : array
        completeness grid
    richness_mass_relation : array
        richness_mass_relation grid
    richness_grid : array
        richness grid
    logm_grid : array
        logm grid
    z_grid : array
        redshift grid
    add_purity : Bool
        add purity ?
    add_completeness : Bool
        add completeness ?
    Returns:
    --------
    dNdOmega : ndarray
        Cluster abundance prediction in redshift and proxy bins
    """  
    dNdOmega = np.zeros([len(logm_bins), len(z_bins)])

    if add_completeness:
        integrand = dNdzdlogMdOmega *  completeness
    else:
        integrand = dNdzdlogMdOmega
            
    for i, logm_bin in enumerate(logm_bins):
        #resize richness-axis
        index_logm_grid_cut, logm_grid_cut = reshape_axis(logm_grid, logm_bin)
        for j, z_bin in enumerate(z_bins):
            #resize redshift-axis
            index_z_grid_cut, z_grid_cut = reshape_axis(z_grid, z_bin)
            integrand_cut = integrand[index_logm_grid_cut,:]
            integrand_cut = integrand_cut[:,index_z_grid_cut]
            integral = simps(simps(integrand_cut, logm_grid_cut, axis=0),  z_grid_cut, axis=0)
            dNdOmega[i,j] = integral
    return dNdOmega, integrand

def abundance_mass(z_bins, logm_bins, 
                   cosmo, hmd, 
                   theta_completeness,
                   logm_grid, z_grid,
                   dNdzdlogMdOmega=None, completeness=None, 
                   add_completeness=True,
                   compute_dNdzdlogMdOmega=True, 
                   compute_completeness=True):

    shape_integrand = [len(logm_grid), len(z_grid)]
    if compute_completeness:
        completeness_ = compute_completeness_grid(logm_grid, z_grid, theta_completeness)
        completeness = completeness_

    if compute_dNdzdlogMdOmega:
        dNdzdlogMdOmega_ = compute_dNdzdlogMdOmega_grid(logm_grid, z_grid, cosmo, hmd)
        dNdzdlogMdOmega = dNdzdlogMdOmega_
        
    dNdOmega, integrand = Cluster_Abundance_MZ(z_bins, logm_bins,
                             dNdzdlogMdOmega, completeness,
                          logm_grid, z_grid, add_completeness)
    
    return dNdOmega, integrand
    
    
    
