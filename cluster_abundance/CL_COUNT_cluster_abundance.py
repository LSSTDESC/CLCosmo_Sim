import numpy as np
import pyccl as ccl
import scipy
from scipy.integrate import quad,simps, dblquad
from scipy import interpolate

import CL_COUNT_modeling_completeness as comp
import CL_COUNT_modeling_purity as pur
import CL_COUNT_modeling_halo_mass_function as hmf
import CL_COUNT_modeling_richness_mass_relation as rm_relation

def binning(edges): return [[edges[i],edges[i+1]] for i in range(len(edges)-1)]

def reshape_axis(axis, bounds):
    
    index = np.arange(len(axis))
    mask = (axis > bounds[0])*(axis < bounds[1])
    index_cut = index[mask]
    axis_cut = axis[mask]
    axis_cut[0], axis_cut[-1] = bounds
    return index_cut, axis_cut

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
        dNdzdlogMdOmega_grid[:,i] = hmf.dndlog10M(logm_grid ,z, cosmo, hmd) * hmf.dVdzdOmega(z, cosmo)
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
    purity_flat = pur.purity(R_flat, Z_flat, theta_purity)
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
    completeness_flat = comp.completeness(Logm_flat, Z_flat, theta_completeness)
    completeness_grid = np.reshape(completeness_flat, [len(logm_grid), len(z_grid)])
    return completeness_grid

def compute_richness_mass_relation_grid(richness_grid, logm_grid, z_grid, theta_rm):
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
    richness_tab = np.zeros([len(richness_grid), len(logm_grid), len(z_grid)])
    logm_tab = np.zeros([len(richness_grid), len(logm_grid), len(z_grid)])
    z_tab = np.zeros([len(richness_grid), len(logm_grid), len(z_grid)])
    for i, richness in enumerate(richness_grid): richness_tab[i,:,:] = richness
    for i, logm in enumerate(logm_grid): logm_tab[:,i,:] = logm
    for i, z in enumerate(z_grid): z_tab[:,:,i] = z
    mu = rm_relation.proxy_mu_f(logm_tab, z_tab, theta_rm)
    sigma = rm_relation.proxy_sigma_f(logm_tab, z_tab, theta_rm)
    pdf = (1/richness_tab)*np.exp(-.5*(np.log(richness_tab)-mu)**2/(sigma**2))/(np.sqrt(2*np.pi*sigma**2))
    #pdf = rm_relation.pdf_richness_mass_relation(richness_tab, logm_tab, z_tab, theta_rm)
    return pdf

def recompute_count_modelling(count_modelling, compute = None, grids = None, params = None):
    r"""
    recompute the tabulated maps for cluster abundance integrand
    Attributes:
    -----------
    count_modelling: dict
        dictionary of tabulated maps of hmf, richness-mass relation, purity, completeness
    compute: dict
        dictionnary of boolean to choose to compute
    grids: dict
        dictionary of tabulated axis for the mass, richness, redshift
    params: dict
        dictionary of params (cosmology, M-lambda relation, purity and completeness
    Returns:
    --------
    count_modelling: dict
        dictionary of tabulated (recomputed) maps of hmf, richness-mass relation, purity, completeness
    """
    
    richness_grid, logm_grid, z_grid = grids['richness_grid'], grids['logm_grid'], grids['z_grid']
    shape_integrand = [len(richness_grid), len(logm_grid), len(z_grid)]
    
    if compute['compute_purity']:
        purity_ = compute_purity_grid(grids['richness_grid'], grids['z_grid'], params['params_purity'])
        purity_new = np.zeros(shape_integrand)
        for i in range(shape_integrand[1]):
            purity_new[:,i,:] = purity_
        count_modelling['purity'] = purity_new
            
    if compute['compute_completeness']:
        completeness_ = compute_completeness_grid(grids['logm_grid'], grids['z_grid'], params['params_completeness'])
        completeness_new = np.zeros(shape_integrand)
        for i in range(shape_integrand[0]):
            completeness_new[i,:,:] = completeness_
        count_modelling['completeness'] = completeness_new
     
    if compute['compute_richness_mass_relation']:
        count_modelling['richness_mass_relation'] = compute_richness_mass_relation_grid(grids['richness_grid'], 
                                                                         grids['logm_grid'], 
                                                                         grids['z_grid'], 
                                                                         params['params_richness_mass_relation'])

    if compute['compute_dNdzdlogMdOmega']:
        dNdzdlogMdOmega_ = compute_dNdzdlogMdOmega_grid(grids['logm_grid'], grids['z_grid'], 
                                                        params['CCL_cosmology'], params['halo_mass_distribution'])
        dNdzdlogMdOmega_new = np.zeros(shape_integrand)
        for i in range(shape_integrand[0]):
            dNdzdlogMdOmega_new[i,:,:] = dNdzdlogMdOmega_
        count_modelling['dNdzdlogMdOmega'] = dNdzdlogMdOmega_new
    
    return count_modelling

def define_count_integrand(count_modelling, adds):
    r"""
    define count integrand with the option of considering purity and/or completeness
    Attributes:
    ----------
    count_modelling: dict
        dictionnary of tabulated hmf, purity, completeness, richness-mass relation
    adds: dict
        dictionary of booleans, choose if purity and/or completeness are considered
    Returns:
    --------
    integrand: array
        integrand on the mass, richness and redshift axis
    """
    dNdzdlogMdOmega = count_modelling['dNdzdlogMdOmega']
    richness_mass_relation = count_modelling['richness_mass_relation']
    purity = count_modelling['purity']
    completeness = count_modelling['completeness']
    
    if adds['add_purity']:
        if adds['add_completeness']:
            integrand = dNdzdlogMdOmega * richness_mass_relation * completeness / purity
        else:
            integrand = dNdzdlogMdOmega * richness_mass_relation / purity
    else:
        if adds['add_completeness']:
            integrand = dNdzdlogMdOmega * richness_mass_relation * completeness
        else:
            integrand = dNdzdlogMdOmega * richness_mass_relation
    return integrand

def reshape_integrand(integrand_count, index_richness_grid_cut, index_z_grid_cut):
    r"""reshape integrand_count on selected indexes for the richness axis and redshift axis
    Attributes:
    -----------
    integrand_count: array
        integrand over the mass, richness and redshift axis
    index_richness_grid_cut: array
        indexes on the richness axis
    index_richness_grid_cut: array
        indexes on the redshift axis
    Returns:
    --------
    integrand_cut: array
        reshaped integrand
    """
    integrand_cut = integrand_count[index_richness_grid_cut,:,:]
    integrand_cut = integrand_cut[:,:,index_z_grid_cut]
    return integrand_cut
    
def Cluster_SurfaceDensity_ProxyZ(bins, integrand_count = None, grids = None): 
 
    r"""
    Attributes:
    -----------
    bins: dict
        dictionnary of redshift and richness bins
    integrand_count: array
        tabulated integrand, to be integrated over the amss, richness and redshift axis
    grids: dict
        dictionnary of tabulated arrays of the mass, richness, redshift axis
    Returns:
    --------
    dNdOmega: array
        angular surface density of clusters in bins of redshift and mass
    """
    richness_grid, logm_grid, z_grid = grids['richness_grid'], grids['logm_grid'], grids['z_grid']
    z_bins, richness_bins = bins['redshift_bins'], bins['richness_bins']
    
    dNdOmega = np.zeros([len(richness_bins), len(z_bins)])
    for i, richness_bin in enumerate(richness_bins):
        #resize richness-axis
        index_richness_grid_cut, richness_grid_cut = reshape_axis(richness_grid, richness_bin)
        for j, z_bin in enumerate(z_bins):
            #resize redshift-axis
            index_z_grid_cut, z_grid_cut = reshape_axis(z_grid, z_bin)
            integrand_cut = reshape_integrand(integrand_count, index_richness_grid_cut, index_z_grid_cut)
            integral = simps(simps(simps(integrand_cut, 
                                         richness_grid_cut, axis=0), 
                                         logm_grid, axis=0), 
                                         z_grid_cut, axis=0)
            dNdOmega[i,j] = integral
    return dNdOmega

# def Cluster_Abundance_MZ(z_bins, logm_bins, dNdzdlogMdOmega, completeness,
#                              logm_grid, z_grid, add_completeness): 
#     r"""
#     Attributes:
#     -----------
#     z_bins : array
#         list of redshift bins
#     richness_bins : array
#         list of richness bins
#     dNdzdlogMdOmega : array
#         dNdzdlogMdOmega grid
#     purity : array
#         purity grid
#     completeness : array
#         completeness grid
#     richness_mass_relation : array
#         richness_mass_relation grid
#     richness_grid : array
#         richness grid
#     logm_grid : array
#         logm grid
#     z_grid : array
#         redshift grid
#     add_purity : Bool
#         add purity ?
#     add_completeness : Bool
#         add completeness ?
#     Returns:
#     --------
#     dNdOmega : ndarray
#         Cluster abundance prediction in redshift and proxy bins
#     """  
#     dNdOmega = np.zeros([len(logm_bins), len(z_bins)])

#     if add_completeness:
#         integrand = dNdzdlogMdOmega *  completeness
#     else:
#         integrand = dNdzdlogMdOmega
            
#     for i, logm_bin in enumerate(logm_bins):
#         #resize richness-axis
#         index_logm_grid_cut, logm_grid_cut = reshape_axis(logm_grid, logm_bin)
#         for j, z_bin in enumerate(z_bins):
#             #resize redshift-axis
#             index_z_grid_cut, z_grid_cut = reshape_axis(z_grid, z_bin)
#             integrand_cut = integrand[index_logm_grid_cut,:]
#             integrand_cut = integrand_cut[:,index_z_grid_cut]
#             integral = simps(simps(integrand_cut, logm_grid_cut, axis=0),  z_grid_cut, axis=0)
#             dNdOmega[i,j] = integral
#     return dNdOmega, integrand

#def abundance_mass(z_bins, logm_bins, 
#                   cosmo, hmd, 
#                   theta_completeness,
#                   logm_grid, z_grid,
#                   dNdzdlogMdOmega=None, completeness=None, 
#                   add_completeness=True,
#                   compute_dNdzdlogMdOmega=True, 
#                   compute_completeness=True):
#
#    shape_integrand = [len(logm_grid), len(z_grid)]
#    if compute_completeness:
#        completeness_ = compute_completeness_grid(logm_grid, z_grid, theta_completeness)
#        completeness = completeness_
#
#    if compute_dNdzdlogMdOmega:
#        dNdzdlogMdOmega_ = compute_dNdzdlogMdOmega_grid(logm_grid, z_grid, cosmo, hmd)
#        dNdzdlogMdOmega = dNdzdlogMdOmega_
#        
#    dNdOmega, integrand = Cluster_Abundance_MZ(z_bins, logm_bins,
#                             dNdzdlogMdOmega, completeness,
#                          logm_grid, z_grid, add_completeness)
#    
#    return dNdOmega, integrand
    
    
    
