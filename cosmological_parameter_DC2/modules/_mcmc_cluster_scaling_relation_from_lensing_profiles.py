import sys
import pyccl as ccl
import numpy as np
from clmm import Cosmology
from multiprocessing import Pool
import emcee
import matplotlib.pyplot as plt
import time
import pickle
import logging
logger = logging.getLogger('Fit blinded data')
logging.basicConfig(
     format="{asctime} - {levelname} - {message}",
     style="{",
     datefmt="%Y-%m-%d %H:%M:%S",
     level=logging.INFO,)

import _load_data
import _analysis_scaling_relation as analysis_list
def save_pickle(dat, filename, **kwargs):
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()
def load(filename, **kwargs):
    """Loads GalaxyCluster object to filename using Pickle"""
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)

import _read_covariance_shear_richness as GammaLambda_Cov

sys.path.append('../../')
import _redshift_richness_bins as analysis
sys.path.append('../../lensing_profile_measurement/')
import _config_lensing_profiles

sys.path.append('../../modeling')
import CL_COUNT_modeling_completeness as comp
import CL_COUNT_DATAOPS_cluster_abundance_covariance as cl_covar 
import CL_COUNT_modeling_purity as pur
import CL_COUNT_modeling_halo_mass_function as hmf
import CL_COUNT_modeling_richness_mass_relation as rm_relation
import CL_MASS_cluster_mass as cl_mass
import CL_COUNT_cluster_abundance as cl_count
import CL_COUNT_class_likelihood as likelihood
import CL_LENSING_cluster_lensing as cl_lensing
CLCount_likelihood = likelihood.Likelihood()
import argparse

def mcmc(analysis_metadata):

    #this is where to put argparser
    print()
    for k in analysis_metadata.keys():
        logger.info(f'analysis_metadata[{k}] = {analysis_metadata[k]}')
    fit_cosmo = analysis_metadata['fit_cosmo']

    #binning in redshift
    if analysis_metadata['redshift_range'] == 'Full': 
        Z_bin = analysis.Z_bin 
        z_corner = analysis.z_corner
    elif analysis_metadata['redshift_range'] == 'Partial':
        Z_bin = analysis.Z_bin[analysis_metadata['redshift_bin_index']]
        z_corner = analysis.z_corner[analysis_metadata['redshift_corner_index']]
        
    if analysis_metadata['richness_range'] == 'Full': 
        Richness_bin = analysis.Obs_bin 
        rich_corner = analysis.rich_corner
    elif analysis_metadata['richness_range'] == 'Partial':
        Richness_bin = analysis.Obs_bin[analysis_metadata['richness_bin_index']]
        rich_corner = analysis.rich_corner[analysis_metadata['richness_corner_index']]
        
    bins = {'redshift_bins':Z_bin, 'richness_bins': Richness_bin}
    fit_cosmo_str = 'No' if fit_cosmo=='False' else 'Yes'
    #############_Data_#############

    obs = _load_data.load_data(analysis_metadata, Z_bin, Richness_bin, z_corner, rich_corner)
    
    if analysis_metadata['type'] == 'N': N_obs = obs
    if analysis_metadata['type'] == 'M' or analysis_metadata['type'] == 'MxN': N_obs, log10Mass, log10Mass_err = obs
    if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN': N_obs, DS_obs, Err_obs, mask_radius=obs

    ##############################################################################################################################
    # #############################################################################################################################
    
    # cosmology
    Omega_m_true = 0.2648
    Omega_b_true = 0.0448
    Omega_c_true = Omega_m_true - Omega_b_true
    sigma8_true = 0.8
    H0_true = 71
    ns_true = 0.963
    True_value = [Omega_m_true, sigma8_true]
    cosmo_clmm = Cosmology(H0 = H0_true, Omega_dm0 = Omega_c_true, Omega_b0 = Omega_b_true, Omega_k0 = 0.0)
    cosmo = ccl.Cosmology(Omega_c = Omega_c_true, Omega_b = Omega_b_true, h = H0_true/100, sigma8 = sigma8_true, n_s=ns_true)
    #halo model
    massdef = ccl.halos.massdef.MassDef(200, 'critical',)
    if analysis_metadata['hmf'] == 'Bocquet16':
        hmd = ccl.halos.hmfunc.MassFuncBocquet16(mass_def=massdef,hydro=False)
    elif analysis_metadata['hmf'] == 'Bocquet20':
        hmd = ccl.halos.hmfunc.MassFuncBocquet20(mass_def=massdef)
    elif analysis_metadata['hmf'] == 'Despali16':  
        hmd = ccl.halos.hmfunc.MassFuncDespali16(mass_def=massdef)
    #purity
    a_nc, b_nc, a_rc, b_rc = np.log(10)*0.8612, np.log(10)*0.3527, 2.2183, -0.6592
    theta_purity = [a_nc, b_nc, a_rc, b_rc]
    #completeness
    a_nc, b_nc, a_mc, b_mc = 1.1321, 0.7751, 13.31, 0.2025
    theta_completeness = [a_nc, b_nc, a_mc, b_mc]
    #rm_relation (first estimation)
    log10m0, z0 = np.log10(10**14.3), .5
    proxy_mu0, proxy_muz, proxy_mulog10m =  3.34197974,  0.08931269,  2.25997571
    proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m =  0.56014799, -0.05721073,  0.05623336
    theta_rm = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]

    richness_grid = np.logspace(np.log10(20), np.log10(200), 150)
    logm_grid = np.linspace(12, 15.5, 151)
    z_grid = np.linspace(.2, 1, 152)

    grids = {'logm_grid': logm_grid, 'z_grid': z_grid, 'richness_grid':richness_grid}
    count_modelling = {'dNdzdlogMdOmega':None,'richness_mass_relation':None, 'completeness':None, 'purity':None}
    params = {'params_purity':theta_purity, 'params_completeness': theta_completeness, 'params_richness_mass_relation': theta_rm,
             'CCL_cosmology': cosmo, 'halo_mass_distribution': hmd, 'params_concentration_mass_relation': 'Duffy08'}

    compute = {'compute_dNdzdlogMdOmega':True,
               'compute_richness_mass_relation':True,
               'compute_completeness':True, 
               'compute_dNdzdlogMdOmega_log_slope':False,
               'compute_purity':True , 
               'compute_halo_bias':True}
    
    logger.info('[load theory]: compute HMF+bias mass-redshift grids at fixed cosmology')
    count_modelling_new = cl_count.recompute_count_modelling(count_modelling, grids = grids, compute = compute, params = params)

    if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN':

        logger.info('[load theory]: compute mass-redshift lensing mass-redshift profiles')

        rup = analysis_metadata['radius_max']
        rlow = analysis_metadata['radius_min']
        r_edges = _config_lensing_profiles.bin_edges
        #linear
        #r = np.array([(r_edges[i]+r_edges[i+1])/2 for i in range(len(r_edges)-1)])
        #quadratic
        r = np.array([(2/3)*(r_edges[i+1]**3-r_edges[i]**3)/(r_edges[i+1]**2-r_edges[i]**2) for i in range(len(r_edges)-1)])
        #
        r_reshaped = r[(r > rlow)*(r < rup)]
        two_halo_bool = True if analysis_metadata['two_halo']=='True' else False
        cluster_lensing = cl_lensing.compute_cluster_lensing(r_reshaped, analysis_metadata['cM_relation'],
                                                             logm_grid, z_grid, cosmo, cosmo_clmm,
                                                            two_halo = two_halo_bool)
        logger.info('[load theory]: Selection bias: Shear-richness covariance (Ben Zhang)')
        cov, cov_err = GammaLambda_Cov.read_covariance(Richness_bin, Z_bin, photoz = analysis_metadata['photoz'], 
                                                       radius = r)
        cov_DS_selection_bias_obs_mask = cov[mask_radius]
        cov_err_DS_selection_bias_obs_mask = cov_err[mask_radius]
        cov_DS_selection_bias_obs = cov_DS_selection_bias_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
        cov_err_DS_selection_bias_obs = cov_err_DS_selection_bias_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
        #convert to good units (physical units)
        GammaCov = cov_DS_selection_bias_obs * (H0_true/100) * 10**11
        GammaCov_err = cov_err_DS_selection_bias_obs * (H0_true/100) * 10**11

        logger.info('[load theory]: log-slope of the halo mass function (Ben Zhang)')
        Gamma1 = GammaLambda_Cov.read_gamma(Richness_bin, Z_bin, Gamma1_name = 'HMF_Despali16_Payerne_scaling_beta1_new.pkl')

    if analysis_metadata['type'] == 'N' or analysis_metadata['type'] == 'WLxN' or analysis_metadata['type'] == 'MxN':
        logger.info('[load theory]: Compute Sij matrix (SSC) from PySSC (Lacasa et al.)')
        f_sky = (439.78987/(360**2/np.pi))
        CLCovar = cl_covar.Covariance_matrix()
        Sij_partialsky = CLCovar.compute_theoretical_Sij(Z_bin, cosmo, f_sky)

    def prior(theta, add_bias='False'):
        if add_bias=='False':
            proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = theta
        else:
            proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m, b = theta
            if b < -0.5: return -np.inf
            if b > 0.5: return -np.inf
        #mean parameter priors
        if proxy_mu0 < 0: return -np.inf
        #
        if proxy_muz < -2: return -np.inf
        if proxy_muz > 2: return -np.inf
        #
        if proxy_mulog10m < 0: return -np.inf
        #dispersion parameter priors
        if proxy_sigma0 < 0: return -np.inf
        #
        if proxy_sigmaz < -2: return -np.inf
        if proxy_sigmaz > 2: return -np.inf
        #
        if proxy_sigmalog10m < -2: return -np.inf
        if proxy_sigmalog10m > 2: return -np.inf

    def lnL(theta):
        
        if analysis_metadata['add_bias_lensing']=='False':
            proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = theta
        else:
            proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m, b = theta

        prior(theta, add_bias=analysis_metadata['add_bias_lensing'])

        theta_rm_new = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]

        if fit_cosmo == True:
            cosmo_new = ccl.Cosmology(Omega_c = Om-Omega_b_true, Omega_b = Omega_b_true, h = H0_true/100, sigma8 = s8, n_s=ns_true)
            hmd_new = ccl.halos.hmfunc.MassFuncBocquet16(cosmo_new, mass_def=massdef,hydro=False)
        else: 
            hmd_new = hmd
            cosmo_new = cosmo

        params_new = {'params_purity':theta_purity, 
                      'params_completeness': theta_completeness, 
                      'params_richness_mass_relation': theta_rm_new,
                      'CCL_cosmology': cosmo_new, 
                      'halo_mass_distribution': hmd_new}

        compute_new = {'compute_dNdzdlogMdOmega':fit_cosmo,
                       'compute_richness_mass_relation':True, 
                       'compute_dNdzdlogMdOmega_log_slope':False,
                       'compute_completeness':False, 
                       'compute_purity':False,
                       'compute_halo_bias':False}

        adds_N = {'add_purity':True, 
                  'add_completeness':True}
        adds_NDS = {'add_purity':False, 
                    'add_completeness':True}

        count_modelling_new = cl_count.recompute_count_modelling(count_modelling, grids = grids, compute = compute_new, params = params_new)
        test_sign = count_modelling_new['richness_mass_relation - sigma'].flatten() < 0
        if len(test_sign[test_sign==True]) != 0: return -np.inf
        integrand_count_new = cl_count.define_count_integrand(count_modelling_new, adds_N)
        integrand_count_ds_new = cl_count.define_count_integrand(count_modelling_new, adds_NDS)
        Omegaredmapper = 439.78987
        Omega = 4*np.pi*(Omegaredmapper/(360**2/np.pi))

        N = Omega*cl_count.Cluster_SurfaceDensity_ProxyZ(bins, integrand_count = integrand_count_new, grids = grids)
        if analysis_metadata['type'] == 'N' or analysis_metadata['type'] == 'WLxN' or analysis_metadata['type'] == 'MxN':
            gaussian=True
            if gaussian:
                NAverageHaloBias = Omega * cl_count.Cluster_NHaloBias_ProxyZ(bins, integrand_count = integrand_count_new,
                                                                             halo_bias = count_modelling_new['halo_bias'], 
                                                                             grids = grids, cosmo = cosmo)
                CLCovar = cl_covar.Covariance_matrix()
                NNSbb = CLCovar.sample_covariance_full_sky(Z_bin, Richness_bin, 
                                                              NAverageHaloBias.T, 
                                                              Sij_partialsky)
                Cov_tot = NNSbb + np.diag(N.T.flatten())
                CLCount_likelihood.lnLikelihood_Binned_Gaussian(N, N_obs.T, Cov_tot)
                lnLCLCount = CLCount_likelihood.lnL_Binned_Gaussian

            else:
                CLCount_likelihood.lnLikelihood_Binned_Poissonian(N, N_obs.T)
                lnLCLCount = CLCount_likelihood.lnL_Binned_Poissonian

        if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN':
            NDS_profiles = Omega * cl_lensing.Cluster_dNd0mega_Lensing_ProxyZ(bins, integrand_count = integrand_count_ds_new, 
                                                                              cluster_lensing = cluster_lensing, 
                                                                              lensing_radius = r_reshaped, grids = grids)
            DS_profiles = NDS_profiles/N
            if analysis_metadata['add_bias_lensing']=='True': 
                DS_profiles = (1+b)*DS_profiles
            if analysis_metadata['shear_richness_cov']=='True': 
                factor = (np.log(10) * Gamma1 / proxy_mulog10m)
                DS_profiles = DS_profiles + factor * GammaCov 
                Err_obs2_tot = Err_obs**2 + factor**2 * GammaCov_err**2
                Err_obs_tot = Err_obs2_tot**.5
                Sum_logSigma = np.sum(np.log(Err_obs_tot.flatten()**2))
            else: 
                Err_obs_tot = Err_obs
                Sum_log = 0
                

            lnLCLWL = -0.5*np.sum(((DS_profiles - DS_obs)/Err_obs_tot)**2) - 0.5*Sum_logSigma

        if analysis_metadata['type'] == 'M' or analysis_metadata['type'] == 'MxN':
            NM = Omega * cl_mass.Cluster_dNd0mega_Mass_ProxyZ(bins, integrand_count = integrand_count_ds_new, grids = grids)
            Mth = NM/N
            log10Mth = np.log10(Mth)
            lnLCLM = -.5*np.sum(((log10Mth - log10Mass)/log10Mass_err)**2)

        if analysis_metadata['type'] == 'WL': return lnLCLWL 
        if analysis_metadata['type'] == 'N' : return lnLCLCount
        if analysis_metadata['type'] == 'M' : return lnLCLM
        if analysis_metadata['type'] == 'WLxN': return lnLCLWL + lnLCLCount
        if analysis_metadata['type'] == 'MxN': return lnLCLM + lnLCLCount

    if analysis_metadata['add_bias_lensing']=='False': 
        initial = [proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]
        labels = [r'\ln \lambda_0', r'\mu_z', r'\mu_m', r'\sigma_{\ln \lambda, 0}', r'\sigma_z', r'\sigma_m']
        ndim=6
    else: 
        b=0
        initial = [proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m, b]
        labels = [r'\ln \lambda_0', r'\mu_z', r'\mu_m', r'\sigma_{\ln \lambda, 0}', r'\sigma_z', r'\sigma_m', 'b']
        ndim=7
    t = time.time()
    logger.info(lnL(initial))
    tf = time.time()
    logger.info(tf-t)
    nwalker = 100
    pos = np.array(initial) + .01*np.random.randn(nwalker, ndim)
    sampler = emcee.EnsembleSampler(nwalker, ndim, lnL)
    sampler.run_mcmc(pos, 200, progress=True);
    flat_samples = sampler.get_chain(discard=0, thin=1, flat=True)
    results={'flat_chains':flat_samples, 'analysis':analysis_metadata, 'label_parameters':labels}
    save_pickle(results, f'../chains/'+analysis_metadata['type']+'/'+analysis_metadata['name']+'.pkl')
