import sys
import pyccl as ccl
import numpy as np
from clmm import Cosmology
from multiprocessing import Pool
import emcee
import matplotlib.pyplot as plt
import time, os
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
    fit_cosmo = True if analysis_metadata['fit_cosmo'] == 'True' else 'False'

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
    fit_cosmo_str = 'No' if fit_cosmo==False else 'Yes'
    #############_Data_#############

    obs = _load_data.load_data_vary_cosmology(analysis_metadata)
    
    if analysis_metadata['type'] == 'N': N_obs = obs
    if analysis_metadata['type'] == 'M' or analysis_metadata['type'] == 'MxN': 
        N_obs, log10Mass, log10Mass_err = obs

    if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN': 
        _, N_obs, DS_obs, Err_obs, mask_radius=obs

    ##############################################################################################################################
    # #############################################################################################################################
    
    # cosmology
    H0_true = 71
    h = H0_true/100
    Omega_b_true = 0.02258 / (h**2)
    Omega_c_true = 0.1109 / (h**2)
    Omega_m_true = Omega_b_true + Omega_c_true
    sigma8_true = 0.8
    ns_true = 0.963
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
    proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m =  0.56014799, 0,  0
    theta_rm = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, 0, 0]

    richness_grid = np.logspace(np.log10(20), np.log10(200), 150)
    logm_grid = np.linspace(12, 15.5, 20)
    z_grid = np.geomspace(.2, 1, 50)

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
        
        #logger.info('[load theory]: Selection bias: Shear-richness covariance (Ben Zhang)')
        #cov, cov_err = GammaLambda_Cov.read_covariance(Richness_bin, Z_bin, photoz = analysis_metadata['photoz'],radius = r)
        #cov_DS_selection_bias_obs_mask = cov[mask_radius]
        #cov_err_DS_selection_bias_obs_mask = cov_err[mask_radius]
        #cov_DS_selection_bias_obs = cov_DS_selection_bias_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
        #cov_err_DS_selection_bias_obs = cov_err_DS_selection_bias_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
        #convert to good units (physical units)
        #GammaCov = cov_DS_selection_bias_obs * (H0_true/100) * 10**11
        #GammaCov_err = cov_err_DS_selection_bias_obs * (H0_true/100) * 10**11
        #logger.info('[load theory]: log-slope of the halo mass function (Ben Zhang)')
        #Gamma1 = GammaLambda_Cov.read_gamma(Richness_bin, Z_bin, Gamma1_name = 'HMF_Despali16_Payerne_scaling_beta1_new.pkl')

    if analysis_metadata['type'] == 'N' or analysis_metadata['type'] == 'WLxN' or analysis_metadata['type'] == 'MxN':
        if analysis_metadata['Gauss+SSC-CC_likelihood'] == 'True':
            logger.info('[load theory]: Compute Sij matrix (SSC) from PySSC (Lacasa et al.)')
            f_sky = (439.78987/(360**2/np.pi))
            CLCovar = cl_covar.Covariance_matrix()
            Sij_partialsky = CLCovar.compute_theoretical_Sij(Z_bin, cosmo, f_sky)

    def prior_no_cosmo(theta, add_bias='False'):
        if add_bias=='False':
            proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0 = theta
        else:
            proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, b = theta
            if b < -0.5: return -1
            if b > 0.5: return -1
        #mean parameter priors
        if proxy_mu0 < 0: return -1
        #
        if proxy_muz < -2: return -1
        if proxy_muz > 2: return -1
        #
        if proxy_mulog10m < 0: return -1
        #dispersion parameter priors
        if proxy_sigma0 < 0: return -1

        return 0

    def prior_cosmo(theta_cosmo, params = 'Om_s8'):

        if params == 'Om_s8':
            Om, s8 = theta_cosmo
            if Om < 0.2: return -1
            if Om > 0.6: return -1
            if s8 < 0.2: return -1
            if s8 > 1: return -1

        if params == 'Om_w':
            Om, w = theta_cosmo
            if Om < 0.2: return -1
            if Om > 0.6: return -1
            if w < -2.5: return -1
            if w > 0: return -1

        return 0

    def lnL(theta):
        
        if not fit_cosmo:
            
            if analysis_metadata['add_bias_lensing']=='False':
                
                proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0 = theta_no_cosmo
            else:
                proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, b = theta_no_cosmo
    
            p = prior_no_cosmo(theta_no_cosmo, add_bias=analysis_metadata['add_bias_lensing'])
            if p == -1: return -np.inf
    
            theta_rm_new = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0]

            hmd_new = hmd
            
            cosmo_new = cosmo

        if fit_cosmo:
                
            if analysis_metadata['add_bias_lensing']=='False':
                
                COSMO1, COSMO2, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0 = theta
                theta_no_cosmo = proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0
                
            else:
                
                COSMO1, COSMO2, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, b = theta
                theta_no_cosmo = proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, b

            p1 = prior_no_cosmo(theta_no_cosmo, add_bias=analysis_metadata['add_bias_lensing'])
            if p1 == -1: return -np.inf
            p2 = prior_cosmo([COSMO1, COSMO2], params = analysis_metadata['cosmo_params'])
            if p2 == -1: return -np.inf

            #change cosmology
            if analysis_metadata['cosmo_params'] == 'Om_s8':
                Om, s8 = COSMO1, COSMO2
                cosmo_new = ccl.Cosmology(Omega_c = Om-Omega_b_true, Omega_b = Omega_b_true, 
                                          h = H0_true/100, sigma8 = s8, n_s=ns_true, w0 = -1, wa = 0)

            if analysis_metadata['cosmo_params'] == 'Om_w':
                Om, w = COSMO1, COSMO2
                cosmo_new = ccl.Cosmology(Omega_c = Om-Omega_b_true, Omega_b = Omega_b_true, 
                                          h = H0_true/100, sigma8 = sigma8_true, n_s=ns_true, w0 = w, wa = 0)

        theta_rm_new = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, 0, 0]

        params_new = {'params_purity':theta_purity, 
                      'params_completeness': theta_completeness, 
                      'params_richness_mass_relation': theta_rm_new,
                      'CCL_cosmology': cosmo_new, 
                      'halo_mass_distribution': hmd}

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

            if analysis_metadata['Gauss+SSC-CC_likelihood'] == 'True':
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

            if fit_cosmo:
                #cosmo_clmm.be_cosmo = cosmo_new
                cluster_lensing_new = cluster_lensing#cl_lensing.compute_cluster_lensing(r_reshaped, 
                                                                  #       analysis_metadata['cM_relation'],
                                                                  #       logm_grid, z_grid, cosmo_new, 
                                                                  #       cosmo_clmm,
                                                                   #      two_halo = two_halo_bool)
            else:
                cluster_lensing_new = cluster_lensing
            
            NDS_profiles = Omega * cl_lensing.Cluster_dNd0mega_Lensing_ProxyZ(bins, integrand_count = integrand_count_ds_new, 
                                                                              cluster_lensing = cluster_lensing_new, 
                                                                              lensing_radius = r_reshaped, grids = grids)
            DS_profiles = NDS_profiles/N

            Err_obs_tot = Err_obs(Om)
            Sum_logSigma = np.sum(np.log(Err_obs_tot.flatten()**2))

            lnLCLWL = -0.5*np.sum(((DS_profiles - DS_obs(Om))/Err_obs_tot)**2) - 0.5*Sum_logSigma

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


    initial_rm = [proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0]
    initial_b = [0]
    initial_Om_s8 = [Omega_m_true, sigma8_true]
    initial_Om_w = [Omega_m_true, -1]
    
    label_rm = [r'\ln \lambda_0', r'\mu_z', r'\mu_m', r'\sigma_{\ln \lambda, 0}']
    label_b = ['b']
    label_Om_s8 = [r'\Omega_m', '\sigma_8']
    label_Om_w = [r'\Omega_m', 'w']

    if not fit_cosmo:
        initial = initial_rm
        labels = label_rm
        ndim=4
        
    else:
        if analysis_metadata['cosmo_params'] == 'Om_s8': 
            initial_cosmo = initial_Om_s8
            label_cosmo = label_Om_s8
            
        if analysis_metadata['cosmo_params'] == 'Om_w': 
            initial_cosmo = initial_Om_w
            label_cosmo = label_Om_w
            
        initial = initial_cosmo + initial_rm
        labels = label_cosmo + label_rm
        ndim=6

    
    t = time.time()
    logger.info(lnL(initial))
    tf = time.time()
    likelihood_time = tf-t
    logger.info('likelihood evaluation time = ' + str(tf-t) + ' seconds')
    nwalker = 100
    npath = 200
    total_time = nwalker * npath * likelihood_time
    logger.info('total run time = ' + str(total_time/(60*60)) + ' hours')
    pos = np.array(initial) + .01*np.random.randn(nwalker, ndim)
    pos = pos[(pos[:,0] > 0.2)*(pos[:,0] < 0.6)*(pos[:,1] > 0.2)*(pos[:,1] < 1)*(pos[:,5] > 0)]
    nwalker = len(pos)

    filename = f"../chains/{analysis_metadata['type']}/{analysis_metadata['name']}.h5"
    backend = emcee.backends.HDFBackend(filename)
    
    if not os.path.exists(filename):
        # First run: reset backend and start fresh
        backend.reset(nwalker, ndim)
        p0 = pos
    else:
        # Resume run from last saved position
        try:
            p0 = backend.get_last_sample().coords
            nwalker = len(p0)
            logger.info("Resuming from last saved state in %s", filename)
        except Exception:
            # If file exists but is empty/corrupted → restart
            backend.reset(nwalker, ndim)
            p0 = pos
            logger.info("Backend file empty or corrupted. Restarting fresh.")
    
    # --- Run sampler ---
    sampler = emcee.EnsembleSampler(nwalker, ndim, lnL, backend=backend)
    
    nsteps = npath  # total steps to run in this execution
    sampler.run_mcmc(p0, nsteps, progress=True)
    #sampler = emcee.EnsembleSampler(nwalker, ndim, lnL)
    #sampler.run_mcmc(pos, npath, progress=True);
    #flat_samples = sampler.get_chain(discard=0, thin=1, flat=True)
    #results={'flat_chains':flat_samples, 'analysis':analysis_metadata, 'label_parameters':labels}
    #save_pickle(results, f'../chains/'+analysis_metadata['type']+'/'+analysis_metadata['name']+'.pkl')
