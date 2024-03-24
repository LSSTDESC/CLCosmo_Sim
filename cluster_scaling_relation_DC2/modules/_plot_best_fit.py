import sys
import pyccl as ccl
import numpy as np
from clmm import Cosmology
from multiprocessing import Pool
import emcee
import time
import matplotlib.pyplot as plt
import pickle
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

sys.path.append('../../modeling')
import CL_COUNT_modeling_completeness as comp
import CL_COUNT_modeling_purity as pur
import CL_COUNT_modeling_halo_mass_function as hmf
import CL_COUNT_modeling_richness_mass_relation as rm_relation
import CL_MASS_cluster_mass as cl_mass
import CL_COUNT_cluster_abundance as cl_count
import CL_COUNT_class_likelihood as likelihood
import CL_LENSING_cluster_lensing as cl_lensing
CLCount_likelihood = likelihood.Likelihood()

code, config_name, index_analysis = sys.argv
analysis_metadata = analysis_list.config[config_name][int(index_analysis)]
fit_cosmo = analysis_metadata['fit_cosmo']

#cosmology
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
massdef = ccl.halos.massdef.MassDef(200, 'critical',)# c_m_relation=None)
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
#rm_relation
log10m0, z0 = np.log10(10**14.3), .5
proxy_mu0, proxy_muz, proxy_mulog10m =  3.091, 0, 1.05*np.log(10)
proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m =  0.594, 0., 0.
theta_rm = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]

richness_grid = np.logspace(np.log10(20), np.log10(200), 150)
logm_grid = np.linspace(12, 15.5, 151)
z_grid = np.linspace(.2, 1, 152)

grids = {'logm_grid': logm_grid, 'z_grid': z_grid, 'richness_grid':richness_grid}
count_modelling = {'dNdzdlogMdOmega':None,'richness_mass_relation':None, 'completeness':None, 'purity':None }
params = {'params_purity':theta_purity, 'params_completeness': theta_completeness, 'params_richness_mass_relation': theta_rm,
         'CCL_cosmology': cosmo, 'halo_mass_distribution': hmd}

compute = {'compute_dNdzdlogMdOmega':True,'compute_richness_mass_relation':True, 'compute_completeness':True, 'compute_purity':True }

count_modelling_new = cl_count.recompute_count_modelling(count_modelling, grids = grids, compute = compute, params = params)

Z_bin = analysis.Z_bin
Richness_bin = analysis.Obs_bin
bins = {'redshift_bins':Z_bin, 'richness_bins': Richness_bin}

#############_Data_#############
# Lensing
data = np.load(analysis_metadata['lensing_data'], allow_pickle=True)
profiles = data['stacked profile']
covariances = data['stacked covariance']
r = profiles['radius'][0]

# Count
table_redmapper = load('../../data/lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl')
N_obs, proxy_edges, z_edges = np.histogram2d(table_redmapper['redshift'], 
                                                        table_redmapper['richness'],
                                                   bins=[analysis.z_corner, analysis.rich_corner])

#Masses
if (analysis_metadata['type']=='M') or (analysis_metadata['type']=='MxN') or (analysis_metadata['type']=='M_DC2xN') or (analysis_metadata['type']=='M_DC2'):
    mass_data = load(analysis_metadata['mass_file'])['masses']
    log10Mass = np.zeros(N_obs.T.shape)
    log10Mass_err = np.zeros(N_obs.T.shape)
    for i, richness_bin in enumerate(Richness_bin):
        for j, redshift_bin in enumerate(Z_bin):
            maskz = (mass_data['z_mean'] > redshift_bin[0])*(mass_data['z_mean'] < redshift_bin[1])
            maskr = (mass_data['obs_mean'] > richness_bin[0])*(mass_data['obs_mean'] < richness_bin[1])
            mask = maskz * maskr
            mass_in_bin = mass_data[mask]
            log10Mass[i,j] = np.array(mass_in_bin['log10M200c_WL'])
            log10Mass_err[i,j] = np.array(mass_in_bin['err_log10M200c_WL'])

#lensing profiles
if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN':
    
    DS_obs = np.zeros([len(r), len(Richness_bin), len(Z_bin)])
    Err_obs = np.zeros([len(r), len(Richness_bin), len(Z_bin)])

    for i, z_bin in enumerate(Z_bin):
        mask_z = (profiles['z_mean'] > z_bin[0])*(profiles['z_mean'] < z_bin[1])
        for j, richness_bin in enumerate(Richness_bin):
            mask_richness = (profiles['obs_mean'] > richness_bin[0])*(profiles['obs_mean'] < richness_bin[1])
            mask_tot = mask_z * mask_richness
            DS_obs[:,j,i] = profiles['gt'][mask_tot][0]
            Err_obs[:,j,i] = covariances['cov_t'][mask_tot][0].diagonal()**.5
    DS_obs_full = DS_obs
    Err_obs_full = Err_obs

    r_grid = np.zeros(DS_obs.shape)
    for i, z_bin in enumerate(Z_bin):
        for j, richness_bin in enumerate(Richness_bin):
            r_grid[:,j,i] = r

    rup = analysis_metadata['radius_max']
    rlow = analysis_metadata['radius_min']
    mask = (r_grid > rlow)*(r_grid < rup)

    DS_obs_mask = DS_obs[mask]
    DS_obs = DS_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])

    Err_obs_mask = Err_obs[mask]
    Err_obs = Err_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])

    r_reshaped = r[(r > rlow)*(r < rup)]
    print('-----     observed lensing profiles loaded + covariances')
    cluster_lensing = cl_lensing.compute_cluster_lensing(r_reshaped, analysis_metadata['cM_relation'],
                                                         logm_grid, z_grid, cosmo, cosmo_clmm,
                                                        two_halo = analysis_metadata['two_halo'])
    print('-----     theoretical lensing profiles computed')
    
    cov, cov_err = GammaLambda_Cov.read_covariance(photoz = analysis_metadata['photoz'], radius = r)
    cov_DS_selection_bias_obs_mask = cov[mask]
    cov_err_DS_selection_bias_obs_mask = cov_err[mask]
    cov_DS_selection_bias_obs = cov_DS_selection_bias_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
    cov_err_DS_selection_bias_obs = cov_err_DS_selection_bias_obs_mask.reshape([len(r[(r > rlow)*(r < rup)]), len(Richness_bin), len(Z_bin)])
    print('-----     shear-richness covariances loaded')
    
    Gamma1 = GammaLambda_Cov.read_gamma()
    print('-----     log-slope of the halo mass function loaded')
    
    GammaCov = np.zeros(cov_DS_selection_bias_obs.shape)
    for i in range(len(Richness_bin)):
        for j in range(len(Z_bin)):
            GammaCov[:,i,j] = cov_DS_selection_bias_obs[:,i,j] * Gamma1[i,j]
    print('-----     \gamma_1 * Cov(\Delta\Sigma, \lambda|m,z)')

    print('start')

    
def Prediction(theta):
    if fit_cosmo == True:
        proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m, Om, s8 = theta
        if proxy_sigma0 < 0: return -np.inf
        if Om > .5: return -np.inf
        if Om < .1: return -np.inf
        if s8 > 1: return -np.inf
        if s8 < .5: return -np.inf
    else: 
        proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = theta
        if proxy_sigma0 < 0: return -np.inf
    
    theta_rm_new = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]
    #theta_rm_new = [log10m0, z0, proxy_mu0, 0, proxy_mulog10m, proxy_sigma0, 0, proxy_sigmalog10m]
    
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
                   'compute_completeness':False, 
                   'compute_purity':False}
        
    adds_N = {'add_purity':True, 'add_completeness':True}
    adds_NDS = {'add_purity':False, 'add_completeness':True}

    count_modelling_new = cl_count.recompute_count_modelling(count_modelling, grids = grids, compute = compute_new, params = params_new)
    test_sign = count_modelling_new['richness_mass_relation - sigma'].flatten() < 0
    if len(test_sign[test_sign==True]) != 0: return -np.inf
    integrand_count_new = cl_count.define_count_integrand(count_modelling_new, adds_N)
    integrand_count_ds_new = cl_count.define_count_integrand(count_modelling_new, adds_NDS)
    Omegaredmapper = 439.78987
    Omega = 4*np.pi*(Omegaredmapper/(360**2/np.pi))
    N = Omega * cl_count.Cluster_SurfaceDensity_ProxyZ(bins, integrand_count = integrand_count_new, grids = grids)
    #
    if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN':
        NDS_profiles = Omega * cl_lensing.Cluster_dNd0mega_Lensing_ProxyZ(bins, integrand_count = integrand_count_ds_new, 
                                                                          cluster_lensing = cluster_lensing, 
                                                                          lensing_radius = r_reshaped, grids = grids)
       
        DS_profiles = NDS_profiles/N
        if analysis_metadata['shear_richness_cov']: DS_profiles = DS_profiles + GammaCov/proxy_mulog10m
    #
    if analysis_metadata['type'] == 'M' or analysis_metadata['type'] == 'MxN':
        NM = Omega * cl_mass.Cluster_dNd0mega_Mass_ProxyZ(bins, integrand_count = integrand_count_ds_new, grids = grids)
        Mth = NM/N
        log10Mth = np.log10(Mth)
        
    if analysis_metadata['type'] == 'WL': 
        return {'WL':DS_profiles}
    if analysis_metadata['type'] == 'N' : 
        return {'N':N}
    if analysis_metadata['type'] == 'M' : 
        return {'M':log10Mth}
    if analysis_metadata['type'] == 'WLxN': 
        return {'WL':DS_profiles,'N':N}
    if analysis_metadata['type'] == 'MxN': 
        return {'M':log10Mth,'N':N}

n_cut = 15000
result_mcmc = np.load(f'../chains/'+analysis_metadata['type']+'/'+analysis_metadata['name']+'.pkl', allow_pickle=True)
chains = result_mcmc['flat_chains'][n_cut:]
params = np.mean(chains, axis=0)
fit = Prediction(params)

if 'N' in analysis_metadata['type']:

    fig, ax = plt.subplots(1, len(Z_bin), figsize=(10,3), sharex = True, sharey = True)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)

    for i, z_bin in enumerate(Z_bin):
        z0 = z_bin[0]
        z1 = z_bin[1]
        ax[i].set_title(f'{z0:.1f} < z < {z1:.1f}', fontsize=10)
        ax[i].plot(np.log(np.mean(Richness_bin, axis=1)), fit['N'][:,i], ls='-', 
                   label = r'Best fit', color = 'C2')
        ax[i].errorbar(np.log(np.mean(Richness_bin, axis=1)), N_obs.T[:,i],  N_obs.T[:,i]**.5,  
                       color = 'dodgerblue', marker = '.', linestyle='none')
        ax[i].set_xlabel(r'$\ln \lambda$', fontsize=15)
        ax[0].set_ylabel('N', fontsize=20)
        ax[0].legend(fontsize=12, frameon=True,ncols=1, loc = 'upper right')
        ax[i].set_yscale('log')
        ax[i].set_ylim(2, 1e3)
        ax[i].set_xlim(np.log(20), np.log(200))
        ax[i].tick_params(axis='both', which="both", labelsize= 12)
        ax[i].grid(True)
    plt.savefig('../fig/best_fit_N_from_'+analysis_metadata['name']+'.png', bbox_inches='tight', dpi=100)
    
if 'WL' in analysis_metadata['type']:
    n_z_bin = len(Z_bin) 
    n_m_bin = len(Richness_bin) 
    scale = 4
    fig, axs = plt.subplots(n_m_bin, n_z_bin, figsize = (10,5), sharex = True, sharey = True)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
    #fig.subplots_adjust(wspace=0, hspace=0)
    for i, z_bin in enumerate(Z_bin):
        for j, m_bin in enumerate(Richness_bin):
            label_z =   f'{z_bin[0]:.1f} < z < {z_bin[1]:.1f}'
            label_M = f'{m_bin[0]:.0f} < ' + r'$\lambda$' +f' < {m_bin[1]:.0f}'
            axs[j,i].errorbar(r, DS_obs_full[:,j,i], Err_obs_full[:,j,i], 
                              linestyle='none', marker = '.', markersize=3, color='dodgerblue')
            axs[j,i].loglog(r[(r > rlow)*(r < rup)], fit['WL'][:,j,i], color = 'C2')
            axs[j,i].set_xscale('log')
            axs[j,i].set_yscale('log')
            axs[j,i].set_ylim(1e12, 3e14)
            axs[j,i].set_xlim(0.4, 12)
            axs[j,i].tick_params(axis='both', which = 'major', labelsize= 10)
            axs[j,i].legend(frameon = False, loc = 'upper right', fontsize = 10)
            axs[j,i].set_xlabel('R [Mpc]', fontsize = 10)
            z0 = z_bin[0]
            z1 = z_bin[1]
            axs[j,i].set_xlabel('R [Mpc]', fontsize = 10)
            axs[j,0].set_ylabel(label_M, fontsize = 10)
            axs[0,i].set_title(label_z, fontsize = 10)

    for ax in fig.get_axes():
        ax.label_outer()
    plt.savefig('../fig/best_fit_WL_from_'+analysis_metadata['name']+'.png', bbox_inches='tight', dpi=100)