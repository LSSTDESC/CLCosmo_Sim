import sys
import pyccl as ccl
import numpy as np
from clmm import Cosmology
from multiprocessing import Pool
import emcee
import time
import pickle
import analysis_list
def save_pickle(dat, filename, **kwargs):
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()
def load(filename, **kwargs):
    """Loads GalaxyCluster object to filename using Pickle"""
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)

import redshift_richness_bins as analysis
sys.path.append('/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/cluster_abundance')
import CL_COUNT_modeling_completeness as comp
import CL_COUNT_modeling_purity as pur
import CL_COUNT_modeling_halo_mass_function as hmf
import CL_COUNT_modeling_richness_mass_relation as rm_relation
import CL_MASS_cluster_mass as cl_mass
import CL_COUNT_cluster_abundance as cl_count
import CL_COUNT_class_likelihood as likelihood
import CL_LENSING_cluster_lensing as cl_lensing
CLCount_likelihood = likelihood.Likelihood()


code, index_analysis = sys.argv
analysis_metadata = analysis_list.analysis[int(index_analysis)]
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
log10m0, z0 = np.log10(10**14.2), .4
proxy_mu0, proxy_muz, proxy_mulog10m =  3.091, 0, 1.05*np.log(10)
proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m =  0.594, 0., 0.026*np.log(10)
theta_rm = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]

richness_grid = np.logspace(np.log10(20), np.log10(200), 100)
logm_grid = np.linspace(12, 15.5, 101)
z_grid = np.linspace(.2, 1, 102)

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
table_redmapper = load('/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/richness_mass_from_DC2/catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl')
N_obs, proxy_edges, z_edges = np.histogram2d(table_redmapper['redshift'], 
                                                        table_redmapper['richness'],
                                                   bins=[analysis.z_corner, analysis.rich_corner])

#Masses
mass_data = load('/pbs/throng/lsst/users/cpayerne/CLMassDC2/notebooks/plots/WL_mean_masses/halomodel_nfw_freec.pkl')['masses']
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

    #cluster_lensing = np.zeros([len(r), len(logm_grid), len(z_grid)])#
    cluster_lensing = cl_lensing.compute_cluster_lensing(r_reshaped, analysis_metadata['cM_relation'],
                                                         logm_grid, z_grid, cosmo, cosmo_clmm,
                                                        two_halo = analysis_metadata['two_halo'])
    print('lensing_profile_computed')

print('start')

    
def lnL(theta):
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
    
    if fit_cosmo == True:
        compute_new = {'compute_dNdzdlogMdOmega':True,
                       'compute_richness_mass_relation':True, 
                       'compute_completeness':False, 
                       'compute_purity':False}
    else:
        compute_new = {'compute_dNdzdlogMdOmega':False,
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
    CLCount_likelihood.lnLikelihood_Binned_Poissonian(N, N_obs.T)
    lnLCLCount = CLCount_likelihood.lnL_Binned_Poissonian
    #
    if analysis_metadata['type'] == 'WL' or analysis_metadata['type'] == 'WLxN':
        NDS_profiles = Omega * cl_lensing.Cluster_dNd0mega_Lensing_ProxyZ(bins, integrand_count = integrand_count_ds_new, 
                                                                          cluster_lensing = cluster_lensing, 
                                                                          lensing_radius = r_reshaped, grids = grids)
       
        DS_profiles = NDS_profiles/N
        lnLCLWL = -.5*np.sum(((DS_profiles - DS_obs)/Err_obs)**2)
    #
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
   

if fit_cosmo == True:
    initial = [proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m, Omega_m_true, sigma8_true]
    ndim = 8
else: 
    initial = [proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]
    ndim = 6
t = time.time()
print(lnL(initial))
tf = time.time()
print(tf-t)
nwalker = 100
pos = np.array(initial) + .01*np.random.randn(nwalker, ndim)
#pos = np.array(initial) + .01*np.random.randn(nwalker, 8)
#with Pool() as pool:
#sampler = emcee.EnsembleSampler(nwalker, 6, lnL, )#pool=pool)
sampler = emcee.EnsembleSampler(nwalker, ndim, lnL,)# pool=pool)
sampler.run_mcmc(pos, 300, progress=True);
flat_samples = sampler.get_chain(discard=0, thin=1, flat=True)
name = analysis_metadata['name']
type_of_fit = analysis_metadata['type']
if fit_cosmo== True: add = '_fit_cosmo'
else: add = ''
save_pickle(flat_samples, f'/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/richness_mass_from_DC2/chains_article/{type_of_fit}_{name}' + add + '.pkl')
#save_pickle(flat_samples, f'/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/richness_mass_from_DC2/chains/manuscript_analysis_full_cosmology_{t}.pkl')