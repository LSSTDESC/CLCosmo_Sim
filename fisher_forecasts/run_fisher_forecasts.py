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

sys.path.append('../modeling')
import CL_COUNT_modeling_completeness as comp
import CL_COUNT_DATAOPS_cluster_abundance_covariance as cl_covar 
import CL_COUNT_modeling_purity as pur
import CL_COUNT_modeling_halo_mass_function as hmf
import CL_COUNT_modeling_richness_mass_relation as rm_relation
import CL_MASS_cluster_mass as cl_mass
import CL_COUNT_cluster_abundance as cl_count
import CL_COUNT_class_likelihood as likelihood
import CL_LENSING_cluster_lensing as cl_lensing

def binning(corner): return [[corner[i],corner[i+1]] for i in range(len(corner)-1)]

z_corner = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1])
Z_bin = np.array(binning(z_corner))
rich_corner = np.array([20, 35, 70, 100, 200])
rich_corner = np.array([r for r in rich_corner])
Obs_bin = np.array(binning(rich_corner))
Obs_bin_center = np.mean(Obs_bin, axis=1)
Z_bin_center = np.mean(Z_bin, axis=1)

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
hmd = ccl.halos.hmfunc.MassFuncBocquet16(mass_def=massdef,hydro=False)
#purity
a_nc, b_nc, a_rc, b_rc = np.log(10)*0.8612, np.log(10)*0.3527, 2.2183, -0.6592
theta_purity = [a_nc, b_nc, a_rc, b_rc]
#completeness
a_nc, b_nc, a_mc, b_mc = 1.1321, 0.7751, 13.31, 0.2025
theta_completeness = [a_nc, b_nc, a_mc, b_mc]
#rm_relation (first estimation)
Omegaredmapper = 440.78987
Omega = 4*np.pi*(Omegaredmapper/(360**2/np.pi))
log10m0, z0 = np.log10(10**14.3), .5
proxy_mu0, proxy_muz, proxy_mulog10m =  3.34197974,  0.08931269,  2.25997571
proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m =  0.56014799, -0.05721073,  0.05623336
theta_rm = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]

richness_grid = np.logspace(np.log10(20), np.log10(200), 250)
logm_grid = np.linspace(12, 15.5, 251)
z_grid = np.linspace(.2, 1, 252)

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
bins = {'redshift_bins':Z_bin, 'richness_bins': Obs_bin}
adds_N = {'add_purity':False, 
              'add_completeness':False}

def N(theta):

    proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = theta

    theta_rm_new = [log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m]

    cosmo_new = ccl.Cosmology(Omega_c = Omega_m_true-Omega_b_true, Omega_b = Omega_b_true, 
                              h = H0_true/100, sigma8 = sigma8_true, n_s=ns_true)
    hmd_new = ccl.halos.hmfunc.MassFuncBocquet16(mass_def=massdef,hydro=False)

    params_new = {'params_purity':theta_purity, 
                  'params_completeness': theta_completeness, 
                  'params_richness_mass_relation': theta_rm_new,
                  'CCL_cosmology': cosmo_new, 
                  'halo_mass_distribution': hmd_new}
    compute_new = {'compute_dNdzdlogMdOmega':True,
                   'compute_richness_mass_relation':True, 
                   'compute_dNdzdlogMdOmega_log_slope':False,
                   'compute_completeness':False, 
                   'compute_purity':False,
                   'compute_halo_bias':False}

    count_modelling_new = cl_count.recompute_count_modelling(count_modelling, grids = grids, compute = compute_new, params = params_new)
    integrand_count_new = cl_count.define_count_integrand(count_modelling_new, adds_N)
    N = Omega*cl_count.Cluster_SurfaceDensity_ProxyZ(bins, integrand_count = integrand_count_new, grids = grids)

    return N

# =========================================================
# =========================================================
# 5-point derivative with explicit flattening
# =========================================================
def compute_dN(theta_fid, eps_rel):
    n_params = len(theta_fid)

    N0 = N(theta_fid)
    N0_flat = N0.T.reshape(-1)
    n_bins = len(N0_flat)

    dN = np.zeros((n_params, n_bins))

    for i in range(n_params):

        # relative step
        step = eps_rel[i] * np.abs(theta_fid[i])
        if step == 0:
            step = eps_rel[i]

        # shifted parameter vectors
        t_m2 = theta_fid.copy()
        t_m1 = theta_fid.copy()
        t_p1 = theta_fid.copy()
        t_p2 = theta_fid.copy()

        t_m2[i] -= 2 * step
        t_m1[i] -= step
        t_p1[i] += step
        t_p2[i] += 2 * step

        # evaluate model (KEEP 2D → then flatten consistently)
        N_m2 = N(t_m2)
        N_m1 = N(t_m1)
        N_p1 = N(t_p1)
        N_p2 = N(t_p2)

        N_m2_flat = N_m2.reshape(-1)
        N_m1_flat = N_m1.reshape(-1)
        N_p1_flat = N_p1.reshape(-1)
        N_p2_flat = N_p2.reshape(-1)

        # 5-point stencil derivative (on flattened vectors)
        dN[i] = (-N_p2_flat + 8*N_p1_flat - 8*N_m1_flat + N_m2_flat) / (12 * step)

    return dN, N0_flat


# =========================================================
# Fisher matrix (Poisson counts)
# =========================================================
def fisher_matrix(theta_fid, eps_rel, Nmin=5):

    dN, N_fid = compute_dN(theta_fid, eps_rel)

    # --- explicit flattening ---
    N_fid_flat = N_fid.reshape(-1)
    dN_flat = dN.reshape(dN.shape[0], -1)

    # mask low-count bins
    mask = N_fid_flat > Nmin

    dN_flat = dN_flat[:, mask]
    N_fid_flat = N_fid_flat[mask]

    # Fisher matrix
    F = dN_flat @ ((1.0 / N_fid_flat)[:, None] * dN_flat.T)

    # symmetrize
    F = F#0.5 * (F + F.T)

    return F


# =========================================================
# Fiducial parameters
# =========================================================
theta_fid = np.array([
    proxy_mu0,
    proxy_muz,
    proxy_mulog10m,
    proxy_sigma0,
    proxy_sigmaz,
    proxy_sigmalog10m
])


# =========================================================
# Relative step sizes
# =========================================================
eps_rel = np.array([
    1e-2,
    1e-2,
    1e-2,
    1e-2,
    1e-2,
    1e-2
])


# =========================================================
# Compute Fisher + covariance
# =========================================================
F = fisher_matrix(theta_fid, eps_rel)

cov = np.linalg.inv(F)
errors = np.sqrt(np.diag(cov))

print("Parameter errors:")
print(errors)


# =========================================================
# Sanity checks
# =========================================================

print("\nSymmetric:", np.allclose(F, F.T))
print("Eigenvalues:", np.linalg.eigvals(F))


# =========================================================
# Step-size stability test
# =========================================================
print("\nStep size stability test:")

for factor in [0.5, 1.0, 2.0]:
    F_test = fisher_matrix(theta_fid, eps_rel * factor)
    cov_test = np.linalg.inv(F_test)
    err_test = np.sqrt(np.diag(cov_test))
    print(f"factor={factor} ->", err_test)