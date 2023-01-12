from scipy.stats import norm
import numpy as np
import pyccl as ccl

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

# class Richness():
#     r"""
#     class for the modelling of galaxy cluster richness (arXiv:1904.07524)
#     """
#     def __init__(self, theta_mu, theta_sigma, theta_pivot):
#         r"""
#         Attributes:
#         -----------
#         theta_mu : 3D array
#             paramters of the mean richness-mass relation
#             mu_0, A_z^mu, A_m^mu
#         theta_sigma : 3D array
#             paramters of the mean richness-mass relation
#             sigma_0, sigma_z^mu, sigma_m^mu
#         theta_pivot : 3D array
#             paramters of the mean richness-mass relation
#             m_0, z_0
#         """
#         self.theta_mu = theta_mu
#         self.theta_sigma = theta_sigma
#         self.theta_pivot = theta_pivot

#     def proxy_mu(self, redshift, logm):
#         r"""
#         Attributes:
#         -----------
#         redshift: array
#             redshift of the cluster
#         logm: array
#             \log10(M) of the cluster
#         Returns:
#         --------
#         mu: array
#             mean richness at mass logm and redshift
#         """
#         log10m0, z0 = self.theta_pivot
#         proxy_mu0, proxy_muz, proxy_mulog10m = self.theta_mu
#         return proxy_mu0 + proxy_muz * np.log((1+redshift)/(1 + z0)) + proxy_mulog10m * (logm-log10m0)

#     def proxy_sigma(self, redshift, logm):
#         r"""
#         Attributes:
#         -----------
#         redshift: array
#             redshift of the cluster
#         logm: array
#             \log10(M) of the cluster
#         Returns:
#         --------
#         sigma: array
#             standard deviation of richness at mass logm and redshift
#         """
#         log10m0, z0 = self.theta_pivot
#         proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = self.theta_sigma
#         return proxy_sigma0 + proxy_sigmaz * np.log((1+redshift)/(1 + z0)) + proxy_sigmalog10m * (logm-log10m0)

#     def pdf(self,proxy, redshift, logm):
#         r"""
#         Attributes:
#         -----------
#         lnlambda : array
#             richness
#         redshift: array
#             redshift of the cluster
#         logm: array
#             \log10(M) of the cluster
#         Returns:
#         --------
#         P: array
#             probability density to find richness lnlambda given mass and redhsift
#         """
#         proxy_mu = self.proxy_mu(redshift, logm)
#         proxy_sigma = self.proxy_sigma(redshift, logm)
#         return np.exp(-.5*(proxy-proxy_mu)**2/proxy_sigma**2)/np.sqrt(2*np.pi*proxy_sigma**2)

#     def cdf(self, proxy_bin, redshift, logm):
#         r"""
#         Attributes:
#         -----------
#         lnlambda_bin : array
#             richness range
#         redshift: array
#             redshift of the cluster
#         logm: array
#             \log10(M) of the cluster
#         Returns:
#         --------
#         P: array
#             cumulative in richness bin given mass and redhsift
#         """
#         proxy_mu = self.proxy_mu(redshift, logm)
#         proxy_sigma = self.proxy_sigma(redshift, logm)
#         return norm.cdf(proxy_bin[1], proxy_mu, proxy_sigma) - norm.cdf(proxy_bin[0], proxy_mu, proxy_sigma)


#     def lnLambda_random(self, redshift, logm):

#         random = np.random.randn(len(logm))

#         mu = self.mu_loglambda_logM_f(redshift, logm)
#         sigma = self.sigma_loglambda_logm_f(redshift, logm)
#         return mu + sigma * random

