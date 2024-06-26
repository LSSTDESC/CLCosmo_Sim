from scipy.stats import norm
import numpy as np

def proxy_mu_f(logm, z, theta_rm):
    r"""proxy mu"""
    log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = theta_rm
    proxy_mu = proxy_mu0 + proxy_muz * np.log((1+z)/(1 + z0)) + proxy_mulog10m * (logm-log10m0)
    return proxy_mu

def proxy_sigma_f(logm, z, theta_rm):
    r"""proxy sigma"""
    log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = theta_rm
    proxy_sigma = proxy_sigma0 + proxy_sigmaz * np.log((1+z)/(1 + z0)) + proxy_sigmalog10m * (logm-log10m0)
    return proxy_sigma

def pdf_richness_mass_relation(richness, logm, z, theta_rm):
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
        richness-mass relation P(lambda|m,z)
    """
    log10m0, z0, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m = theta_rm
    proxy_mu = proxy_mu_f(logm, z, theta_rm)
    proxy_sigma = proxy_sigma_f(logm, z, theta_rm)
    return (1/richness)*np.exp(-(np.log(richness)-proxy_mu)**2/(2*proxy_sigma**2))/np.sqrt(2*np.pi*proxy_sigma**2)

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
    
