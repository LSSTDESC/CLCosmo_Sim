import pyccl as ccl
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
from scipy.integrate import quad,simps, dblquad
from astropy import units as u
import itertools
from astropy import constants as const
from scipy import interpolate
import emcee
import iminuit
from iminuit import Minuit
from scipy.special import gamma
import pickle as pkl
from astropy.cosmology import FlatLambdaCDM

class Richness():
    r"""
    class for the modelling of galaxy cluster richness (arXiv:1904.07524)
    """
    def __init__(self, theta_mu, theta_sigma, theta_pivot):
        r"""
        Attributes:
        -----------
        theta_mu : 3D array
            paramters of the mean richness-mass relation
            mu_0, A_z^mu, A_m^mu
        theta_sigma : 3D array
            paramters of the mean richness-mass relation
            sigma_0, sigma_z^mu, sigma_m^mu
        theta_pivot : 3D array
            paramters of the mean richness-mass relation
            m_0, z_0
        """
        self.theta_mu = theta_mu
        self.theta_sigma = theta_sigma
        self.theta_pivot = theta_pivot

    def mu_loglambda_logM_f(self, redshift, logm):
        r"""
        Attributes:
        -----------
        redshift: array
            redshift of the cluster
        logm: array
            \log10(M) of the cluster
        Returns:
        --------
        mu: array
            mean richness at mass logm and redshift
        """
        m0, z0 = self.theta_pivot
        loglambda0, A_z_mu, A_logm_mu = self.theta_mu
        sigma_lambda0, A_z_sigma, A_logm_sigma = self.theta_sigma
        return loglambda0 + A_z_mu * np.log10((1+redshift)/(1 + z0)) + A_logm_mu * (logm-np.log10(m0))

    def sigma_loglambda_logm_f(self, redshift, logm):
        r"""
        Attributes:
        -----------
        redshift: array
            redshift of the cluster
        logm: array
            \log10(M) of the cluster
        Returns:
        --------
        sigma: array
            standard deviation of richness at mass logm and redshift
        """
        m0, z0 = self.theta_pivot
        loglambda0, A_z_mu, A_logm_mu = self.theta_mu
        sigma_lambda0, A_z_sigma, A_logm_sigma = self.theta_sigma
        return sigma_lambda0 + A_z_sigma * np.log10((1+redshift)/(1 + z0)) + A_logm_sigma * (logm-np.log10(m0))
    
    def pdf(self,lnLambda, redshift, logm):
        r"""
        Attributes:
        -----------
        lnlambda : array
            richness
        redshift: array
            redshift of the cluster
        logm: array
            \log10(M) of the cluster
        Returns:
        --------
        P: array
            probability density to find richness lnlambda given mass and redhsift
        """
        mu = self.mu_loglambda_logM_f(redshift, logm)
        sigma = self.sigma_loglambda_logm_f(redshift, logm)
        return np.exp(-.5*(lnLambda-mu)**2/sigma**2)/np.sqrt(2*np.pi*sigma**2)

    def cdf(self, lnlambda_bin, redshift, logm):
        r"""
        Attributes:
        -----------
        lnlambda_bin : array
            richness range
        redshift: array
            redshift of the cluster
        logm: array
            \log10(M) of the cluster
        Returns:
        --------
        P: array
            cumulative in richness bin given mass and redhsift
        """
        mu = self.mu_loglambda_logM_f(redshift, logm)
        sigma = self.sigma_loglambda_logm_f(redshift, logm)
        return norm.cdf(lnlambda_bin[1], mu, sigma) - norm.cdf(lnlambda_bin[0], mu, sigma)

    r"""
    def lnLambda_random(self, redshift, logm):

        random = np.random.randn(len(logm))

        mu = self.mu_loglambda_logM_f(redshift, logm)
        sigma = self.sigma_loglambda_logm_f(redshift, logm)
        return mu + sigma * random
    """
