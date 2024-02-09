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
    hmf = hmd.__call__(cosmo, 10**np.array(logm), 1./(1. + z))
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



