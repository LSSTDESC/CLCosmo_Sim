from scipy.stats import norm
import numpy as np

def purity(richness, z, theta_purity):
    r"""
    Attributes:
    -----------
    richness : array
        cluster richness
    z : float
        cluster redshift
    theta_purity: array
        parameters of purity
    Returns:
    --------
    purity : array
        purity of cluster detection
    """
    richness0, p0, p1, np_ = theta_purity
    richness_scale = np.log(richness0) + p0 + p1*(1 + z)
    richness_rescale = richness/richness_scale
    return richness_rescale**np_/(richness_rescale**np_+1)

