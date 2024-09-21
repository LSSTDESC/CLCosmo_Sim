import sys
import numpy as np
sys.path.append('../../cluster-lensing-cov/')
import clens.util.constants as cn
from clens.util.parameters import CosmoParameters
from clens.util.scaling_relation import FiducialScalingRelation, Costanzi21ScalingRelation, Murata18ScalingRelation
from clens.util.survey import Survey
import clens.lensing.cov_DeltaSigma as cov_DeltaSigma_module

sys.path.append('../')
import _redshift_richness_bins

sys.path.append('../lensing_profile_measurement/')
import _config_lensing_profiles

code, index = sys.argv

    
class FiducialScalingRelation(object):
    def __init__(self, proxy_mu0, proxy_muz, proxy_mulog10m, proxy_sigma0, proxy_sigmaz, proxy_sigmalog10m ):
        self.log10M0 = np.log10(10**14.3)
        self.z0 = .5
        self.proxy_mu0 = proxy_mu0
        self.proxy_muz = proxy_muz
        self.proxy_mulog10m = proxy_mulog10m
        self.proxy_sigma0 = proxy_sigma0
        self.proxy_sigmaz = proxy_sigmaz
        self.proxy_sigmalog10m = proxy_sigmaz

    def lnlambda_lnM(self, lnM, z):
        M200m = np.exp(lnM)
        M = M200m/1.18
        log10M = np.log10(M)
        lnlambda = self.proxy_mu0 + self.proxy_muz * np.log((1+z)/(1 + self.z0)) + self.proxy_mulog10m * (log10M-self.log10M0)
        return lnlambda

    def scatter(self, lnM, z):
        M200m = np.exp(lnM)
        M = M200m/1.18
        log10M = np.log10(M)
        self.sigma_lambda = self.proxy_sigma0 + self.proxy_sigmaz * np.log((1+z)/(1 + self.z0)) + self.proxy_sigmalog10m * (log10M-self.log10M0)
        return self.sigma_lambda 

rp_min = _config_lensing_profiles.bin_edges[0]
rp_max = _config_lensing_profiles.bin_edges[-1]
n_rp = len(_config_lensing_profiles.bin_edges)-1

def demo_cov(z_bin, richness_bin):
    co = CosmoParameters()
    #sr = FiducialScalingRelation()
    #sr = Costanzi21ScalingRelation()
    #sr = Murata18ScalingRelation()
    sr=FiducialScalingRelation(3.345953364933381,
                                 0.06378674560099672,
                                 2.2274595352385975,
                                 0.5634900828247924,
                                 -0.04529427946062734,
                                 0.09764409931399064)
    su = Survey(zs_min=0.2, zs_max=3, top_hat=False, n_src_arcmin=27, sigma_gamma=0.3)
    fsky = 440 / 41253.
    cds = cov_DeltaSigma_module.CovDeltaSigma(co=co, su=su, sr=sr, fsky=fsky)
    cds.calc_cov(rp_min=rp_min, rp_max=rp_max, n_rp=n_rp, 
                 zh_min=z_bin[0], zh_max=z_bin[1], 
                 lambda_min=richness_bin[0], lambda_max=richness_bin[1], diag_only=False)
    return cds
Obs_bin = _redshift_richness_bins.Obs_bin
Z_bin = _redshift_richness_bins.Z_bin
demo_cov(Z_bin[index], Obs_bin[index])