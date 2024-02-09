import numpy as np
import pickle
from statistics import median

def binning(corner): return [[corner[i],corner[i+1]] for i in range(len(corner)-1)]

def load(filename, **kwargs):
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)

redMaPPer_clusters = load('/pbs/throng/lsst/users/cpayerne/ThesisAtCCin2p3/Galaxy_Cluster_Catalogs_details/cosmoDC2/RedMapper_galaxy_clusters.pkl')

# z_corner = np.array([0.2, 0.35, 0.5, 0.65, 0.8 ])
# Z_bin = binning(z_corner)
# rich_corner = np.linspace(20, 120, 5)
# Obs_bin = binning(rich_corner)

#z_corner = np.array([0.2, 0.35, 0.5, 0.65, 0.8 ])
#z_corner = np.linspace(0.2, 1, 5)
#z_corner = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
#z_corner = np.array([0.2, 0.4, 0.6, 0.8, 1.])
z_corner = np.array([0.2, 0.4, 0.6, 0.8, 1.])
z_corner = np.array([0.2, 0.35, 0.50, 0.65, 0.8, 1])
z_corner = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1])
#z_corner = np.array([0.2, 0.35, 0.5, 0.65, 0.8, 1])
Z_bin = binning(z_corner)
rich_corner = np.array([20, 35, 55, 80, 120])
rich_corner = np.array([20, 30, 50, 80, 120])
rich_corner = np.array([20, 30, 50, 100, 200])
rich_corner = np.array([20, 30, 45, 60, 100, 200])
rich_corner = np.array([20, 35, 55, 80, 111, 200])
rich_corner = np.array([20, 35, 70, 100, 200])
#rich_corner = np.logspace(np.log10(20), np.log10(120), 5)
rich_corner = np.array([r for r in rich_corner])
Obs_bin = binning(rich_corner)


mask_z = (redMaPPer_clusters['redshift'] > z_corner[0])*(redMaPPer_clusters['redshift'] < z_corner[-1])
mask_obs = (redMaPPer_clusters['richness'] > rich_corner[0])*(redMaPPer_clusters['richness'] < rich_corner[-1])
z0 = median(redMaPPer_clusters['redshift'][mask_z*mask_obs])
richness0 = median(redMaPPer_clusters['richness'][mask_z*mask_obs])
