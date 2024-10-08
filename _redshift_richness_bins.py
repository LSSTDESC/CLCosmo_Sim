import numpy as np

def binning(corner): return [[corner[i],corner[i+1]] for i in range(len(corner)-1)]

z_corner = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1])
#z_corner = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1])
#z_corner = np.array([0.2, 0.4, 0.6, 0.8,  1])
Z_bin = np.array(binning(z_corner))
rich_corner = np.array([20, 35, 70, 100, 200])
#rich_corner = np.array([20, 30, 50, 80, 120, 200])
rich_corner = np.array([r for r in rich_corner])
Obs_bin = np.array(binning(rich_corner))
