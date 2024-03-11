import numpy as np
import pickle
from statistics import median

def binning(corner): return [[corner[i],corner[i+1]] for i in range(len(corner)-1)]

z_corner = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1])
Z_bin = binning(z_corner)
rich_corner = np.array([20, 35, 70, 100, 200])
rich_corner = np.array([r for r in rich_corner])
Obs_bin = binning(rich_corner)
