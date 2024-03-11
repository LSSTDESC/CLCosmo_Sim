import numpy as np
import glob
import pickle
from astropy.table import QTable, Table, vstack, join
def save_pickle(dat, filename, **kwargs):
     file = open(filename,'wb')
     pickle.dump(dat, file)
     file.close()
files = glob.glob('./ind_profile_redmapper_per_cluster_index/ind_profile_*')
for i, f in enumerate(files):
    data_f = np.load(f, allow_pickle=True)
    if i==0:
        data_f0 = data_f
        continue
    else:
        data_f0 = vstack([data_f0, data_f])
save_pickle(data_f0, 'ind_profile_redmapper.pkl')