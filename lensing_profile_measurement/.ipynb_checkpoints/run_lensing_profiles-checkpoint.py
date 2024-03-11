import numpy as np
import clmm
import sys
import compute_lensing_profile_utils
import glob
import time
import pickle
from astropy.table import QTable, Table, vstack, join, hstack

def save_pickle(dat, filename, **kwargs):
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()

def mask(table, z_cl):
    masks = (table['mag_i'] < 24.5)*(table['mag_r'] < 28)*(table['z'] > z_cl + .1)
    return table[masks]

cosmo = clmm.Cosmology(H0 = 71.0, Omega_dm0 = 0.265 - 0.0448, Omega_b0 = 0.0448, Omega_k0 = 0.0)
name_lens_cat = '/pbs/throng/lsst/users/cpayerne/CLMassDC2/data/lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl'
#name_lens_cat = '/pbs/throng/lsst/users/cpayerne/CLMassDC2/data/lens_catalog_SkySim5000.pkl'
ra_name, dec_name, z_name = 'ra', 'dec', 'redshift'
obs_name = 'richness'
lens_cat = edit.load_pickle(name_lens_cat)
lens_cat_to_extract = lens_cat[(lens_cat['richness'] > 20)*(lens_cat['redshift'] > .2)]

where_bckd_catalog = '/sps/lsst/users/cpayerne/CLMassDC2/cosmoDC2/redmapper_clusters_new/'
file = glob.glob(where_bckd_catalog + 'l*')
cluster_id_saved = []
for f in file:
    cluster_id_saved.append(f.split('.pkl')[0].split('halo_')[1])
mask_id = np.isin(cluster_id_saved, lens_cat_to_extract['cluster_id'])
file = np.array(file)[mask_id]
print(len(file))
ind_profile = {n:[] for n in names}
for i, name_file in enumerate(file):
    if i%300==0: print(i)
    #cluster infos
    cluster_id = int(name_file.split('.pkl')[0].split('halo_')[1])
    lens = lens_cat[lens_cat['cluster_id'] == cluster_id][0]
    ra, dec, z = lens['ra'], lens['dec'], lens['redshift']
    obs = lens[obs_name]
    #bckgd galaxy catalog
    table = edit.load_pickle(name_file)
    table = mask(table, z)
    #add masks ?
    cluster_data = [cluster_id, ra, dec, z, obs]
    data_prf = compute_lensing_profile_utils.compute_lening_profile(cluster_id, ra, dec, z, label_pz, table, cosmo)
    data_to_save = cluster_data + data_prf
    for s, n in enumerate(names): ind_profile[n].append(data_to_save[s])
#break

save_pickle(Table(ind_profile), '/pbs/throng/lsst/users/cpayerne/CLMassDC2/data/data_new_version/ind_profile_redmapper.pkl')
