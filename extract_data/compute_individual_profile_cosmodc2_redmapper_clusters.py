import numpy as np
import clmm
import sys
sys.path.append('/pbs/throng/lsst/users/cpayerne/LikelihoodsClusterAbundance/modules/')
import edit
import glob
import time
from astropy.table import QTable, Table, vstack, join, hstack

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

bin_edges = clmm.dataops.make_bins(0.5, 10, 15, method='evenlog10width')

label_pz = ['true', 'flex', 'bpz']
label_prf = ['DSt', 'DSx', 'W_l', 'radius']
names_cl=['id', ra_name, dec_name, z_name, obs_name]
label_prf_full = [label_prf_ + '_' + label_pz_ for label_pz_ in label_pz for label_prf_ in label_prf]
names = names_cl + label_prf_full
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
    cl = clmm.galaxycluster.GalaxyCluster('halo', ra, dec, z, clmm.gcdata.GCData(Table(table)))
    theta1, g_t, g_x = cl.compute_tangential_and_cross_components(is_deltasigma=False, cosmo=cosmo)
    cluster_data = [cluster_id, ra, dec, z, obs]
    data_prf = []
    for label in label_pz:
        
        if label=='true': sigma_c = cosmo.eval_sigma_crit(z, cl.galcat['z'])
        elif label=='flex': sigma_c = cl.galcat['sigmac_photoz_flex']
        elif label=='bpz': sigma_c = cl.galcat['sigmac_photoz_bpz']

        cl.galcat['dst'] = sigma_c*cl.galcat['et']
        cl.galcat['dsx'] = sigma_c*cl.galcat['ex']
        cl.galcat['w_ls'] = sigma_c**(-2.)
        ce = clmm.ClusterEnsemble('id', [])

        p = ce.make_individual_radial_profile(cl, 'Mpc', bins=bin_edges, error_model='ste',
                                           cosmo=cosmo, tan_component_in='dst', cross_component_in='dsx',
                                           tan_component_out='gt', cross_component_out='gx',
                                           tan_component_in_err=None, cross_component_in_err=None,
                                           weights_in='w_ls', weights_out='W_l')
        data = ce.data[0]
        data_prf.extend([data['gt'], data['gx'], data['W_l'], data['radius']])
    data_to_save = cluster_data + data_prf
    for s, n in enumerate(names): ind_profile[n].append(data_to_save[s])
#break

edit.save_pickle(Table(ind_profile), '/pbs/throng/lsst/users/cpayerne/CLMassDC2/data/data_new_version/ind_profile_redmapper.pkl')
