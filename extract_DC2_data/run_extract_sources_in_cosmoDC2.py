import numpy as np
import GCRCatalogs
GCRCatalogs.set_root_dir_by_site('in2p3')
import healpy
import glob, sys
import clmm
import time
import _utils_photometric_redshifts
import _utils_extract_sources_in_cosmoDC2
import _config_extract_sources_in_cosmoDC2 

sys.path.append('../lensing_profile_measurement')
import _config_lensing_profiles
import _utils_lensing_profiles
#test

from astropy.table import QTable, Table, vstack, join, hstack
import pickle,sys
def load(filename, **kwargs):
    """Loads GalaxyCluster object to filename using Pickle"""
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)
def save_pickle(dat, filename, **kwargs):
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()

from clmm import Cosmology
from scipy.integrate import simps

#cosmoDC2 cosmology
cosmo = Cosmology(H0 = 71.0, Omega_dm0 = 0.265 - 0.0448, Omega_b0 = 0.0448, Omega_k0 = 0.0)

#connection with qserv
import mysql
from mysql.connector import Error
start, end = int(sys.argv[1]), int(sys.argv[2])

#select galaxy clusters
lens_catalog_name=_config_extract_sources_in_cosmoDC2.lens_catalog_name
lens_catalog=np.load(lens_catalog_name, allow_pickle=True)
where_to_save=_config_extract_sources_in_cosmoDC2.where_to_save
mask_select = (lens_catalog['richness'] > 20)*(lens_catalog['redshift'] > .2)
lens_catalog = lens_catalog[mask_select]
lens_catalog_truncated=lens_catalog
lens_catalog_truncated = lens_catalog_truncated[np.arange(start, end)]
#load source catalogs
photoz=True
if photoz == True:
    gc_bpz  = "cosmoDC2_v1.1.4_image_with_photozs_v1"
    gc_flex = "cosmoDC2_v1.1.4_image_with_photozs_flexzboost_v1"
    #list of healpix in cosmoDC2
    healpix_dc2 = GCRCatalogs.load_catalog("cosmoDC2_v1.1.4_image").get_catalog_info()['healpix_pixels']
    z_bins  = GCRCatalogs.load_catalog(gc_flex).photoz_pdf_bin_centers
    z_bins[0] = 1e-7
    photoz_gc=[gc_bpz, gc_flex]
    photoz_label=['_bpz', '_flex']

def qserv_query(lens_z, lens_distance, ra, dec, rmax = 10):
    r"""
    quantities wanted + cuts for qserv
    Attributes:
    -----------
    z: float
        lens redshift
    lens_distance: float
        distance to the cluster
    ra: float
        lens right ascension
    dec: float
        lens declinaison
    rmax: float
        maximum radius
    """
    zmax = 3.
    zmin = 0.#zmin = lens_z + .05
    theta_max = (rmax/lens_distance) * (180./np.pi)
    query = "SELECT data.coord_ra as ra, data.coord_dec as dec, data.redshift as z, "
    query += "data.galaxy_id as galaxy_id, "
    query += "data.mag_i, data.mag_r, data.mag_y, "
    query += "data.shear_1 as shear1, data.shear_2 as shear2, data.convergence as kappa, "
    query += "data.ellipticity_1_true as e1_true_uncorr, data.ellipticity_2_true as e2_true_uncorr " 
    query += "FROM cosmoDC2_v1_1_4_image.data as data "
    query += f"WHERE data.redshift >= {zmin} AND data.redshift < {zmax} "
    query += f"AND scisql_s2PtInCircle(coord_ra, coord_dec, {ra}, {dec}, {theta_max}) = 1 "
    query += f"AND data.mag_i <= 24.6 "
    query += f"AND data.mag_r <= 28.0 "
    query += ";" 
    return query

def query_photoz():
    
    return ['photoz_pdf', 'photoz_mean','photoz_mode','photoz_odds','galaxy_id']

if _config_extract_sources_in_cosmoDC2.compute_individual_lensing_profile:
    ind_profile = {n:[] for n in _config_lensing_profiles.names}

for n, lens in enumerate(lens_catalog_truncated):
    
    if len(lens_catalog_truncated) == 0: continue

    z, ra, dec=lens['redshift'], lens['ra'], lens['dec']
    cluster_id=lens['cluster_id']
    lens_distance=cosmo.eval_da(z)
    print(f'====> cluster n°{cluster_id}, at redshift z = {z:.2f} and distance DA = {lens_distance:.2f} Mpc')
    name_cat = 'lensing_catalog_halo_' + str(cluster_id)
    name_full_cat = where_to_save + name_cat + '.pkl'
    print(f'====> name for save: '+name_full_cat)
    if _config_extract_sources_in_cosmoDC2.save_catalog:
        if name_full_cat in glob.glob(where_to_save + '*'):
            print('already saved')
        continue
    conn = mysql.connector.connect(host='ccqserv201', user='qsmaster', port=30040)
    cursor = conn.cursor(dictionary=True, buffered=True)
    #extract background galaxies with qserv (only true photoz)
    print('==> extracting true redshift + shapes infos (Qserv)')
    dat_extract=_utils_extract_sources_in_cosmoDC2.extract(qserv_query = qserv_query(z, lens_distance, ra, dec, rmax = 10.1),
                                        conn_qserv=conn, cosmo=cosmo)
    conn.close()
    print('==> Qserv: Number of extracted galaxies = ' + str(len(dat_extract['galaxy_id'])))
    #extract photometric redshifts with GCRCatalogs
    print('==> extracting photoz redshift infos (GCRCatalogs)')
    id_gal=dat_extract['galaxy_id']
    ras=dat_extract['ra']
    decs=dat_extract['dec']
    #find all different healpix pixels
    healpix = np.unique(healpy.ang2pix(32, ras, decs, nest=False, lonlat=True))
    healpix = healpix[np.isin(healpix, healpix_dc2)]
    print(f'==> healpix pixels covered (nside: 32): {healpix}')
    table_photoz = Table()
    table_photoz['galaxy_id'] = id_gal
    if photoz == True:
        timei = time.time()
        for k, gc_ in enumerate(photoz_gc):
            print(f'----> extraction in {gc_bpz}')
            pz_table = Table(names = ['sigmac_photoz', 'p_background', 'photoz_dispersion', 
                                  'sigmac_estimate_0', 'sigmac_estimate_1', 'sigmac_estimate_2', 
                                  'z_estimate_0', 'z_estimate_1', 'z_estimate_2', 
                                  'galaxy_id', 'photoz_mean', 'photoz_mode', 'photoz_odds'])
            photoz_gc_ = GCRCatalogs.load_catalog(gc_)
            for i, hp in enumerate(healpix):
                #browse healpix pixels
                print(f'-----> heapix pixel = ' + str(hp))
                chunk = photoz_gc_.get_quantities(query_photoz(), native_filters=[f'healpix_pixel=={hp}'], return_iterator=True)
                for j in range(3):
                    #browse chunk data
                    print('chunk = ' + str(j))
                    try: 
                        dat_extract_photoz_chunk = Table(next(chunk))
                    except: continue
                    print(f'number of galaxies in the full healpix = ' + str(len(dat_extract_photoz_chunk['galaxy_id'])))
                    #use only selected galaxies
                    dat_extract_photoz_chunk_truncated = dat_extract_photoz_chunk[np.isin(dat_extract_photoz_chunk['galaxy_id'],
                                                                                          table_photoz['galaxy_id'])]
                    if len(dat_extract_photoz_chunk_truncated['galaxy_id']) == 0: continue
                    print('number of galaxies tageted in the healpix = ' + str(len(dat_extract_photoz_chunk_truncated['galaxy_id'])))
                    pzbins_table=np.array([z_bins for i in range(len(dat_extract_photoz_chunk_truncated['photoz_pdf'].data))])
                    #compute WL weights with 
                    pz_quantities_chunk = _utils_photometric_redshifts.compute_photoz_quantities(z, dat_extract_photoz_chunk_truncated['photoz_pdf'], 
                                                                           pzbins_table, n_samples_per_pdf=3, cosmo=cosmo,
                                                                           use_clmm=False)
                    pz_quantities_chunk['galaxy_id'] = dat_extract_photoz_chunk_truncated['galaxy_id']
                    pz_quantities_chunk['photoz_mean'] = dat_extract_photoz_chunk_truncated['photoz_mean']
                    pz_quantities_chunk['photoz_mode'] = dat_extract_photoz_chunk_truncated['photoz_mode']
                    pz_quantities_chunk['photoz_odds'] = dat_extract_photoz_chunk_truncated['photoz_odds']
                    pz_table = vstack([pz_table, pz_quantities_chunk])

            #rename photoz columns
            colnames = pz_table.colnames
            for name in colnames:
                if name!='galaxy_id':
                    pz_table.rename_column(name, name + photoz_label[k])
            print(len(pz_table['galaxy_id']))
            dat_extract = join(Table(dat_extract), pz_table, keys='galaxy_id')

        timef = time.time()
        print(str(timef-timei) + ' (s)')
    
    if _config_extract_sources_in_cosmoDC2.save_catalog:
        save_pickle(dat_extract, where_to_save + name_cat + '.pkl')
    
    if _config_extract_sources_in_cosmoDC2.compute_individual_lensing_profile:
        cluster_data = [cluster_id, ra, dec, z, lens['richness']]
        names = _config_lensing_profiles.names
        #do magnitude selection (magnitude)
        dat_extract_mag_cut = _config_lensing_profiles.source_selection_magnitude(dat_extract)
        #do bckd source selection (photometric redshift)
        data_to_save = cluster_data
        for label_ in _config_lensing_profiles.label_pz:
            if label_ != 'true':
                mask_z = dat_extract_mag_cut['p_background'+'_'+label_] > _config_lensing_profiles.p_background_min
                mask_z *= dat_extract_mag_cut['photoz_mean'+'_'+label_] > z + _config_lensing_profiles.Dz_background
            else: 
                mask_z = dat_extract_mag_cut['z'] > z + 0.1
            dat_extract_mag_cut_bkgd_pz = dat_extract_mag_cut[mask_z]
            bin_edges = _config_lensing_profiles.bin_edges
            dat_prf_pz = _utils_lensing_profiles.compute_lensing_profile(cluster_id, ra, dec, z, bin_edges, label_, dat_extract_mag_cut_bkgd_pz, cosmo)
            data_to_save = data_to_save + dat_prf_pz
        for s, n in enumerate(names): ind_profile[n].append(data_to_save[s])
                
path = '../data/ind_profile_redmapper_per_cluster_index/'
save_pickle(Table(ind_profile), path+f'ind_profile_redmapper_{start}_{end}.pkl')
