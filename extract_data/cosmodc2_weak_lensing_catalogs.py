import numpy as np
import GCRCatalogs
import healpy
import glob
import clmm
import time
from astropy.table import QTable, Table, vstack, join, hstack
#from astropy.io.misc.hdf5 import write_table_hdf5, read_table_hdf5
import pickle,sys
sys.path.append('/pbs/throng/lsst/users/cpayerne/LikelihoodsClusterAbundance/modules/')
import edit
import extract_in_cosmodc2_utils as cosmodc2
import photoz_utils

from clmm import Cosmology
from scipy.integrate import simps

#cosmoDC2 cosmology
cosmo = Cosmology(H0 = 71.0, Omega_dm0 = 0.265 - 0.0448, Omega_b0 = 0.0448, Omega_k0 = 0.0)

#connection with qserv
import mysql
from mysql.connector import Error
start, end = int(sys.argv[1]), int(sys.argv[2])

#select galaxy clusters
lens_catalog_name='/pbs/throng/lsst/users/cpayerne/CLMassDC2/data/lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl'
#lens_catalog_name='/pbs/throng/lsst/users/cpayerne/CLMassDC2/data/lens_catalog_SkySim5000.pkl'

lens_catalog=edit.load_pickle(lens_catalog_name)

# select subsample of clusters SkySim5000
# mask = (lens_catalog['baseDC2/sod_halo_mass']/cosmo['h'] > 1e14)*(lens_catalog['redshift'] < 1.2)*(lens_catalog['redshift'] > .2)
# lens_catalog = lens_catalog[mask]
# lens_catalog['cluster_id'] = lens_catalog['halo_id']

where_to_save='/sps/lsst/users/cpayerne/CLMassDC2/cosmoDC2/redmapper_clusters_new/'
#where_to_save='/sps/lsst/users/cpayerne/CLMassDC2/cosmoDC2/dm_halos/'

#select subsample of clusters redMaPPer#
mask_select = (lens_catalog['richness'] > 20)*(lens_catalog['redshift'] > .2)
lens_catalog = lens_catalog[mask_select]
#mask_n=np.arange(start, end)
lens_catalog_truncated=lens_catalog#[mask_n]

#file_already_saved = glob.glob(where_to_save + 'l*')
#cluster_id_saved = []
#for f in file_already_saved:
#    cluster_id_saved.append(int(f.split('.pkl')[0].split('halo_')[1]))
#mask_saved = np.isin(lens_catalog_truncated['cluster_id'], cluster_id_saved)
#lens_catalog_truncated = lens_catalog_truncated[np.invert(mask_saved)]
lens_catalog_truncated = lens_catalog_truncated[np.arange(start, end)]
print(len(lens_catalog_truncated))
print('-----')
#print(min(lens_catalog_truncated['richness']))
#print(min(lens_catalog_truncated['redshift']))
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
    zmin = lens_z + .05
    theta_max = (rmax/lens_distance) * (180./np.pi)
    query = "SELECT data.coord_ra as ra, data.coord_dec as dec, data.redshift as z, "
    query += "data.galaxy_id as galaxy_id, "
    query += "data.mag_i, data.mag_r, data.mag_y, "
    query += "data.shear_1 as shear1, data.shear_2 as shear2, data.convergence as kappa, "
    query += "data.ellipticity_1_true as e1_true, data.ellipticity_2_true as e2_true " 
    query += "FROM cosmoDC2_v1_1_4_image.data as data "
    query += f"WHERE data.redshift >= {zmin} AND data.redshift < {zmax} "
    query += f"AND scisql_s2PtInCircle(coord_ra, coord_dec, {ra}, {dec}, {theta_max}) = 1 "
    query += f"AND data.mag_i <= 24.6 "
    query += f"AND data.mag_r <= 28.0 "
    query += ";" 
    return query

def query_photoz():
    
    return ['photoz_pdf', 'photoz_mean','photoz_mode','photoz_odds','galaxy_id']

for n, lens in enumerate(lens_catalog_truncated):
    
    if len(lens_catalog_truncated) == 0: continue

    #print('halo index in list' + str(mask_n[n]))
    #cluster metadata
    z, ra, dec=lens['redshift'], lens['ra'], lens['dec']
    cluster_id=lens['cluster_id']
    lens_distance=cosmo.eval_da(z)
    name_cat = 'lensing_catalog_halo_' + str(cluster_id)
    name_full_cat = where_to_save + name_cat + '.pkl'
    if name_full_cat in glob.glob(where_to_save + '*'):
        print('already saved')
        continue
    conn = mysql.connector.connect(host='ccqserv201', user='qsmaster', port=30040)
    cursor = conn.cursor(dictionary=True, buffered=True)

    #extract background galaxies with qserv (only true photoz)
    print('==> extracting true redshift + shapes infos (Qserv)')
    dat_extract=cosmodc2.extract(qserv_query = qserv_query(z, lens_distance, ra, dec, rmax = 10.1),
                                        conn_qserv=conn, cosmo=cosmo)
    conn.close()
    print('Qserv: Number of background galaxies = ' + str(len(dat_extract['galaxy_id'])))
    if photoz == True:
    
        #extract photometric redshifts with GCRCatalogs
        print('===> extracting photoz redshift infos (GCRCatalogs)')
        id_gal=dat_extract['galaxy_id']
        ras=dat_extract['ra']
        decs=dat_extract['dec']
        
        #find all different healpix pixels
        healpix = np.unique(healpy.ang2pix(32, ras, decs, nest=False, lonlat=True))
        healpix = healpix[np.isin(healpix, healpix_dc2)]
        print(healpix)
        table_photoz = Table()
        table_photoz['galaxy_id'] = id_gal
       
        timei = time.time()
        for k, gc_ in enumerate(photoz_gc):
            
            pz_table = Table(names = ['sigmac_photoz', 'p_background', 'photoz_dispersion', 
                                  'sigmac_estimate_0', 'sigmac_estimate_1', 'sigmac_estimate_2', 
                                  'z_estimate_0', 'z_estimate_1', 'z_estimate_2', 
                                  'galaxy_id', 'photoz_mean', 'photoz_mode', 'photoz_odds'])
            
            photoz_gc_ = GCRCatalogs.load_catalog(gc_)
            for i, hp in enumerate(healpix):
                #browse healpix pixels
                print('heapix pixel = ' + str(hp))
                chunk = photoz_gc_.get_quantities(query_photoz(), native_filters=[f'healpix_pixel=={hp}'], return_iterator=True)
                
                for j in range(3):
                    
                    #browse chunk data
                    print('chunk = ' + str(j))
                    dat_extract_photoz_chunk = Table(next(chunk))
                    print('full healpix = ' + str(len(dat_extract_photoz_chunk['galaxy_id'])))
                    
                    #use only selected galaxies
                    dat_extract_photoz_chunk_truncated = dat_extract_photoz_chunk[np.isin(dat_extract_photoz_chunk['galaxy_id'],
                                                                                          table_photoz['galaxy_id'])]
                    if len(dat_extract_photoz_chunk_truncated['galaxy_id']) == 0: continue
                    
                    print('truncated healpix = ' + str(len(dat_extract_photoz_chunk_truncated['galaxy_id'])))
                    pzbins_table=np.array([z_bins for i in range(len(dat_extract_photoz_chunk_truncated['photoz_pdf'].data))])
                    
                    #compute WL weights with 
                    pz_quantities_chunk = photoz_utils.compute_photoz_quantities(z, dat_extract_photoz_chunk_truncated['photoz_pdf'], 
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
    edit.save_pickle(dat_extract, where_to_save + name_cat + '.pkl')
    #break
