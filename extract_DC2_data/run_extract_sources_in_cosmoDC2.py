import numpy as np
import GCRCatalogs
GCRCatalogs.set_root_dir_by_site('in2p3')
import healpy
import glob, sys
import astropy.units as u
from astropy.io import fits as fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy
import clmm
import time
import _utils_extract_sources_in_cosmoDC2 
import _config_extract_sources_cosmoDC2
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
import argparse

path_to_data = '../../CLCosmo_Sim_database/data/'


def collect_argparser():
    parser = argparse.ArgumentParser(description="Let's extract source cat")
    parser.add_argument("--index_split", type=int, required=True)
    parser.add_argument("--n_splits", type=int, required=True)
    parser.add_argument("--lens_catalog_name", type=str, required=False, default=path_to_data+'lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl')
    parser.add_argument("--compute_individual_lensing_profile", type=str, required=True)
    parser.add_argument("--save_catalog", type=str, required=False, default='False')
    parser.add_argument("--where_to_save_catalog", type=str, required=False, default='/sps/lsst/users/cpayerne/CLMassDC2/cosmoDC2/redmapper_clusters_new/')
    return parser.parse_args()

#select galaxy clusters
_config_extract_sources_in_cosmoDC2 = collect_argparser()
lens_catalog_name=_config_extract_sources_in_cosmoDC2.lens_catalog_name
lens_catalog=np.load(lens_catalog_name, allow_pickle=True)
mask_select=(lens_catalog['richness'] > 20)*(lens_catalog['redshift'] > .2)
lens_catalog=lens_catalog[mask_select]
lens_catalog_truncated=lens_catalog
use_only_member_galaxies = True

n_cl=len(lens_catalog_truncated)
index_cl=np.arange(n_cl)
split_lists=np.array_split(index_cl, _config_extract_sources_in_cosmoDC2.n_splits)
lens_catalog_truncated=lens_catalog_truncated[np.array(split_lists[_config_extract_sources_in_cosmoDC2.index_split])]
n_cl_to_extract=len(lens_catalog_truncated)

print(len(lens_catalog_truncated))
if lens_catalog_name==path_to_data + 'random_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl': suff_='_random'
elif lens_catalog_name==path_to_data + 'lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl': suff_=''
path_where_to_save_profiles = path_to_data + f'ind_profile_redmapper_per_cluster'+suff_+'_index/'
name_save = path_where_to_save_profiles+f'ind_profile_redmapper'+suff_+f'_split={_config_extract_sources_in_cosmoDC2.index_split}_nsplits={_config_extract_sources_in_cosmoDC2.n_splits}_ncl={n_cl_to_extract}.pkl'
print(name_save)
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
    query += "data.galaxy_id as galaxy_id, data.halo_id as halo_id, "
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

if _config_extract_sources_in_cosmoDC2.compute_individual_lensing_profile=='True':
    ind_profile = {n:[] for n in _config_lensing_profiles.names}

for n, lens in enumerate(lens_catalog_truncated):
    
    if len(lens_catalog_truncated) == 0: continue

    z, ra, dec=lens['redshift'], lens['ra'], lens['dec']
    cluster_id=lens['cluster_id']
    #cluster member galaxies
    id_member_galaxy = lens_catalog['id_member'][n]
    ra_member_galaxy = lens_catalog['ra_member'][n]
    dec_member_galaxy = lens_catalog['dec_member'][n]
    
    #ra, dec = np.mean(ra_member_galaxy), np.mean(dec_member_galaxy)

    lens_distance=cosmo.eval_da(z)
    print(f'====> cluster n°{cluster_id}, at redshift z = {z:.2f} and distance DA = {lens_distance:.2f} Mpc')
    print(f'====> ra = {ra:.2f} deg and dec = {dec:.2f} deg')
    l=lens['richness']
    print(f'====> richness = {l:.2f}')
    print(f'====> ra mean = {np.mean(ra_member_galaxy):.2f} deg and dec mean = {np.mean(dec_member_galaxy):.2f} deg')
    conn = mysql.connector.connect(host='ccqserv201', user='qsmaster', port=30040)
    cursor = conn.cursor(dictionary=True, buffered=True)
    #extract background galaxies with qserv (only true photoz)
    print('==> extracting true redshift + shapes infos (Qserv)')
    dat_extract=_utils_extract_sources_in_cosmoDC2.extract(qserv_query = qserv_query(z, lens_distance, ra, dec, rmax = 2),
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
            pz_table=_utils_extract_sources_in_cosmoDC2.extract_photoz(z, pz_catalog = gc_, z_bins_pz = z_bins, healpix_list = healpix, 
                                                               id_gal_to_extract = id_gal, _query_photoz = query_photoz(), cosmo = cosmo)
            #rename photoz columns
            colnames = pz_table.colnames
            for name in colnames:
                if name!='galaxy_id':
                    pz_table.rename_column(name, name + photoz_label[k])
            print(len(pz_table['galaxy_id']))
            dat_extract = join(Table(dat_extract), pz_table, keys='galaxy_id')

        timef = time.time()
        print(str(timef-timei) + ' (s)')
    
    if _config_extract_sources_in_cosmoDC2.compute_individual_lensing_profile=='True':
        cluster_data = [cluster_id, ra, dec, z, lens['richness']]
        names = _config_lensing_profiles.names
        dat_extract_mag_cut = _config_extract_sources_cosmoDC2.source_selection_magnitude(dat_extract)
        data_to_save = cluster_data
        for label_ in _config_lensing_profiles.label_pz:
            if label_ != 'true':
                mask_z = dat_extract_mag_cut['p_background'+'_'+label_] > _config_extract_sources_cosmoDC2.p_background_min
                mask_z = dat_extract_mag_cut['photoz_mean'+'_'+label_] > z + _config_extract_sources_cosmoDC2.Dz_background
            else: 
                mask_z = dat_extract_mag_cut['z'] > z + _config_extract_sources_cosmoDC2.Dz_background
            
            if not use_only_member_galaxies:
                dat_extract_mag_cut_z_cut = dat_extract_mag_cut[mask_z]
            else:
                name = path_to_data + 'matched_pairs_Mfofcut.fits'
                dat = fits.open(name)
                dat_open= dat[1].data
                halo_id = int(dat_open['cat2_id'][dat_open['cat1_id'] == str(cluster_id)])
                mask_member_galaxy = np.isin(dat_extract_mag_cut['halo_id'], halo_id)
                print(np.sum(mask_member_galaxy))
                dat_extract_mag_cut_z_cut = dat_extract_mag_cut[mask_z & mask_member_galaxy]
                print(len(dat_extract_mag_cut_z_cut))
            
            if len(dat_extract_mag_cut_z_cut)!=0:
                bin_edges = _config_lensing_profiles.bin_edges
                dat_prf_pz = _utils_lensing_profiles.compute_lensing_profile(cluster_id, ra, dec, z, 
                                                                             bin_edges, label_, 
                                                                             dat_extract_mag_cut_z_cut, 
                                                                             cosmo)
            else:
                bin_edges = _config_lensing_profiles.bin_edges
                bins = [[bin_edges[i], bin_edges[i+1]] for i in range(len(bin_edges)-1)]
                bin_centers = np.mean(bins, axis=1)
                dat_prf_pz = [0*bin_centers, 0*bin_centers, 0*bin_centers, bin_centers]
                
            print(len(dat_extract_mag_cut_z_cut))
            print(dat_prf_pz)
            data_to_save = data_to_save + dat_prf_pz
        for s, n in enumerate(names): ind_profile[n].append(data_to_save[s])
                
#save_pickle(Table(ind_profile), name_save)

            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
