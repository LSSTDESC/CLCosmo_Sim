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
import pyccl as ccl
sys.path.append('../')
import _utils_extract_sources_in_cosmoDC2 
import _config_extract_sources_cosmoDC2
sys.path.append('../../lensing_profile_measurement')
import _config_lensing_profiles
import _utils_lensing_profiles
sys.path.append('../../')
import _redshift_richness_bins
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
#connection with qserv
import mysql
from mysql.connector import Error
import argparse

path_to_data = '../../../CLCosmo_Sim_database/data_vary_fuducial_cosmology/'


def collect_argparser():
    parser = argparse.ArgumentParser(description="Let's extract source cat")
    parser.add_argument("--index_split", type=int, required=True)
    parser.add_argument("--n_splits", type=int, required=True)
    parser.add_argument("--lens_catalog_name", type=str, required=False, default=path_to_data+'lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl')
    parser.add_argument("--compute_individual_lensing_profile", type=str, required=True)
    parser.add_argument("--save_catalog", type=str, required=False, default='False')
    parser.add_argument("--where_to_save_catalog", type=str, required=False, default='/sps/lsst/users/cpayerne/CLMassDC2/cosmoDC2/redmapper_clusters_new/')
    return parser.parse_args()

H0_true = 71
h = H0_true/100
Omega_b_true = 0.02258 / (h**2)
Omega_c_true = 0.1109 / (h**2)
Omega_m_true = Omega_b_true + Omega_c_true
sigma8_true = 0.8
ns_true = 0.963

cosmo_clmm_fid_distance = Cosmology(H0 = H0_true, Omega_dm0 = Omega_c_true, Omega_b0 = Omega_b_true, Omega_k0 = 0.0)
cosmo_clmm_fid_distance.be_cosmo = ccl.Cosmology(Omega_c = Omega_c_true, Omega_b = Omega_b_true,
                                        h = H0_true/100, sigma8 = sigma8_true, 
                                        n_s=ns_true, w0=-1, wa=0)

#select galaxy clusters
_config_extract_sources_in_cosmoDC2 = collect_argparser()
lens_catalog_name=_config_extract_sources_in_cosmoDC2.lens_catalog_name
lens_catalog=np.load(lens_catalog_name, allow_pickle=True)
mask_select=(lens_catalog['richness'] > 20)*(lens_catalog['redshift'] > .2)*(lens_catalog['redshift'] <= 1)
lens_catalog=lens_catalog[mask_select]
lens_catalog_truncated=lens_catalog
use_only_member_galaxies = False

n_cl=len(lens_catalog_truncated)
index_cl=np.arange(n_cl)
split_lists=np.array_split(index_cl, _config_extract_sources_in_cosmoDC2.n_splits)
lens_catalog_truncated=lens_catalog_truncated[np.array(split_lists[_config_extract_sources_in_cosmoDC2.index_split])]
n_cl_to_extract=len(lens_catalog_truncated)

if lens_catalog_name=='../../../CLCosmo_Sim_database/data/random_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl': 
    suff_='_random'
elif lens_catalog_name=='../../../CLCosmo_Sim_database/data/lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl': 
    suff_=''
path_where_to_save_profiles = '../../../CLCosmo_Sim_database/data_vary_fiducial_cosmology/ind_gammat_profile_redmapper_per_cluster_index/'

name_save = path_where_to_save_profiles+f'ind_gammat_profile_redmapper'+suff_+f'_split={_config_extract_sources_in_cosmoDC2.index_split}_nsplits={_config_extract_sources_in_cosmoDC2.n_splits}_ncl={n_cl_to_extract}.pkl'
#load source catalogs
healpix_dc2 = GCRCatalogs.load_catalog("cosmoDC2_v1.1.4_image").get_catalog_info()['healpix_pixels']
photoz=False
if photoz == True:
    gc_bpz  = "cosmoDC2_v1.1.4_image_with_photozs_v1"
    gc_flex = "cosmoDC2_v1.1.4_image_with_photozs_flexzboost_v1"
    #list of healpix in cosmoDC2
    z_bins  = GCRCatalogs.load_catalog(gc_flex).photoz_pdf_bin_centers
    z_bins[0] = 1e-7
    photoz_gc=[gc_bpz, gc_flex]
    photoz_label=['_bpz', '_flex']

def qserv_query(lens_z, theta_max, ra, dec):
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
    n = ['DSt', 'DSx', 'W_l', 'radius']
    name_per_profile = ['id', 'ra', 'dec', 'z', 'richness'] + n
    ind_profile = {n:[] for n in name_per_profile}

for ncl, lens in enumerate(lens_catalog_truncated):
    
    if len(lens_catalog_truncated) == 0: continue

    z, ra, dec=lens['redshift'], lens['ra'], lens['dec']
    cluster_id=lens['cluster_id']
    l=lens['richness']
    #cluster member galaxies
    id_member_galaxy = lens_catalog['id_member'][ncl]
    ra_member_galaxy = lens_catalog['ra_member'][ncl]
    dec_member_galaxy = lens_catalog['dec_member'][ncl]
    lens_distance = cosmo_clmm_fid_distance.eval_da(z)

    print(f'====> cluster n°{cluster_id}, at redshift z = {z:.2f} and distance DA = {lens_distance:.2f} Mpc')
    print(f'====> ra = {ra:.2f} deg and dec = {dec:.2f} deg')
    print(f'====> richness = {l:.2f}')
    print(f'====> ra mean = {np.mean(ra_member_galaxy):.2f} deg and dec mean = {np.mean(dec_member_galaxy):.2f} deg')

    lens_z_bin_0 = _redshift_richness_bins.z_corner[_redshift_richness_bins.z_corner <= z][-1]
    lens_z_bin_1 = _redshift_richness_bins.z_corner[_redshift_richness_bins.z_corner >= z][0]
    lens_z_mid = (lens_z_bin_0 + lens_z_bin_1)/2

    lens_distance_z_bin_0=cosmo_clmm_fid_distance.eval_da(lens_z_bin_0)
    lens_distance_z_bin_1=cosmo_clmm_fid_distance.eval_da(lens_z_bin_1)
    lens_distance_z_mid=cosmo_clmm_fid_distance.eval_da(lens_z_mid)
    theta_max = (11/lens_distance_z_mid) * (180/np.pi)
    bin_edges_R = clmm.dataops.make_bins(0.5, 10, 15, method='evenlog10width')
    bin_edges = (bin_edges_R / lens_distance_z_mid) * (10800/np.pi)
    print('bin edges in arcmin = ' ,bin_edges)

    conn = mysql.connector.connect(host='ccqserv201', user='qsmaster', port=30040)
    cursor = conn.cursor(dictionary=True, buffered=True)
    #extract background galaxies with qserv (only true photoz)
    print('==> extracting true redshift + shapes infos (Qserv)')
    dat_extract=_utils_extract_sources_in_cosmoDC2.extract(qserv_query = qserv_query(z, theta_max, ra, dec),
                                        conn_qserv=conn, cosmo=None)
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
        
    if _config_extract_sources_in_cosmoDC2.compute_individual_lensing_profile=='True':
        dat_extract_mag_cut = _config_extract_sources_cosmoDC2.source_selection_magnitude(dat_extract)
        for label_ in _config_lensing_profiles.label_pz:
            if label_ != 'true':
                mask_z = dat_extract_mag_cut['p_background'+'_'+label_] > _config_extract_sources_cosmoDC2.p_background_min
                mask_z = dat_extract_mag_cut['photoz_mean'+'_'+label_] > lens_z_mid + _config_extract_sources_cosmoDC2.Dz_background
            else: 
                mask_z = dat_extract_mag_cut['z'] > lens_z_mid + _config_extract_sources_cosmoDC2.Dz_background
            
            if not use_only_member_galaxies:
                dat_extract_mag_cut_z_cut = dat_extract_mag_cut[mask_z]
                
            else:
                name = path_to_data + 'matched_pairs_Mfofcut.fits'
                dat = fits.open(name)
                dat_open= dat[1].data
                halo_id = int(dat_open['cat2_id'][dat_open['cat1_id'] == str(cluster_id)])
                mask_member_galaxy = np.isin(dat_extract_mag_cut['halo_id'], halo_id)
                dat_extract_mag_cut_z_cut = dat_extract_mag_cut[mask_z & mask_member_galaxy]
            print('==> size of the background galaxy catalog = '  + str(len(dat_extract_mag_cut_z_cut)))
            def compute_profile():
                if len(dat_extract_mag_cut_z_cut)!=0:
                    dat_prf_pz = _utils_lensing_profiles.compute_tangential_shear_lensing_profile(cluster_id, ra, dec, z, 
                                                                             bin_edges, label_, 
                                                                             dat_extract_mag_cut_z_cut)
                else:
                    bins = [[bin_edges[i], bin_edges[i+1]] for i in range(len(bin_edges)-1)]
                    bin_centers = np.mean(bins, axis=1)
                    dat_prf_pz = [0*bin_centers, 0*bin_centers, 0*bin_centers, bin_centers]
                return dat_prf_pz

            prf_fid = compute_profile()
            print(prf_fid[0])
                
            data_to_save = [cluster_id, ra, dec, z, lens['richness']] + prf_fid
        for s, n in enumerate(name_per_profile): 
            ind_profile[n].append(data_to_save[s])
        
    #if ncl > 10: break

print(len(ind_profile['richness']))
print(ind_profile.keys())

save_pickle(Table(ind_profile), name_save)

            

    

    

    

    

    

    

    

    
