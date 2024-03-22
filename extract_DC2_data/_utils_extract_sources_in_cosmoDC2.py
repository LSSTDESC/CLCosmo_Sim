import sys, os
import numpy as np
from astropy.table import QTable, Table, vstack, join
import pickle 
import pandas as pd
import clmm
import cmath

r"""
extract background galaxy catalog with qserv for:
cosmodc2:
- true shapes
- true redshift
and GCRCatalogs:
- photoz addons
"""
def _fix_axis_ratio(q_bad):
    # back out incorrect computation of q using Johnsonb function
    e_jb = np.sqrt((1 - q_bad**2)/(1 + q_bad**2))
    q_new = np.sqrt((1 - e_jb)/(1 + e_jb)) # use correct relationship to compute q from e_jb 
    return q_new

def _fix_ellipticity_disk_or_bulge(ellipticity):
    # back out incorrect computation of q using Johnsonb function 
    q_bad = (1-ellipticity)/(1+ellipticity) #use default e definition to calculate q
    # q_bad incorrectly computed from e_jb using q_bad = sqrt((1 - e_jb^2)/(1 + e_jb^2))
    q_new = _fix_axis_ratio(q_bad)
    e_new = (1 - q_new)/(1 + q_new)  # recompute e using default (1-q)/(1+q) definition
    return e_new

def correct_shear_ellipticity(ellipticity_uncorr_e1, ellipticity_uncorr_e2):
    ellipticity_uncorr_norm = (ellipticity_uncorr_e1**2+ellipticity_uncorr_e2**2)**.5
    complex_ellipticity_uncorr = ellipticity_uncorr_e1 + 1j*ellipticity_uncorr_e2
    phi = np.array([cmath.phase(c) for c in complex_ellipticity_uncorr])
    ellipticity_corr_norm = _fix_ellipticity_disk_or_bulge(ellipticity_uncorr_norm)
    ellipticity_corr = ellipticity_corr_norm*np.exp(1j*phi)
    ellipticity_corr_e1, ellipticity_corr_e2 = ellipticity_corr.real, ellipticity_corr.imag
    return ellipticity_corr_e1, ellipticity_corr_e2

def extract_photoz(id_gal, healpix=None, GCRcatalog=None):
    r"""
    extract background galaxy catalog with GCRcatalog (healpix subdivision)
    Attributes:
    -----------
    id_gal: array
        background galaxy id
    healpix: array
        list of healpix pixels where to find galaxies
    GCRcatalog: GCRcatalog
        background galaxy GCRcatalog object
    Returns:
    --------
    tab_astropy_ordered: Table
        photoz informations 
    """
    Table_id_gal = Table()
    Table_id_gal['galaxy_id'] = id_gal
    quantities_photoz=['photoz_pdf','photoz_mean','photoz_mode',
                       'photoz_odds','galaxy_id']
    z_bins = GCRcatalog.photoz_pdf_bin_centers
    z_bins[0] = 1e-7
    for n, hp in enumerate(np.array(healpix)):
        tab = GCRcatalog.get_quantities(quantities_photoz, 
                                        native_filters=['healpix_pixel=='+str(hp)])
        tab_astropy = Table()
        tab_astropy['galaxy_id']   = tab['galaxy_id']
        tab_astropy['photoz_pdf']  = tab['photoz_pdf']
        tab_astropy['photoz_mean'] = tab['photoz_mean']
        tab_astropy['photoz_mode'] = tab['photoz_mode']
        tab_astropy['photoz_odds'] = tab['photoz_odds']
        mask_id=np.isin(np.array(tab_astropy['galaxy_id']), id_gal)
        tab_astropy=tab_astropy[mask_id]
        if n==0: 
            table_photoz=tab_astropy
        else: 
            tab_=vstack([table_photoz,tab_astropy])
            table_photoz = tab_
            
    n_gal = len(table_photoz['galaxy_id'])
    table_photoz['pzbins'] = np.array([z_bins for i in range(n_gal)])
    table_photoz_ordered = join(table_photoz, Table_id_gal, keys='galaxy_id')
    return table_photoz_ordered

def extract(lens_redshift=None,
            qserv_query=None, GCRcatalog=None, conn_qserv=None,
           cosmo=None):
    r"""
    extract background galaxy catalog
    Attributes:
    -----------
    lens_redshift: float
        lens redshift
    lens_ra: float
        lens right ascension
    lens_dec: float
        lens declinaison
    Returns:
    --------
    cl: GalaxyCluster
        background galaxy catalog
    """
    #qserv
    query_mysql = qserv_query
    tab = pd.read_sql_query(query_mysql, conn_qserv)
    try: 
        tab = QTable.from_pandas(tab)
    except: 
        print('no data')
        return None
    #compute reduced shear and ellipticities
    tab['g1'], tab['g2'] = clmm.utils.convert_shapes_to_epsilon(tab['shear1'],tab['shear2'], 
                                                                shape_definition = 'shear',
                                                                kappa = tab['kappa'])
    ellipticity_uncorr_e1 = tab['e1_true_uncorr']
    ellipticity_uncorr_e2 = tab['e2_true_uncorr']
    ellipticity_corr_e1, ellipticity_corr_e2 = correct_shear_ellipticity(ellipticity_uncorr_e1, ellipticity_uncorr_e2)
    tab['e1_true'] = ellipticity_corr_e1
    tab['e2_true'] = ellipticity_corr_e2
    tab['e1'], tab['e2'] = clmm.utils.compute_lensed_ellipticity(tab['e1_true'], tab['e2_true'], 
                                                                 tab['shear1'], tab['shear2'], 
                                                                 tab['kappa'])
    return tab
