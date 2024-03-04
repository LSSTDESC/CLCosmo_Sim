import sys, os
import numpy as np
from astropy.table import QTable, Table, vstack, join
import pickle 
import pandas as pd
import clmm
sys.path.append('/pbs/throng/lsst/users/cpayerne/CLMassDC2/data/data_extraction/')
import photoz_utils

r"""
extract background galaxy catalog with qserv for:
cosmodc2:
- true shapes
- true redshift
and GCRCatalogs:
- photoz addons
"""

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
    tab['e1'], tab['e2'] = clmm.utils.compute_lensed_ellipticity(tab['e1_true'], tab['e2_true'], 
                                                                 tab['shear1'], tab['shear2'], 
                                                                 tab['kappa'])
    return tab
