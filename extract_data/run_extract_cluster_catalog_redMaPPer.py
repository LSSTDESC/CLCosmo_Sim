import GCRCatalogs
import matplotlib.pyplot as plt
import pickle
import sys
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.table import Table
def load(filename, **kwargs):
    """Loads GalaxyCluster object to filename using Pickle"""
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)
def save_pickle(dat, filename, **kwargs):
     file = open(filename,'wb')
     pickle.dump(dat, file)
     file.close()
catalog = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_redmapper_v0.8.1')
quantity = ['cluster_id','ra', 'dec', 'redshift', 'redshift_err', 'richness', 'richness_err']
dat = catalog.get_quantities(quantity)
save_pickle(Table(dat), 'lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl')