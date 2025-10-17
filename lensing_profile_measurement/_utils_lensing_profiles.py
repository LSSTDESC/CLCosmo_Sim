import clmm
import numpy as np
from astropy.table import QTable, Table, vstack, join, hstack
from scipy.stats import binned_statistic

def compute_excess_surface_density_lensing_profile(cluster_id, ra, dec, z, bin_edges, label, table, cosmo):
    
        cl = clmm.galaxycluster.GalaxyCluster('halo', ra, dec, z, clmm.gcdata.GCData(Table(table)))
        theta1, g_t, g_x = cl.compute_tangential_and_cross_components(is_deltasigma=False, cosmo=None)
        data_prf = []

        if label=='true': sigma_c = cosmo.eval_sigma_crit(z, cl.galcat['z'])
        elif label=='flex': sigma_c = cl.galcat['sigmac_photoz_flex']
        elif label=='bpz': sigma_c = cl.galcat['sigmac_photoz_bpz']

        cl.galcat['dst'] = sigma_c*cl.galcat['et']
        cl.galcat['dsx'] = sigma_c*cl.galcat['ex']
        cl.galcat['w_ls'] = sigma_c**(-2.)
        #print(max(cl.galcat['theta'] * cosmo.eval_da(z)))
        ce = clmm.ClusterEnsemble('id', [])
        p = ce.make_individual_radial_profile(cl, 'Mpc', bins=bin_edges, error_model='ste',
                                           cosmo=cosmo, tan_component_in='dst', cross_component_in='dsx',
                                           tan_component_out='gt', cross_component_out='gx',
                                           tan_component_in_err=None, cross_component_in_err=None,
                                           weights_in='w_ls', weights_out='W_l')
        data = ce.data[0]
        data_prf = [data['gt'], data['gx'], data['W_l'], data['radius']]
        return data_prf

def compute_tangential_shear_lensing_profile(cluster_id, ra, dec, z, bin_edges, label, table):
    
        cl = clmm.galaxycluster.GalaxyCluster('halo', ra, dec, z, clmm.gcdata.GCData(Table(table)))
        theta1, g_t, g_x = cl.compute_tangential_and_cross_components(is_deltasigma=False, cosmo=None)
        data_prf = []

        cl.galcat['dst'] = cl.galcat['et']
        cl.galcat['dsx'] = cl.galcat['ex']
        cl.galcat['w_ls'] = 1
        print(max(cl.galcat['theta']) * (10800/np.pi))
        ce = clmm.ClusterEnsemble('id', [])
        p = ce.make_individual_radial_profile(cl, 'arcmin', bins=bin_edges, error_model='ste',
                                           cosmo=None, tan_component_in='dst', cross_component_in='dsx',
                                           tan_component_out='gt', cross_component_out='gx',
                                           tan_component_in_err=None, cross_component_in_err=None,
                                           weights_in='w_ls', weights_out='W_l')
        data = ce.data[0]
        data_prf = [data['gt'], data['gx'], data['W_l'], data['radius']]
        return data_prf
    
        
        