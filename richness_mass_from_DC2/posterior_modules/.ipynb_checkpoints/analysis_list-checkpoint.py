import numpy as np

#analysis = {}

analysis_0 = {'name':'nfw_Duffy_1h_true_z',
              'type': 'WLxN',
              'fit_cosmo':False,
              'density_profile':'nfw',
              'cM_relation':'Duffy08',
              'two_halo':False,
              'hmf':'Despali16',
              'radius_max':5.5,
              'radius_min':1,
         'lensing_data':'/pbs/throng/lsst/users/cpayerne/CLMassDC2/notebooks/data_for_notebooks/stacked_esd_profiles_redmapper_true.pkl' }


analysis = [analysis_0]
#, analysis_1, analysis_2,analysis_3, analysis_4,
                   # analysis_10, analysis_11, analysis_12,analysis_13, analysis_14]