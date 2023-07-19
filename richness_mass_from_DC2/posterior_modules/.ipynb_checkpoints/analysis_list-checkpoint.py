import numpy as np

analysis = {}

analysis_0 = {'type': 'WLxN','fit_cosmo':True}
analysis_1 = {'type': 'MxN','fit_cosmo':True}
analysis_2 = {'type': 'N','fit_cosmo':True}
analysis_3 = {'type': 'M','fit_cosmo':True}
analysis_4 = {'type': 'WL','fit_cosmo':True}

analysis_10 = {'type': 'WLxN','fit_cosmo':False}
analysis_11 = {'type': 'MxN','fit_cosmo':False}
analysis_12 = {'type': 'N','fit_cosmo':False}
analysis_13 = {'type': 'M','fit_cosmo':False}
analysis_14 = {'type': 'WL','fit_cosmo':False}

analysis['type'] = [analysis_0, analysis_1, analysis_2,analysis_3, analysis_4,
                    analysis_10, analysis_11, analysis_12,analysis_13, analysis_14]