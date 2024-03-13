import numpy as np

analysis = {}

def add_name(analysis):
    #kind of fit
    name = 'MCMC_'
    if analysis['fit_cosmo']: fit_ = 'fit_m-r+cosmo_'
    else: fit_ = 'fit_m-r_'
    
    #from which observable
    name = name+ fit_
    name = name + analysis['type'] + '_'
    
    #if lensing
    if (analysis['type'] == 'WL') or (analysis['type'] == 'WLxN'):
        
        add_lensing_infos = '1-halo=' + analysis['density_profile'] + '+c-M=' + analysis['cM_relation']
        if analysis['two_halo']:
            add_lensing_infos = add_lensing_infos + '_' + '+2-halo'
        if analysis['shear_richness_cov']:
            add_lensing_infos = add_lensing_infos + '+cov(g,richness)'
        add_lensing_infos = add_lensing_infos + '_photoz='+analysis['photoz']
        name = name + add_lensing_infos
    #kind of hmf
    name = name + '_hmf='+analysis['hmf']
    analysis['name'] = name
    return analysis
    
lensing_data_truez = '../../data/stacked_esd_profiles_redmapper_true.pkl'
lensing_data_bpz = '../../data/stacked_esd_profiles_redmapper_BPZ.pkl'
lensing_data_flex = '../../data/stacked_esd_profiles_redmapper_flex.pkl'
#impact_c_m
analysis_Duffy08 = {
              'type': 'WLxN',
              'fit_cosmo':False,
              'density_profile':'nfw',
              'cM_relation':'Duffy08',
              'two_halo':False,
              'hmf':'Despali16',
              'radius_max':5.5,
              'radius_min':1,
              'photoz':'Truez',
              'shear_richness_cov':False,
              'lensing_data':lensing_data_truez,
'name_plot':'nfw - Duffy08'}
analysis_Diemer15 = analysis_Duffy08.copy()
analysis_Diemer15['cM_relation'] = 'Diemer15'
analysis_Diemer15['name_plot'] = 'nfw - Diemer15'
analysis_Bhattacharya13 = analysis_Duffy08.copy()
analysis_Bhattacharya13['cM_relation'] = 'Bhattacharya13'
analysis_Bhattacharya13['name_plot'] = 'nfw - Bhattacharya13'
analysis_Prada12 = analysis_Duffy08.copy()
analysis_Prada12['cM_relation'] = 'Prada12'
analysis_Prada12['name_plot'] = 'nfw - Prada12'
cM = [analysis_Duffy08, analysis_Diemer15, analysis_Prada12, analysis_Bhattacharya13]
analysis['c_M'] = [add_name(a) for a in cM]

#hmf
analysis_Bocquet16 = analysis_Duffy08.copy()
analysis_Bocquet16['hmf'] = 'Bocquet16'
analysis_Bocquet16['name_plot'] = 'nfw - Bocquet16'
analysis_Bocquet20 = analysis_Duffy08.copy()
analysis_Bocquet20['hmf'] = 'Bocquet20'
analysis_Bocquet20['name_plot'] = 'nfw - Bocquet20'
analysis_hmf = [analysis_Bocquet16, analysis_Bocquet20]
analysis['hmf'] = [add_name(a) for a in analysis_hmf]

#photoz
analysis_BPZ = analysis_Duffy08.copy()
analysis_BPZ['lensing_data'] = lensing_data_bpz
analysis_BPZ['name_plot'] = 'BPZ'
analysis_flex = analysis_Duffy08.copy()
analysis_flex['lensing_data'] = lensing_data_flex
analysis_flex['name_plot'] = 'FleXZBoost'
photoz = [analysis_BPZ, analysis_flex]
analysis['photoz'] = [add_name(a) for a in photoz]


#GammaCov covariance
analysis_GammaCov1 = analysis_Duffy08.copy()
analysis_GammaCov1['shear_richness_cov'] = True
analysis_GammaCov1['name_plot'] = r'nfw + Cov($\Delta\Sigma,\lambda$)'
GammaCov = [analysis_GammaCov1]
analysis['GammaLambda'] = [add_name(a) for a in GammaCov]

analysis_list = []
for i, k in enumerate(analysis.keys()):
    analysis_list += analysis[k]
    
for i in range(len(analysis_list)):
    analysis_list[i]['type'] = 'WLxN'
    analysis_list[i] = add_name(analysis_list[i])
    
               