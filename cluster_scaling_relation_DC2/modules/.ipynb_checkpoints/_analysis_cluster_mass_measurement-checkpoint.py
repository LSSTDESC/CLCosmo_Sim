import numpy as np

def add_name(analysis):
    name = 'cluster-masses_'
    add_lensing_infos = '1-halo=' + analysis['halo_profile']
    if analysis['cM'] != None: add_lensing_infos = add_lensing_infos + '+c-M=' + analysis['cM']
    else: add_lensing_infos = add_lensing_infos + '+c-M=None' 
    if analysis['two_halo_term']: add_lensing_infos = add_lensing_infos + '_' + '+2-halo'
    add_lensing_infos = add_lensing_infos + '_photoz='+analysis['photoz']
    name = name + add_lensing_infos
    analysis['name'] = name
    analysis['name_save'] = '../../cluster_mass_measurement/' + analysis['name'] + '.pkl'
    return analysis

analysis = {}
#Modeling_choices

#impact halo profile
analysis_1h_nfw = {'data_path': '../../data/stacked_esd_profiles_redmapper_true.pkl',
                                'photoz':'Truez',
                                'halo_profile': 'nfw',
                                'cM': None,
                                'two_halo_term': False,
                                'r_min': 1.5,
                                'r_max': 5.5,}
add_name(analysis_1h_nfw)
#einasto
analysis_1h_einasto        = analysis_1h_nfw.copy()
analysis_1h_einasto['halo_profile'] = 'Einasto'
add_name(analysis_1h_einasto)
#hernquist
analysis_1h_hernquist        = analysis_1h_nfw .copy()
analysis_1h_hernquist['halo_profile'] = 'Hernquist'
add_name(analysis_1h_hernquist)
analysis['halo_model'] = [analysis_1h_nfw, analysis_1h_einasto, analysis_1h_hernquist]

#impact c(M) relation
# diemer

analysis_1h_nfw_Diemer15_true        = analysis_1h_nfw.copy()
analysis_1h_nfw_Diemer15_true['cM'] = 'Diemer15'
add_name(analysis_1h_nfw_Diemer15_true)
# ## Duffy08
analysis_1h_nfw_Duffy08_true        = analysis_1h_nfw.copy()
analysis_1h_nfw_Duffy08_true['cM'] = 'Duffy08'
add_name(analysis_1h_nfw_Duffy08_true)

# ## Bhattacharya13
analysis_1h_nfw_Bhattacharya13_true = analysis_1h_nfw.copy()
analysis_1h_nfw_Bhattacharya13_true['cM'] = 'Bhattacharya13'
add_name(analysis_1h_nfw_Bhattacharya13_true)

### Duffy08
analysis_1h_nfw_Prada12_true        = analysis_1h_nfw.copy()
analysis_1h_nfw_Prada12_true['cM'] = 'Prada12'
add_name(analysis_1h_nfw_Prada12_true)

analysis['cM'] = [analysis_1h_nfw_Diemer15_true, analysis_1h_nfw_Duffy08_true, analysis_1h_nfw_Bhattacharya13_true, analysis_1h_nfw_Prada12_true]
# ""
##impact of two halo term
analysis_2h_nfw_Diemer15_true = analysis_1h_nfw.copy()
analysis_2h_nfw_Diemer15_true['r_max'] = 15
analysis_2h_nfw_Diemer15_true['two_halo_term']=True
analysis_2h_nfw_Diemer15_true['cM']='Diemer15'
add_name(analysis_2h_nfw_Diemer15_true)

analysis['2h'] = [analysis_2h_nfw_Diemer15_true]
# ""
# #impact of photoz
analysis_1h_nfw_Diemer15_bpz = analysis_1h_nfw.copy()
analysis_1h_nfw_Diemer15_bpz['data_path'] = '../../data/stacked_esd_profiles_redmapper_BPZ.pkl'
analysis_1h_nfw_Diemer15_bpz['photoz'] = 'BPZ'
add_name(analysis_1h_nfw_Diemer15_bpz)

analysis_1h_nfw_Diemer15_flex = analysis_1h_nfw.copy()
analysis_1h_nfw_Diemer15_flex['data_path'] = '../../data/stacked_esd_profiles_redmapper_flex.pkl'
analysis_1h_nfw_Diemer15_flex['photoz'] = 'flex'
add_name(analysis_1h_nfw_Diemer15_flex)

analysis['photoz'] = [analysis_1h_nfw_Diemer15_bpz, analysis_1h_nfw_Diemer15_flex]
# ""

