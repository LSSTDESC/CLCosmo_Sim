import numpy as np
import _add_name_save

analysis = {}
#Modeling_choices

#impact halo profile
analysis_1h_nfw = {'data_path': '../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl',
                                'photoz':'Truez',
                                'density_profile': 'nfw',
                                'cM_relation': None,
                                'two_halo': False,
                                'radius_min': 1.5,
                                'radius_max': 5.5,}

#einasto
analysis_1h_einasto        = analysis_1h_nfw.copy()
analysis_1h_einasto['density_profile'] = 'Einasto'
#hernquist
analysis_1h_hernquist        = analysis_1h_nfw.copy()
analysis_1h_hernquist['density_profile'] = 'Hernquist'
analysis['halo_model'] = [analysis_1h_nfw, analysis_1h_einasto, analysis_1h_hernquist]

#impact c(M) relation
# diemer
analysis_1h_nfw_Diemer15_true        = analysis_1h_nfw.copy()
analysis_1h_nfw_Diemer15_true['cM_relation'] = 'Diemer15'
# ## Duffy08
analysis_1h_nfw_Duffy08_true        = analysis_1h_nfw.copy()
analysis_1h_nfw_Duffy08_true['cM_relation'] = 'Duffy08'

# ## Bhattacharya13
analysis_1h_nfw_Bhattacharya13_true = analysis_1h_nfw.copy()
analysis_1h_nfw_Bhattacharya13_true['cM_relation'] = 'Bhattacharya13'

### Duffy08
analysis_1h_nfw_Prada12_true        = analysis_1h_nfw.copy()
analysis_1h_nfw_Prada12_true['cM_relation'] = 'Prada12'

analysis['cM_relation'] = [analysis_1h_nfw_Diemer15_true, 
                           analysis_1h_nfw_Duffy08_true, 
                           analysis_1h_nfw_Bhattacharya13_true, 
                           analysis_1h_nfw_Prada12_true]
# ""
##impact of two halo term
analysis_2h_nfw_Diemer15_true = analysis_1h_nfw.copy()
analysis_2h_nfw_Diemer15_true['r_max'] = 15
analysis_2h_nfw_Diemer15_true['two_halo']=True
analysis_2h_nfw_Diemer15_true['cM_relation']='Diemer15'

analysis['2h'] = [analysis_2h_nfw_Diemer15_true]
# ""
# #impact of photoz
analysis_1h_nfw_Diemer15_bpz = analysis_1h_nfw_Diemer15_true.copy()
analysis_1h_nfw_Diemer15_bpz['data_path'] = '../../data/stacked_esd_profiles_redmapper_BPZ_full_coverage.pkl'
analysis_1h_nfw_Diemer15_bpz['photoz'] = 'BPZ'

analysis_1h_nfw_Diemer15_flex = analysis_1h_nfw_Diemer15_true.copy()
analysis_1h_nfw_Diemer15_flex['data_path'] = '../../data/stacked_esd_profiles_redmapper_flex_full_coverage.pkl'
analysis_1h_nfw_Diemer15_flex['photoz'] = 'flex'

analysis['photoz'] = [analysis_1h_nfw_Diemer15_bpz, analysis_1h_nfw_Diemer15_flex]
# ""
for analysis_name in analysis.keys():
    for k in range(len(analysis[analysis_name])):
        analysis[analysis_name][k]['name'], analysis[analysis_name][k]['name_save'] = _add_name_save.add_name_save_cluster_mass_measurement(analysis[analysis_name][k])