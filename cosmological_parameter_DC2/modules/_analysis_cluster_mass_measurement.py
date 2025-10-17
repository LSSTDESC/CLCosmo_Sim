import numpy as np
import _add_name_save

analysis = {}
#Modeling_choices

#impact halo profile
path_to_data = '../../CLCosmo_Sim_database/data_vary_fiducial_cosmology/'
analysis_1h_nfw = {'data_path': path_to_data + 'stacked_esd_profiles_redmapper_vary_Omega_m.pkl',
                                'photoz':'Truez',
                                'density_profile': 'nfw',
                                'cM_relation': 'None',
                                'two_halo': 'False',
                                'radius_min': 1.,
                                'radius_max': 3.5,}

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

analysis['cM_relation'] = [analysis_1h_nfw,
                           analysis_1h_nfw_Diemer15_true, 
                           analysis_1h_nfw_Duffy08_true, 
                           analysis_1h_nfw_Bhattacharya13_true, 
                           analysis_1h_nfw_Prada12_true]
# ""
for analysis_name in analysis.keys():
    for k in range(len(analysis[analysis_name])):
        analysis[analysis_name][k]['name'], analysis[analysis_name][k]['name_save'] = _add_name_save.add_name_save_cluster_mass_measurement(analysis[analysis_name][k])
