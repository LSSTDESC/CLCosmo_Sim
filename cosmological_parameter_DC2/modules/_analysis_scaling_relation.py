import numpy as np
import copy, glob, sys
import _add_name_save
import clmm
sys.path.append('../../')
import _redshift_richness_bins as analysis_
sys.path.append('../../lensing_profile_measurement/')
import _config_lensing_profiles
analysis = {}
path_to_data = '../../../CLCosmo_Sim_database/data/'
lensing_data_truez = path_to_data + 'stacked_esd_profiles_redmapper_true_full_coverage.pkl'
#baseline
edges = _config_lensing_profiles.bin_edges
analysis_Duffy08_baseline = {'type': 'WLxN',
                            'fit_cosmo':'True',
                            'cosmo_params': 'Om_s8',
                            'density_profile':'nfw',
                            'cM_relation':'Duffy08',
                            'two_halo':'False',
                            'hmf':'Despali16',
                            'radius_max':3.5,
                            'radius_min':1.,
                             'radius_edges':edges,
                             'radius_centers':np.array([(edges[i+1] + edges[i])/2 for i in range(len(edges)-1)]),
                            'photoz':'Truez',
                            'shear_richness_cov':'False',
                            'redshift_range':'Full',
                            'richness_range':'Full',
                            'redshift_bins': analysis_.Z_bin,
                            'richness_bins': analysis_.Obs_bin,
                             'redshift_corner': analysis_.z_corner,
                            'richness_corner': analysis_.rich_corner,
                            'redshift_bin_index': np.arange(len(analysis_.Z_bin)),
                            'richness_bin_index': np.arange(len(analysis_.Obs_bin)),
                            'redshift_corner_index': np.arange(len(analysis_.z_corner)),
                            'richness_corner_index': np.arange(len(analysis_.rich_corner)),
                            'lensing_data':lensing_data_truez,
                            'add_bias_lensing':'False',
                            'Gauss+SSC-CC_likelihood':'True',
                            'name_plot':'baseline'}

# analysis_Duffy08_baseline_low_z = analysis_Duffy08_baseline.copy()
# analysis_Duffy08_baseline_low_z['redshift_range'] = 'Partial'
# analysis_Duffy08_baseline_low_z['redshift_bin_index'] = np.arange(len(analysis_.Z_bin)-1)
# analysis_Duffy08_baseline_low_z['redshift_corner_index'] = np.arange(len(analysis_.z_corner)-1)
# analysis_Duffy08_baseline_low_z['name_plot'] = 'baseline - low z'

# analysis_Duffy08_baseline_low_richness = analysis_Duffy08_baseline.copy()
# analysis_Duffy08_baseline_low_richness['richness_range'] = 'Partial'
# analysis_Duffy08_baseline_low_richness['richness_bin_index'] = np.arange(len(analysis_.Obs_bin)-1)
# analysis_Duffy08_baseline_low_richness['richness_corner_index'] = np.arange(len(analysis_.rich_corner)-1)
# analysis_Duffy08_baseline_low_richness['name_plot'] = 'baseline - low richness'

baseline = [
            analysis_Duffy08_baseline,] 
            #analysis_Duffy08_baseline_low_z,
            #analysis_Duffy08_baseline_low_richness]
analysis['baseline'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in baseline]

analysis_list = []
for ff, k in enumerate(analysis.keys()):
    analysis_list += analysis[k]
analysis_WLxN = analysis_list

analysis_N = copy.deepcopy(analysis_WLxN)
for analysis_ in analysis_N:
    analysis_['type'] = 'N'
    _add_name_save.add_name_save_cluster_scaling_relation(analysis_)

analysis_WL = copy.deepcopy(analysis_WLxN)
for analysis_ in analysis_WL:
    analysis_['type'] = 'WL'
    _add_name_save.add_name_save_cluster_scaling_relation(analysis_)

analysis_MxN = copy.deepcopy(analysis_WLxN)
for analysis_ in analysis_MxN:
    analysis_['type'] = 'MxN'
    name_cl_mass, name_save_cl_mass = _add_name_save.add_name_save_cluster_mass_measurement(analysis_)
   # print(name_cl_mass)
    file_mass = glob.glob('../../cluster_mass_measurement/*')
    file_mass = [f.split('/cluster_mass_measurement/')[1].split('.pkl')[0] for f in file_mass]
    if name_cl_mass not in file_mass:
        print(name_cl_mass)
    analysis_['mass_file'] = name_save_cl_mass
    _add_name_save.add_name_save_cluster_scaling_relation(analysis_)

analysis_M = copy.deepcopy(analysis_WLxN)
for analysis_ in analysis_M:
    analysis_['type'] = 'M'
    name_cl_mass, name_save_cl_mass = _add_name_save.add_name_save_cluster_mass_measurement(analysis_)
    analysis_['mass_file'] = name_save_cl_mass
    _add_name_save.add_name_save_cluster_scaling_relation(analysis_)

config = {}
config['WL'] = analysis_WL
config['WLxN'] = analysis_WLxN
config['M'] = analysis_M
config['MxN'] = analysis_MxN
config['N'] = analysis_N