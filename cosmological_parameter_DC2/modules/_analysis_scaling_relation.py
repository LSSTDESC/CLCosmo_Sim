import numpy as np
import copy, glob, sys
import _add_name_save
sys.path.append('../../')
import _redshift_richness_bins as analysis_
analysis = {}
path_to_data = '../../../CLCosmo_Sim_database/data/'
lensing_data_truez = path_to_data + 'stacked_esd_profiles_redmapper_true_full_coverage.pkl'
lensing_data_bpz = path_to_data + 'stacked_esd_profiles_redmapper_BPZ_full_coverage.pkl'
lensing_data_flex = path_to_data + 'stacked_esd_profiles_redmapper_flex_full_coverage.pkl'
#baseline
analysis_Duffy08_baseline = {'type': 'WLxN',
                            'fit_cosmo':'False',
                            'density_profile':'nfw',
                            'cM_relation':'Duffy08',
                            'two_halo':'False',
                            'hmf':'Despali16',
                            'radius_max':3.5,
                            'radius_min':1.,
                            'photoz':'Truez',
                            'shear_richness_cov':'False',
                            'redshift_range':'Full',
                            'richness_range':'Full',
                            'redshift_bin_index': np.arange(len(analysis_.Z_bin)),
                            'richness_bin_index': np.arange(len(analysis_.Obs_bin)),
                            'redshift_corner_index': np.arange(len(analysis_.z_corner)),
                            'richness_corner_index': np.arange(len(analysis_.rich_corner)),
                            'lensing_data':lensing_data_truez,
                            'add_bias_lensing':'False',
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

#baseline = [
            #analysis_Duffy08_baseline, 
            #analysis_Duffy08_baseline_low_z,
            #analysis_Duffy08_baseline_low_richness]
#analysis['baseline'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in baseline]

# impact c-m
# analysis_None = analysis_Duffy08_baseline.copy()
# analysis_None['cM_relation'] = 'None'
# analysis_None['name_plot'] = 'nfw - cM=None'
# analysis_Diemer15 = analysis_Duffy08_baseline.copy()
# analysis_Diemer15['cM_relation'] = 'Diemer15'
# analysis_Diemer15['name_plot'] = 'nfw - Diemer15'
# analysis_Bhattacharya13 = analysis_Duffy08_baseline.copy()
# analysis_Bhattacharya13['cM_relation'] = 'Bhattacharya13'
# analysis_Bhattacharya13['name_plot'] = 'nfw - Bhattacharya13'
# analysis_Prada12 = analysis_Duffy08_baseline.copy()
# analysis_Prada12['cM_relation'] = 'Prada12'
# analysis_Prada12['name_plot'] = 'nfw - Prada12'
# cM = [analysis_None, analysis_Diemer15, analysis_Prada12, analysis_Bhattacharya13]
# analysis['c_M'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in cM]

# #impact halo profile
# analysis_NFW = analysis_Duffy08_baseline.copy()
# analysis_NFW['density_profile'] = 'nfw'
# analysis_NFW['cM_relation'] = 'None'

# analysis_Hernquist = analysis_Duffy08_baseline.copy()
# analysis_Hernquist['density_profile'] = 'hernquist'
# analysis_Hernquist['cM_relation'] = 'None'

# analysis_Einasto = analysis_Duffy08_baseline.copy()
# analysis_Einasto['density_profile'] = 'einasto'
# analysis_Einasto['cM_relation'] = 'None'

# halomodel = [analysis_NFW, analysis_Hernquist, analysis_Einasto]
# analysis['halo_model'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in halomodel]

# #two-halo term
# analysis_2h = analysis_Duffy08_baseline.copy()
# analysis_2h['two_halo'] = 'True'
# analysis_2h['radius_max'] = 10
# analysis_2h['name_plot'] = 'nfw - Diemer15 (2h)'
# analysis_2hshort = analysis_Duffy08_baseline.copy()
# analysis_2hshort['two_halo'] = 'True'
# analysis_2hshort['radius_max'] = 3.5
# analysis_2hshort['name_plot'] = 'nfw - Diemer15 (2h) short'
# twoh = [analysis_2h, analysis_2hshort]
# analysis['2h']  = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in twoh]

# #hmf
# analysis_Bocquet16 = analysis_Diemer15.copy()
# analysis_Bocquet16['hmf'] = 'Bocquet16'
# analysis_Bocquet16['name_plot'] = 'nfw - Bocquet16'
# analysis_hmf = [analysis_Bocquet16]
# analysis['hmf'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in analysis_hmf]

# #photoz
# analysis_BPZ = analysis_Duffy08_baseline.copy()
# analysis_BPZ['lensing_data'] = lensing_data_bpz
# analysis_BPZ['name_plot'] = 'BPZ'
# analysis_BPZ['photoz'] = 'BPZ'

# analysis_BPZb = analysis_Duffy08_baseline.copy()
# analysis_BPZb['lensing_data'] = lensing_data_bpz
# analysis_BPZb['name_plot'] = 'BPZ(1+b)'
# analysis_BPZb['photoz'] = 'BPZ'
# analysis_BPZb['add_bias_lensing'] = 'True'

# analysis_flex = analysis_Duffy08_baseline.copy()
# analysis_flex['lensing_data'] = lensing_data_flex
# analysis_flex['name_plot'] = 'FleXZBoost'
# analysis_flex['photoz'] = 'flex'

# analysis_flexb = analysis_Duffy08_baseline.copy()
# analysis_flexb['lensing_data'] = lensing_data_flex
# analysis_flexb['name_plot'] = 'FleXZBoost(1+b)'
# analysis_flexb['photoz'] = 'flex'
# analysis_flexb['add_bias_lensing'] = 'True'

# photoz = [analysis_BPZ,analysis_BPZb,analysis_flex,analysis_flexb]
# analysis['photoz'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in photoz]

# impact radius cut
# analysis_4 = analysis_Duffy08_baseline.copy()
# analysis_4['radius_max'] = 4
# analysis_4['name_plot'] = r'baseline - Rmax = 4'
# analysis_45 = analysis_Duffy08_baseline.copy()
# analysis_45['radius_max'] = 4.5
# analysis_45['name_plot'] = r'baseline - Rmax = 4.5'
# analysis_5 = analysis_Duffy08_baseline.copy()
# analysis_5['radius_max'] = 5
# analysis_5['name_plot'] = r'baseline - Rmax = 5'
# analysis_55 = analysis_Duffy08_baseline.copy()
# analysis_55['radius_max'] = 5.5
# analysis_55['name_plot'] = r'baseline - Rmax = 5.5'
# analysis_8 = analysis_Duffy08_baseline.copy()
# analysis_8['radius_max'] = 8
# analysis_8['name_plot'] = r'baseline - Rmax = 8'
# analysis_10 = analysis_Duffy08_baseline.copy()
# analysis_10['radius_max'] = 10
# analysis_10['name_plot'] = r'baseline - Rmax = 10'
# rmax = [analysis_4, analysis_45, analysis_5, analysis_55, analysis_8, analysis_10]
# analysis['Rmax'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in rmax]

# #GammaCov covariance
analysis_GammaCov1 = analysis_Duffy08_baseline.copy()
analysis_GammaCov1['shear_richness_cov'] = 'True'
analysis_GammaCov1['photoz'] = 'Truez'
analysis_GammaCov1['name_plot'] = r'nfw + Cov($\Delta\Sigma,\lambda$)'

analysis_GammaCov2 = analysis_Duffy08_baseline.copy()
analysis_GammaCov2['shear_richness_cov'] = 'True'
analysis_GammaCov2['photoz'] = 'BPZ'
analysis_GammaCov2['lensing_data'] = lensing_data_bpz
analysis_GammaCov2['name_plot'] = r'nfw + Cov($\Delta\Sigma,\lambda$) + BPZ'

analysis_GammaCov2bis = analysis_Duffy08_baseline.copy()
analysis_GammaCov2bis['shear_richness_cov'] = 'True'
analysis_GammaCov2bis['photoz'] = 'BPZ'
analysis_GammaCov2bis['lensing_data'] = lensing_data_bpz
analysis_GammaCov2bis['add_bias_lensing'] = 'True'
analysis_GammaCov2bis['name_plot'] = r'nfw + Cov($\Delta\Sigma,\lambda$) + BPZ (1+b)'

analysis_GammaCov3 = analysis_Duffy08_baseline.copy()
analysis_GammaCov3['shear_richness_cov'] = 'True'
analysis_GammaCov3['photoz'] = 'flex'
analysis_GammaCov3['lensing_data'] = lensing_data_flex
analysis_GammaCov3['name_plot'] = r'nfw + Cov($\Delta\Sigma,\lambda$) + FlexZBoost'

analysis_GammaCov3bis = analysis_Duffy08_baseline.copy()
analysis_GammaCov3bis['shear_richness_cov'] = 'True'
analysis_GammaCov3bis['photoz'] = 'flex'
analysis_GammaCov3bis['lensing_data'] = lensing_data_flex
analysis_GammaCov3bis['add_bias_lensing'] = 'True'
analysis_GammaCov3bis['name_plot'] = r'nfw + Cov($\Delta\Sigma,\lambda$) + FlexZBoost (1+b)'

GammaCov = [analysis_GammaCov1, analysis_GammaCov2, analysis_GammaCov3, analysis_GammaCov2bis, analysis_GammaCov3bis]
analysis['GammaLambda'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in GammaCov]

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

