import numpy as np
import copy, glob
import _add_name_save
analysis = {}

lensing_data_truez = '../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl'
lensing_data_bpz = '../../data/stacked_esd_profiles_redmapper_BPZ_full_coverage.pkl'
lensing_data_flex = '../../data/stacked_esd_profiles_redmapper_flex_full_coverage.pkl'
#impact_c_m
analysis_Duffy08 = {
              'type': 'WLxN',
              'fit_cosmo':False,
              'density_profile':'nfw',
              'cM_relation':'Duffy08',
              'two_halo':False,
              'hmf':'Despali16',
              'radius_max':5.5,
              'radius_min':1.,
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
analysis['c_M'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in cM]

#two-halo term
analysis_2h = analysis_Diemer15.copy()
analysis_2h['two_halo'] = True
analysis_2h['radius_max'] = 15
analysis_2h['name_plot'] = 'nfw - Diemer15 (2h)'
twoh = [analysis_2h]
analysis['2h']  = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in twoh]

#hmf
analysis_Bocquet16 = analysis_Diemer15.copy()
analysis_Bocquet16['hmf'] = 'Bocquet16'
analysis_Bocquet16['name_plot'] = 'nfw - Bocquet16'
analysis_hmf = [analysis_Bocquet16]
analysis['hmf'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in analysis_hmf]

#photoz
analysis_BPZ = analysis_Diemer15.copy()
analysis_BPZ['lensing_data'] = lensing_data_bpz
analysis_BPZ['name_plot'] = 'BPZ'
analysis_BPZ['photoz'] = 'BPZ'
analysis_flex = analysis_Diemer15.copy()
analysis_flex['lensing_data'] = lensing_data_flex
analysis_flex['name_plot'] = 'FleXZBoost'
analysis_flex['photoz'] = 'flex'
photoz = [analysis_BPZ, analysis_flex]
analysis['photoz'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in photoz]

#GammaCov covariance
analysis_GammaCov1 = analysis_Diemer15.copy()
analysis_GammaCov1['shear_richness_cov'] = True
analysis_GammaCov1['name_plot'] = r'nfw + Cov($\Delta\Sigma,\lambda$)'
GammaCov = [analysis_GammaCov1]
analysis['GammaLambda'] = [_add_name_save.add_name_save_cluster_scaling_relation(a) for a in GammaCov]

print('--- WLxN ---')
analysis_list = []
for ff, k in enumerate(analysis.keys()):
    analysis_list += analysis[k]
analysis_WLxN = analysis_list

print('--- N ---')
analysis_N = copy.deepcopy(analysis_WLxN)
for analysis_ in analysis_N:
    analysis_['type'] = 'N'
    _add_name_save.add_name_save_cluster_scaling_relation(analysis_)

print('--- WL ---')
analysis_WL = copy.deepcopy(analysis_WLxN)
for analysis_ in analysis_WL:
    analysis_['type'] = 'WL'
    _add_name_save.add_name_save_cluster_scaling_relation(analysis_)

print('--- MxN ---')
analysis_MxN = copy.deepcopy(analysis_WLxN)
for analysis_ in analysis_MxN:
    analysis_['type'] = 'MxN'
    name_cl_mass, name_save_cl_mass = _add_name_save.add_name_save_cluster_mass_measurement(analysis_)
    file_mass = glob.glob('../../cluster_mass_measurement/*')
    file_mass = [f.split('/cluster_mass_measurement/')[1].split('.pkl')[0] for f in file_mass]
    if name_cl_mass not in file_mass:
        print(name_cl_mass)
    analysis_['mass_file'] = name_save_cl_mass
    _add_name_save.add_name_save_cluster_scaling_relation(analysis_)

print('--- M ---')
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
               