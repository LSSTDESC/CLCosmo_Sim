def add_name_save_cluster_mass_measurement(analysis):
    name = 'cluster-masses_'
    add_lensing_infos = '1-halo=' + analysis['density_profile']
    add_lensing_infos = add_lensing_infos + '+c-M=' + analysis['cM_relation']
    if analysis['two_halo']=='True': add_lensing_infos = add_lensing_infos + '_' + '+2-halo'
    add_lensing_infos = add_lensing_infos+ '_rmin='+str(analysis['radius_min'])+'-rmax='+str(analysis['radius_max'])
    add_lensing_infos = add_lensing_infos + '_photoz='+analysis['photoz']
    name = name + add_lensing_infos
    name_save= '../../cluster_mass_measurement/' + name + '.pkl'
    return name, name_save

def add_name_save_cluster_scaling_relation(analysis):
    #kind of fit
    name = 'MCMC_'
    if analysis['fit_cosmo']=='True': fit_ = 'fit_m-r+cosmo_'
    else: fit_ = 'fit_m-r_'
    
    #from which observable
    name = name+ fit_
    name = name + analysis['type'] 
    
    #if lensing
    if (analysis['type'] == 'WL') or (analysis['type'] == 'WLxN') or (analysis['type'] == 'MxN') or (analysis['type'] == 'M'):
        add_lensing_infos = '_1-halo=' + analysis['density_profile'] + '+c-M=' + analysis['cM_relation']
        if analysis['two_halo']=='True':
            add_lensing_infos = add_lensing_infos + '_+2-halo'
        if analysis['add_bias_lensing']=='True':
            add_lensing_infos = add_lensing_infos + '+lensing_bias'
        add_lensing_infos = add_lensing_infos+ '_rmin='+str(analysis['radius_min'])+'-rmax='+str(analysis['radius_max'])
        if analysis['shear_richness_cov']=='True':
            add_lensing_infos = add_lensing_infos + '+cov(g,richness)'
        add_lensing_infos = add_lensing_infos + '_photoz='+analysis['photoz']
        name = name + add_lensing_infos
    #kind of hmf
    if analysis['redshift_range']=='Partial': name += '_low_z_sample'
    if analysis['richness_range']=='Partial': name += '_low_richness_sample'
    name = name + '_hmf='+analysis['hmf']
    
    analysis['name'] = name
    return analysis
