def add_name_save_cluster_mass_measurement(analysis):
    name = 'cluster-masses_'
    add_lensing_infos = '1-halo=' + analysis['density_profile']
    if analysis['cM_relation'] != None: add_lensing_infos = add_lensing_infos + '+c-M=' + analysis['cM_relation']
    else: add_lensing_infos = add_lensing_infos + '+c-M=None' 
    if analysis['two_halo']: add_lensing_infos = add_lensing_infos + '_' + '+2-halo'
    add_lensing_infos = add_lensing_infos + '_photoz='+analysis['photoz']
    name = name + add_lensing_infos
    name_save= '../../cluster_mass_measurement/' + name + '.pkl'
    return name, name_save

def add_name_save_cluster_scaling_relation(analysis):
    #kind of fit
    name = 'MCMC_'
    if analysis['fit_cosmo']: fit_ = 'fit_m-r+cosmo_'
    else: fit_ = 'fit_m-r_'
    
    #from which observable
    name = name+ fit_
    name = name + analysis['type'] + '_'
    
    #if lensing
    if (analysis['type'] == 'WL') or (analysis['type'] == 'WLxN') or (analysis['type'] == 'MxN') or (analysis['type'] == 'M'):
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