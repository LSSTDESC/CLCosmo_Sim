import numpy as np
from astropy.table import Table
import random
import pickle
from scipy.integrate import quad

def load(filename, **kwargs):
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)

def binning(corner): return [[corner[i],corner[i+1]] for i in range(len(corner)-1)]

def load_profile(profile_name = None, r_in = 'r', gt_in = 'gt',
                 gx_in = 'gx', weight = 'norm_sum', rmin = 1, rmax = 5):
    r"""
    Attributes:
    -----------
    load individual profiles on a given radial range
    profile_name: Table
        individual profiles
    r_in: str
        name of radius colomn (in)
    gt_in:str
        name of gt column (in)
    gx_in:str
        name of gx column (in)
    weight:str
        name of weight column
    rmin:float
        minimum radius
    rmax:float
        maximum radius
    """
    profile = load(profile_name)
    gt_av_cut,gx_av_cut, r_av_cut, norm_sum_cut = [], [], [], []
    mask = (profile[r_in][0] > rmin) * (profile[r_in][0] < rmax)
    for p in profile:
        gt_av_cut.append(p[gt_in][mask])
        gx_av_cut.append(p[gx_in][mask])
        r_av_cut.append(p[r_in][mask])
        norm_sum_cut.append(p[weight][mask])
    profile[gt_in] = np.array(gt_av_cut)
    profile[gx_in] = np.array(gx_av_cut)
    profile[r_in] = np.array(r_av_cut)
    profile[weight] = np.array(norm_sum_cut)
    return profile

def mean_value(profile = None, r_in = 'r',gt_in = 'gt', gx_in = 'gx',r_out = 'r',
               gt_out = 'gt',gx_out = 'gx',weight = 'W_l'):
    r"""
    compute stacked profile from individual profiles
    Attributes:
    -----------
    profile_name: Table
        individual profiles
    r_in: str
        name of radius colomn (in)
    gt_in:str
        name of gt column (in)
    gx_in:str
        name of gx column (in)
    r_out: str
        name of radius column (out)
    gt_out: str
        name of gt column (out)
     gx_out: str
        name of gx column (out)
    weight:str
        name of weight column
    """
    gt_w=np.average(profile[gt_in],weights=profile[weight],axis=0)
    gx_w=np.average(profile[gx_in],weights=profile[weight],axis=0)
    r_w=np.average(profile[r_in],weights=None,axis=0)
    return gt_w, gx_w, r_w

def stacked_profile(profile = None,r_in = '1',gt_in = '1', gx_in = '1',
                    r_out = '1',gt_out = '1', gx_out = '1',
                    weight = '1',
                    z_name = '1', z_err = '1',
                    obs_name = '1', obs_err = '1',
                    Z_bin = None, Obs_bin = None, add_columns_to_bin = []):
    r"""
    compute stacked profiles from individual profiles in bins of redshift 
    and a given observable
    Attributes:
    -----------
    profile_name: Table
        individual profiles
    r_in: str
        name of radius colomn (in)
    gt_in:str
        name of gt column (in)
    gx_in:str
        name of gx column (in)
    r_out: str
        name of radius column (out)
    gt_out: str
        name of gt column (out)
     gx_out: str
        name of gx column (out)
    weight:str
        name of weight column
    z_name: str
        name of cluster redshift
    z_err: str
        name of redshift error
    obs_name: str
        name of observable
    obs_err: str
        name of err observable
    Z_bin: list
        list of redshift bins
    Obs_bin: list
        list of observable bins
    add_columns_to_bin: list
        if provided, bin columns
    """
    colname = ['z_mean','obs_mean','obs_rms',
               'radius','gt','gx', 
               'gt_individual', 'radius_individual',
               'n_stack','cluster_id',
               'z_individual', 'obs_individual',
               'z_bin', 'obs_bin']
    
    colname = colname + add_columns_to_bin
    data = {name : [] for name in colname}
    for z_bin in Z_bin:
        condition_z = (profile[z_name] < z_bin[1])*(profile[z_name] > z_bin[0])
        for obs_bin in Obs_bin:
            condition = condition_z * (profile[obs_name] < obs_bin[1]) * (profile[obs_name] > obs_bin[0])
            if condition.shape == (len(profile), 1):
                condition = [c[0] for c in condition]
            p = profile[condition]
            print(len(p))
            if len(p) == 0: continue
            obs_mean, obs_rms = np.average(p[obs_name]), np.std(p[obs_name])/np.sqrt(len(p))
            z_mean = np.average(p[z_name])
            gt, gx, r = mean_value(profile = p, r_in = r_in,gt_in = gt_in, gx_in = gx_in,
               r_out = r_out,gt_out = gt_out,gx_out = gx_out,weight = weight)
            n = len(p)
            gt_individual = p[gt_in]
            radius_individual = p[r_in]
            array = [z_mean, obs_mean, obs_rms, 
                     r, gt, gx, gt_individual, radius_individual, 
                     n, p['cluster_id'], p[z_name], p[obs_name], z_bin, obs_bin]
            array = array + [p[name_add] for name_add in add_columns_to_bin]
            for i, name in enumerate(colname):
                data[name].append(array[i])
    data = Table(data)
    return data

def av_sigma_2(zl, cosmo=None):
    def sigma_c_2(zs):
        return cosmo.eval_sigma_crit(zl, zs)**2.
    def Chang(z): #Chang distribution of source redshift
        a = 1.24
        b = 1.01
        z_0 = 0.51
        return np.exp(-(z/z_0)**b)*z**a
    norm = quad(Chang,0,np.inf)[0]
    def __INTEGRAND__(z): return Chang(z)*sigma_c_2(z)
    return quad(__INTEGRAND__,zl + 0.1,100)[0]/norm

def shape_noise_ds(profile = 1,z_name = 'z_mean',r_corner=None, shapenoise=.1, ns_arcmin2=None, cosmo=None):
    
    #def ns_arcmin2(z):
    #    dy = 35-10
    #    dx=0.9-0.3
    #    return (z-0.3) * dy/dx + 35
    res=[]
    for p in profile:
        dS_arcmin2 = 3437.75**2*np.pi*np.array([r_corner[i+1]**2 - r_corner[i]**2 for i in range(len(r_corner)-1)])/(cosmo.eval_da(p[z_name])**2)
        Ngal = dS_arcmin2*ns_arcmin2#ns_arcmin2(p['z_mean'])*dS_arcmin2
        av_s_2 = av_sigma_2(p[z_name], cosmo=cosmo)
        res.append(av_s_2*shapenoise**2/(Ngal*p['n_stack']))
    return np.array(res)

    

def bootstrap_covariance(profile = 1,
                    r_in = 'radius',
                    gt_in = 'gt', gx_in = 'gx',
                    r_out = 'radius',
                    gt_out = 'gt', gx_out = 'gx',
                    weight = 'W_l',
                    z_name = 'z_mean', obs_name = 'obs_mean',
                    n_boot = 100,
                    Z_bin = None, Obs_bin = None):
    
    colname = ['z_mean','obs_mean','obs_rms', 'cov_t','inv_cov_t', 'cov_x', 'gt_boot', 'gx_boot', 'gt_err', 'gx_err', 'gt_random']
    data = {name : [] for name in colname}
    for z_bin in Z_bin:
        condition_z = (profile[z_name] < z_bin[1])*(profile[z_name] > z_bin[0])
        for obs_bin in Obs_bin:
            condition = condition_z*(profile[obs_name] < obs_bin[1])*(profile[obs_name] > obs_bin[0])
            if condition.shape == (len(profile), 1):
                condition = [c[0] for c in condition]
            p = profile[condition]
            if len(p) == 0: continue
            obs_mean, obs_rms = np.mean(p[obs_name]), np.std(p[obs_name])
            z_mean = np.mean(p[z_name])
            n_cluster=len(p)
            gt, gx = [], []
            n_boot_used = int(n_boot)#*( 1 + 1000/n_cluster ))
            for n in range(n_boot_used):
                p_boot = p[np.random.choice(np.arange(n_cluster), n_cluster)]
                gt_boot, gx_boot, r_boot = mean_value(profile = p_boot, 
                                       r_in = r_in, gt_in = gt_in, gx_in = gx_in,
                                       r_out = r_out, gt_out = gt_out, gx_out = gx_out, weight = weight)
                gt.append(np.array(gt_boot)), gx.append(np.array(gx_boot))

            gt, gx = np.array(gt), np.array(gx)
            Xt, Xx = np.stack((gt.astype(float)), axis = 1), np.stack((gx.astype(float)), axis = 1)
            cov_t, cov_x = np.cov(Xt, bias=False), np.cov(Xx)
            array = [z_mean, obs_mean, obs_rms, cov_t, np.linalg.inv(cov_t),
                     cov_x, gt, gx, np.sqrt(cov_t.diagonal()), np.sqrt(cov_x.diagonal()), gt]
            gt, gx=0,0
            for i, name in enumerate(colname):
                data[name].append(array[i])
                 
    data = Table(data)
    
    return data

def jacknife_covariance(profile = 1,
                    r_in = '1',
                    gt_in = '1', gx_in = '1',
                    r_out = '1',
                    gt_out = '1', gx_out = '1',
                    weight = '1',
                    z_name = '1', obs_name = '1',
                    n_jack = 1,
                    ra = 'ra', dec = 'dec',
                    Z_bin = 1, Obs_bin = 1):
    
    n_jack_ra = round(np.sqrt(n_jack))
    n_jack_dec = n_jack_ra
    ra_max, ra_min = np.max(profile[ra]), np.min(profile[ra])
    dec_max, dec_min = np.max(profile[dec]), np.min(profile[dec])
    ra_corner = np.linspace(ra_min, ra_max, n_jack_ra + 1)
    dec_corner = np.linspace(dec_min, dec_max, n_jack_dec + 1)
    Ra_bin  = binning(ra_corner)
    Dec_bin  = binning(dec_corner)
    colname = ['z_mean','obs_mean','obs_rms', 'cov_t', 'cov_x', 'gt_boot', 'gx_boot', 'gt_err', 'gx_err', 'Hartlap']
    data = {name : [] for name in colname}
    for z_bin in Z_bin:
        
        condition_z = (profile[z_name] < z_bin[1])*(profile[z_name] > z_bin[0])
        
        for obs_bin in Obs_bin:
            
            condition = condition_z*(profile[obs_name] < obs_bin[1])*(profile[obs_name] > obs_bin[0])
            
            condition = [c[0] for c in condition]
            
            p = profile[condition]
            
            if len(p) == 0: continue
                
            obs_mean, obs_rms = np.mean(p[obs_name]), np.std(p[obs_name])
            z_mean = np.mean(p[z_name])
            gt_JK, gx_JK = [], []

            for ra_bin in Ra_bin:
                
                for dec_bin in Dec_bin:
                    mask_jacknife = (p[ra] > ra_bin[0])*(p[ra] < ra_bin[1])*(p[dec] > dec_bin[0])*(p[dec] < dec_bin[1])
                    profile_jacknife = p[np.invert(mask_jacknife)]
                    gt_jk, gx_jk, r_jk = mean_value(profile = profile_jacknife, 
                                           r_in = r_in,
                                           gt_in = gt_in, 
                                           gx_in = gx_in,
                                           r_out = r_out,
                                           gt_out = gt_out,
                                           gx_out = gx_out,
                                           weight = weight)
                    gt_JK.append(np.array(gt_jk))
                    gx_JK.append(np.array(gx_jk))
            gt, gx = np.array(gt_JK), np.array(gx_JK)
            Xt, Xx = np.stack((gt.astype(float)), axis = 1), np.stack((gx.astype(float)), axis = 1)
            cov_t, cov_x = np.cov(Xt, bias = False), np.cov(Xx, bias = False)
            cov_t, cov_x = ((n_jack-1)**2/n_jack)*cov_t, ((n_jack-1)**2/n_jack)*cov_x
            H = (n_jack - cov_t.shape[0] - 2)/(n_jack - 1)
            array = [z_mean, obs_mean, obs_rms, cov_t, 
                     cov_x, gt, gx, np.sqrt(cov_t.diagonal()), np.sqrt(cov_x.diagonal()), H]
            
            for i, name in enumerate(colname):
                data[name].append(array[i])
                 
    data = Table(data)
    
    return data

def sample_covariance(profile = 1,
                    r_in = '1',
                    gt_in = '1', gx_in = '1',
                    r_out = '1',
                    gt_out = '1', gx_out = '1',
                    weight = '1',
                    z_name = '1', obs_name = '1',
                    Z_bin = 1, Obs_bin = 1):

    colname = ['z_mean','obs_mean','obs_rms', 'cov_t', 'cov_x', 'gt_boot', 'gx_boot', 'gt_err', 'gx_err']
    
    data = {name : [] for name in colname}
    
    for z_bin in Z_bin:
        
        condition_z = (profile[z_name] < z_bin[1])*(profile[z_name] > z_bin[0])
        
        for obs_bin in Obs_bin:
            
            condition = condition_z*(profile[obs_name] < obs_bin[1])*(profile[obs_name] > obs_bin[0])
            
            condition = [c[0] for c in condition]
            
            p = profile[condition]
            
            if len(p) == 0: continue
                
            obs_mean, obs_rms = np.mean(p[obs_name]), np.std(p[obs_name])
            z_mean = np.mean(p[z_name])
            gt_individual, gx_individual, w_unit = [], [], []
            
            for i, p_individual in enumerate(p):
                
                gt_individual.append(np.array(p_individual[gt_in]))
                gx_individual.append(np.array(p_individual[gx_in]))
                w_unit.append(np.array([1 if w_ > 0 else 0 for w_ in p_individual[weight]]))
                
            p['w_unit'] = np.array(w_unit)
                
            gt_mean = np.sum(p[gt_in]*p['w_unit'], axis = 0)/np.sum(p['w_unit'], axis = 0)
            
            cov_diag_sample = np.sum((p[gt_in] - gt_mean)**2*p['w_unit'], axis = 0)/(np.sum(p['w_unit'], axis = 0) - 1)
            
            cov_diag_mean = cov_diag_sample/len(p)
            
            cov_t = np.zeros([len(cov_diag_mean),len(cov_diag_mean)])
            cov_x = np.zeros([len(cov_diag_mean),len(cov_diag_mean)])
            for i in range(len(cov_diag_mean)):
                cov_t[i,i] = cov_diag_mean[i]

            #gt, gx, w = np.array(gt_individual), np.array(gx_individual), np.array(weights_individual)
            #Xt, Xx = np.stack((gt.astype(float)), axis = 1), np.stack((gx.astype(float)), axis = 1)
            #cov_t, cov_x = np.cov(Xt, aweights =  None)/len(p), np.cov(Xx)/len(p)
            
            array = [z_mean, obs_mean, obs_rms, cov_t, cov_diag_mean, 1, 1, cov_diag_sample, cov_diag_sample]
            
            for i, name in enumerate(colname):
                data[name].append(array[i])
                 
    data = Table(data)
    
    return data

