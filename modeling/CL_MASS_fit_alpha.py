import numpy as np
from scipy.optimize import curve_fit

def numerical_hessian(func, p, eps=1e-4):
    """
    Compute numerical Hessian matrix of func at point p.
    Central differences for better accuracy.
    """
    n = len(p)
    hessian = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            p_i1 = p.copy()
            p_i2 = p.copy()
            p_j1 = p.copy()
            p_j2 = p.copy()

            # Move parameter i
            p_i1[i] += eps
            p_i2[i] -= eps

            # Move parameter j
            p_j1[j] += eps
            p_j2[j] -= eps

            if i == j:
                # Second derivative wrt same parameter
                fpp = func(p_i1)
                fmm = func(p_i2)
                f0  = func(p)
                hessian[i, i] = (fpp - 2*f0 + fmm) / eps**2

            else:
                # Mixed derivative wrt i and j
                pp = p.copy()
                pm = p.copy()
                mp = p.copy()
                mm = p.copy()

                pp[i] += eps; pp[j] += eps
                pm[i] += eps; pm[j] -= eps
                mp[i] -= eps; mp[j] += eps
                mm[i] -= eps; mm[j] -= eps

                fpp = func(pp)
                fpm = func(pm)
                fmp = func(mp)
                fmm = func(mm)

                hessian[i, j] = (fpp - fpm - fmp + fmm) / (4 * eps**2)
                hessian[j, i] = hessian[i, j]

    return hessian

def store_log10M(Z_bin, Obs_bin, param_list, name_k, z_name='z_mean', 
                                                     richness_name='obs_mean',
                                                     log10m_name='log10M200c_WL',
                                                     log10m_err_name='err_log10M200c_WL'):
    
    Log10M = np.zeros((len(Obs_bin), len(Z_bin), len(param_list)))
    Log10M_err = np.zeros((len(Obs_bin), len(Z_bin), len(param_list)))
    
    for i, richness_bin in enumerate(Obs_bin):
        
        for j, z_bin in enumerate(Z_bin):
        
            log10m, log10m_err = [], []
            
            for k, param_ in enumerate(param_list):
                
               mass_file = np.load(name_k.format(k), allow_pickle=True)['masses']
               mask = (np.array(mass_file[z_name]) > z_bin[0])
               mask *= (np.array(mass_file[z_name]) < z_bin[1])
               mask *= (np.array(mass_file[richness_name]) > richness_bin[0])
               mask *= (np.array(mass_file[richness_name]) < richness_bin[1])
               log10m.append(mass_file[mask][log10m_name][0])
               log10m_err.append(mass_file[mask][log10m_err_name][0])
    
            Log10M[i,j,:] = np.array(log10m)
            Log10M_err[i,j,:] = np.array(log10m_err)
    
    return Log10M, Log10M_err

def fit_linear_model(log10m_array, Om_array, sigma_log10m_array, Om_fid, Om_range):
    # Apply mask to select data within the given Om_range
    mask = (Om_array > Om_range[0]) & (Om_array < Om_range[1])

    log10m_array = np.asarray(log10m_array)[mask]
    Om_array = np.asarray(Om_array)[mask]
    sigma_log10m_array = np.asarray(sigma_log10m_array)[mask]

    # Define the linear model: log10m = log10m0 + alpha * (Om - Om_fid)
    def model(Om, log10m0, alpha):
        return log10m0 + alpha * (Om - Om_fid)

    # Fit the model to the data
    popt, pcov = curve_fit(
        model,
        Om_array,
        log10m_array,
        sigma=sigma_log10m_array,
        absolute_sigma=True)

    log10m0, alpha = popt
    log10m0_err, alpha_err = np.sqrt(np.diag(pcov))

    return log10m0, alpha, log10m0_err, alpha_err

def fit_Mass_Model(Log10M, Log10M_err, Z_bin, Obs_bin, 
                   Om_list, Om_fid, 
                   Om_range = [0.1, 0.6], merge_redshift=False):

    log10M0 = np.zeros(Log10M[:,:,0].shape)
    log10M0_err = np.zeros(Log10M[:,:,0].shape)

    Alpha = np.zeros(Log10M[:,:,0].shape)
    Alpha_err = np.zeros(Log10M[:,:,0].shape)

    sigma_log10M0 = np.zeros(Log10M[:,:,0].shape)
    sigma_log10M0_err = np.zeros(Log10M[:,:,0].shape)

    sigma_Alpha = np.zeros(Log10M[:,:,0].shape)
    sigma_Alpha_err = np.zeros(Log10M[:,:,0].shape)

    for i, richness_bin in enumerate(Obs_bin):

        if merge_redshift == False:
            
            for j, z_bin in enumerate(Z_bin):
                
                log10m_array, sigma_log10m_array = Log10M[i,j,:], Log10M_err[i,j,:]
    
                log10m0, alpha, log10m0_err, alpha_err = fit_linear_model(log10m_array, Om_list, sigma_log10m_array, 
                                                                          Om_fid, Om_range)
                
                sigma_log10m0, sigma_alpha, sigma_log10m0_err, sigma_alpha_err = fit_linear_model(sigma_log10m_array, Om_list,
                                                                                            0.02*np.ones(len(sigma_log10m_array)), 
                                                                                              Om_fid, Om_range)
                Alpha[i,j] = alpha
                Alpha_err[i,j] = alpha_err
                log10M0[i,j] = log10m0
                log10M0_err[i,j] = log10m0_err
    
                sigma_Alpha[i,j] = sigma_alpha
                sigma_Alpha_err[i,j] = sigma_alpha_err
                sigma_log10M0[i,j] = sigma_log10m0
                sigma_log10M0_err[i,j] = sigma_log10m0_err

        else: 

            from scipy.optimize import minimize
            
            def model_chi2(p, Om, Om_fid, Log10M, Log10M_err):

                alpha = p[0]
                log10M0_array = p[1:]
                res = 0.0  
            
                for j, z_bin in enumerate(Z_bin):
            
                    log10m0_ = log10M0_array[j]
                    ymodel = log10m0_ + alpha * (Om - Om_fid)
                    ydata = Log10M[i, j, :]
                    yerr = Log10M_err[i, j, :]
                    res += 0.5*np.sum(((ymodel - ydata) / yerr)**2)
            
                return res

            p0 = np.array([0.0] + [14.2 for i in range(len(Z_bin))])

            result_mean = minimize(model_chi2,
                              x0=p0,
                              args=(Om_list, Om_fid, Log10M, Log10M_err),
                              method='L-BFGS-B')

            p0 = np.array([0.0] + [0.04 for i in range(len(Z_bin))])

            result_sigma = minimize(model_chi2,
                              x0=p0,
                              args=(Om_list, Om_fid, Log10M_err, 0.2*Log10M_err),
                              method='L-BFGS-B')

            alpha_fit   = result_mean.x[0]
            log10m0_fit = result_mean.x[1:]
            def f_local(p):
                return model_chi2(p, Om_list, Om_fid, Log10M, Log10M_err)

            # compute Hessian
            H = numerical_hessian(f_local, result_mean.x, eps=1e-4)
            Cov = np.linalg.inv(H)
            err_sigma = np.sqrt(np.diag(Cov))
            
            Alpha[i,:] = alpha_fit * np.ones(len(Z_bin))
            Alpha_err[i,:] = err_sigma[0] * np.ones(len(Z_bin))
            log10M0[i,:] = log10m0_fit * np.ones(len(Z_bin))
            log10M0_err[i,:] = err_sigma[1:] * np.ones(len(Z_bin))

            sigma_alpha_fit   = result_sigma.x[0]
            sigma_log10m0_fit = result_sigma.x[1:]
            cov_sigma = np.array(result_sigma.hess_inv.todense())
            err_sigma = cov_sigma.diagonal()**.5

            sigma_Alpha[i,:] = sigma_alpha_fit * np.ones(len(Z_bin))
            sigma_Alpha_err[i,:] = err_sigma[0] * np.ones(len(Z_bin))
            sigma_log10M0[i,:] = sigma_log10m0_fit * np.ones(len(Z_bin))
            sigma_log10M0_err[i,:] = err_sigma[1:] * np.ones(len(Z_bin))
                
    return log10M0, Alpha, log10M0_err, Alpha_err, sigma_log10M0, sigma_Alpha, sigma_log10M0_err, sigma_Alpha_err