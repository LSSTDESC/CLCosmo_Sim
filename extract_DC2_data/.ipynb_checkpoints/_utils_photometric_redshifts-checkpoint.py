import numpy as np
from astropy.table import Table
from scipy.integrate import simps
import clmm.theory as m
from clmm.dataops import compute_background_probability
from scipy.integrate import simps
from scipy.interpolate import interp1d

def draw_z_from_pdf(pdf, pzbins, n_samples=1, use_clmm=False):
    r"""
    Attributes:
    -----------
    pdf: array
        tabulated photometric distribution
    pzbins: array
        tabulated redshift axis
    n_samples: int
        number of samples from the pdf
    use_clmm: Bool
        use clmm or not
    Returns:
    --------
    z_sample: array
        redshift samples from pdfs
    """
    n_pdf = len(pdf)
    z_array=pzbins[0]
    z_sample = np.zeros([n_pdf, n_samples])
    #create cdf from pdf
    cdf_unormed = np.cumsum(pdf,  axis=1)
    cdf_random = np.random.random(z_sample.shape)
    
    for i in range(n_pdf):
        
        #normalize cdf
        cdf = cdf_unormed[i,:]/max(cdf_unormed[i,:])
        #inverse method
        z_sample[i,:]=np.interp(cdf_random[i,:], cdf, z_array)
    
    return z_sample

def compute_photoz_sigmac(z_lens, pdf_norm, pzbins, cosmo=None, use_clmm=False):
    r"""
    Attributes:
    -----------
    z_lens: float
        lens redshift
    pdf: array
        photoz distrib
    pzbins: array
        z_bin centers
    Returns:
    --------
    sigma_c: array
        photometric weak lensing sigma_c
    """
    if use_clmm==True:
        sigma_c=compute_critical_surface_density(cosmo, z_lens, 
                                                 z_source=None, use_pdz=True, 
                                                 pzbins=pzbins, pzpdf=pdf_norm)
        return sigma_c
    else:
        sigmacrit_1 = cosmo.eval_sigma_crit(z_lens, pzbins[0,:])**(-1.)
        sigmacrit_1_integrand = (pdf_norm*sigmacrit_1.T)
        return simps(sigmacrit_1_integrand, pzbins[0,:], axis=1)**(-1.)

def compute_p_background(z_lens, pdf_norm, pzbins, use_clmm=False):
    r"""
    Attributes:
    -----------
    z_lens: float
        lens redshift
    pdf: array
        photoz distrib
    pzbins: array
        z_bin centers
    Returns:
    --------
    p: array
        probability background
    """

    if use_clmm==False:
        
        n_pdf = len(pdf_norm)
        z_array=pzbins[0]
        cdf_unormed = np.cumsum(pdf_norm,  axis=1)
        p_back = np.zeros(n_pdf)
        
        for i in range(n_pdf):

            #normalize cdf
            cdf = cdf_unormed[i,:]/max(cdf_unormed[i,:])
            #inverse method
            cdf_z_lens=np.interp(z_lens,  z_array, cdf,)
            p_back[i] = 1 - cdf_z_lens
            
        return p_back
    else: 
        return compute_background_probability(z_lens, z_source=None, use_pdz=True, 
                                              pzpdf=pdf_norm, pzbins=pzbins, validate_input=True)

def compute_photoz_dispersion(pdf_norm, pzbins):
    r"""
    Attributes:
    -----------
    pdf: array
        photoz distrib
    pzbins: array
        z_bin centers
    Returns:
    --------
    dispersion: array
        standard deviation of photoz distrib
    """
#     if pdf.shape!=(len(pdf), len(pzbins)):
#         pdf_new=np.zeros([len(pdf), len(pzbins[0])])
#         for i, pdf_ in enumerate(pdf):
#             pdf_new[i,:] = pdf_
#         pdf = pdf_new
#     norm=simps(pdf, pzbins, axis=1)
#     pdf_norm=(pdf.T*(1./norm)).T
    #mean
    pdf_times_z = pdf_norm * pzbins[0,:]
    mean_z = simps(pdf_times_z, pzbins[0,:], axis=1)
    #mean2
    pdf_times_z2 = pdf_norm * pzbins[0,:] ** 2
    mean_z2 = simps(pdf_times_z2, pzbins[0,:], axis=1)
    #dispersion
    return np.sqrt(mean_z2 - mean_z ** 2) 

def compute_photoz_quantities(z_lens, pdf, pzbins, n_samples_per_pdf=3, 
                              p_background=True, sigmac_eff=True, sigmac_estimate=True,
                              photoz_dispersion=True,
                              cosmo=None, use_clmm=False):
    r"""
    compute photometric redshift dependant quantities (sigmac, p_background, sampled redshifts)
    Attributes:
    -----------
    z_lens: float
        lens redshift
    pdf: array
        photoz distrib
    pzbins: array
        z_bin centers
    n_samples_per_pdf: int
        number of samples from the pdf
    cosmo: Cosmology CLMM
        CLMM cosmology
    use_clmm: Bool
        use_clmm or not
    Returns:
    --------
    data: Table
        photoz_quantities for WL
    """
    dat_to_save = Table()
    
    #pre-compute pdf norm
    if pdf.shape!=(len(pdf), len(pzbins)):
        pdf_new=np.zeros([len(pdf), len(pzbins[0])])
        for i, pdf_ in enumerate(pdf):
            pdf_new[i,:] = pdf_
        pdf = pdf_new
    norm=simps(pdf, pzbins, axis=1)
    pdf_norm=(pdf.T*(1./norm)).T
    
    if sigmac_eff==True:
        dat_to_save['sigmac_photoz'] = compute_photoz_sigmac(z_lens, pdf_norm, pzbins, cosmo=cosmo, use_clmm=use_clmm)
    
    if sigmac_estimate==True: 
        z_samples = draw_z_from_pdf(pdf_norm, pzbins, n_samples_per_pdf, use_clmm=use_clmm)
        for i in range(n_samples_per_pdf):
            dat_to_save['z_estimate_' + str(i)] = z_samples[:,i]
            dat_to_save['sigmac_estimate_' + str(i)] = cosmo.eval_sigma_crit(z_lens, z_samples[:,i])
            
    if p_background==True:
        dat_to_save['p_background'] = compute_p_background(z_lens, pdf_norm, pzbins, use_clmm=False)
    
    if photoz_dispersion==True:
        dat_to_save['photoz_dispersion'] = compute_photoz_dispersion(pdf_norm, pzbins)

    return dat_to_save
    
 #sigma_c point estimate
    #z_samples = draw_z_from_pdf(pdf_norm, pzbins, n_samples_per_pdf, use_clmm=use_clmm)
    #sigma_c_estimate = np.zeros([len(pdf_norm), n_samples_per_pdf])
    #for i in range(n_samples_per_pdf):
     #   sigma_c_estimate[:,i] = cosmo.eval_sigma_crit(z_lens, z_samples[:,i])

    #data[:,2] = err_photoz
    
   # start = len(name_base) - 1
   # for i in range(n_samples_per_pdf): 
   #     data[:,start + i + 1] = sigma_c_estimate[:,i]
   #     data[:,start + n_samples_per_pdf + i + 1] = z_samples[:,i]
        #label
#     name = []
#     if p_background==True: name = name + ['p_background']
#     if sigmac_eff==True: name = name + ['sigmac_photoz']
#     if sicmac_estimate==True:
#         name = name + ['sigmac_photoz_estimate_' + str(k) for k in range(n_samples_per_pdf)]
#         name = name + ['z_estimate_' + str(k) for k in range(n_samples_per_pdf)]
#     if photoz_dispersion==True: name = name+['photoz_dispersion']
#def cdf_from_pdf(pdf, z_array):
#     r"""
#     Attributes:
#     -----------
#     pdf: array
#         tabulated photometric distribution
#     z_array: array
#         tabulated redshift axis
#     Returns:
#     --------
#     cdf: array
#         tabulated photometric cumulative distribution
#     """
#     cdf = np.array([simps(pdf[np.arange(j+1)], 
#                     z_array[np.arange(j+1)]) 
#                     for j in range(len(z_array))])
#     return cdf/max(cdf)

# def inverse_cdf_fct(cdf, z_array):
#     r"""
#     Attributes:
#     -----------
#     cdf: array
#         tabulated photometric cumulative distribution
#     z_array: array
#         tabulated redshift axis
#     Returns:
#     --------
#     cdf_1: fct
#         inverse cumulative distribution (interpolated)
#     """
#     return interp1d(cdf, z_array, kind='linear')

# def draw_z_from_inverse_cdf(cdf_1, n_samples=1):
#     r"""
#     Attributes:
#     -----------
#     cdf_1: fct
#         inverse cumulative distribution (interpolated)
#     n_samples: int
#         number of samples from the pdf
#     Returns:
#     --------
#     z_samples: array
#         sampled redshifts
#     """
#     z_samples = cdf_1(np.random.random(n_samples) * 1)
#     return z_samples

# def draw_z_from_pdf(pdf, pzbins, n_samples=1, use_clmm=False):
#     r"""
#     Attributes:
#     -----------
#     pdf: array
#         tabulated photometric distribution
#     pzbins: array
#         tabulated redshift axis
#     n_samples: int
#         number of samples from the pdf
#     use_clmm: Bool
#         use clmm or not
#     Returns:
#     --------
#     z_sample: array
#         redshift samples from pdfs
#     """
#     n_pdf = len(pdf)
#     z_array=pzbins[0]
#     z_sample = np.zeros([n_pdf, n_samples])
#     for i in range(n_pdf):
#         #create cdf from pdf
#         cdf = cdf_from_pdf(pdf[i], z_array)
#         #create inverse cdf from cdf
#         cdf_1 = inverse_cdf_fct(cdf, z_array)
#         #draw sample from cdf_1 (inverse method)
#         z_sample[i,:]=np.array(draw_z_from_inverse_cdf(cdf_1, 
#                                                        n_samples=n_samples))

#     return z_sample
