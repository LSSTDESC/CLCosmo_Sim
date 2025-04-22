# CLCosmo_Sim repository

Authors: Payerne, Constantin (constantin.payerne@gmail.com) ;  Zhang, Zhuowen ;  Aguena, Michel ;  Combet, Céline ;  Guillemin, Thibault ;  Ricci, Marina ; Amouroux, Nathan ;  Avestruz, Camille ;  Barroso, Eduardo J. ;  Farahi, Arya ;  Kovacs, Eve ;  Murray, Calum ; Rau, Markus M. ;  Rykoff, Eli S. ;  Schmidt, Samuel J. and  the LSST Dark Energy Science Collaboration

This repository is dedicated to the analysis of the redMaPPer mass-richness relation in the LSST DESC DC2 simulations from a cluster weak gravitational lensing and abundance perspective. The work is presented in the paper ["Analysis of the weak lensing mass-richness relation of redMaPPer clusters in the LSST DESC DC2 simulations"](https://ui.adsabs.harvard.edu/abs/2025arXiv250208444P/abstract). 

We detail below the main features of the LSST DESC CLCosmo_Sim code.

This repository is associated to the [DESC project 380](https://portal.lsstdesc.org/DESCPub/app/PB/show_project?pid=380).

# Recquirements
CLCosmo_Sim has the following dependencies:

- [NumPy](https://www.numpy.org/) 
- [SciPy](https://scipy.org/) 
- [Astropy](https://www.astropy.org/) 
- [Matplotlib](https://matplotlib.org/) 
- [emcee](https://emcee.readthedocs.io/en/stable/)
- [cluster-toolkit](https://cluster-toolkit.readthedocs.io/en/latest/)
- [CCL](https://ccl.readthedocs.io/en/latest/) 
- [PySSC](https://pyssc.readthedocs.io/en/latest/)
- [NumCosmo](https://numcosmo.github.io/) 

# Modeling
In the /modeling directory
## Modeling cluster abundance and cluster lensing

- `CL_COUNT_modeling_halo_mass_function`: modeling of cosmological functions, i.e. the halo mass function, halo bias and comoving volume. We use the [Core Cosmology Library](https://ccl.readthedocs.io/en/latest/) by [Chisari et al. (2019)](https://arxiv.org/abs/1812.05995).

- `CL_COUNT_cluster_abundance`: provides prediction for cluster abundance in the proxy-redshift space (binned) with observable-mass relation, completeness and purity parametrisation. The count in the $i$-th redshift bin and in the $j$-th richness bin is given by

$$
N_{ij}=  \Omega \int_{z_{i}}^{z_{i+1}} dz \int_{\lambda_{j}}^{\lambda_{j+1}} d\lambda \int_{m_{\rm min}}^{+\infty}dm \frac{dn(m, z)}{dm}\frac{d^2V(z)}{dzd\Omega} \Phi(\lambda, m, z) P(\lambda|m,z)
$$

- `CL_COUNT_DATAOPS_cluster_abundance_covariance`: provides prediction for cluster abundance covariance in the proxy-redshift space (binned). The covariance accounts for Poisson shot noise and Super-Sample Covariance (SSC). To compute SSC, we use [PySSC](https://pyssc.readthedocs.io/en/latest/) by [Lacasa et al. (2018)](https://www.aanda.org/articles/aa/full_html/2018/03/aa30281-16/aa30281-16.html).

$$
\Sigma_{\mathrm{N}}[ij,kl] = N_{ij}\delta^K_{ik}\delta^K_{jl}+ N_{ij}N_{kl}\langle b\rangle_{ij}\langle b\rangle_{kl} S_{jl}
$$

- `CL_LENSING_cluster_lensing`: provides prediction for cluster lensing signal in the proxy-redshift space (binned) with observable-mass relation, completeness and purity parametrisation. The model is given by

$$
\Delta\Sigma_{ij}(R) = \frac{1}{N_{ij}}\int_{z_{i}}^{z_{i+1}} dz \int_{\lambda_{j}}^{\lambda_{j+1}} d\lambda \int_{m_{\rm min}}^{+\infty}dm \frac{d^2N(m, z)}{dzdm} c(m,z) P(\lambda|m,z)\Delta\Sigma(R|m,z)
$$

- `CL_MASS_cluster_lensing`: provides prediction for cluster mean mass in the proxy-redshift space (binned) with observable-mass relation, completeness and purity parametrisation.

$$
M_{ij} = \frac{1}{N_{ij}}\int_{z_{i}}^{z_{i+1}} dz \int_{\lambda_{j}}^{\lambda_{j+1}} d\lambda \int_{m_{\rm min}}^{+\infty}dm \frac{d^2N(m, z)}{dzdm} c(m,z) P(\lambda|m,z)m.
$$

- `CL_COUNT_modeling_purity` and `CL_COUNT_modeling_completeness`: modeling purity and completeness, respectively. We consider the form proposed by [Aguena & Lima (2016)](https://arxiv.org/abs/1611.05468).
- `CL_COUNT_modeling_richness_mass_relation`: modeling the observable-mass relation (mean, disperion, probability density function). We consider the log-normal model

$$
P(\lambda|m,z) \propto \frac{1}{\lambda}\exp\left(-\frac{[\ln\lambda - \langle \ln \lambda|m, z\rangle]^2}{2\sigma_{\ln\lambda|m,z}^2}\right)
$$

where we use

$$
\langle \ln \lambda|m, z\rangle = \ln\lambda_0 + \mu_z\ln\left(\frac{1 + z}{1 + z_0}\right) + \mu_m\log_{10}\left(\frac{m}{m_0}\right)
$$

and

$$
\sigma_{\ln \lambda|m, z} = \sigma_{\ln\lambda_0} + \sigma_z\ln\left(\frac{1 + z}{1 + z_0}\right) + \sigma_m\log_{10}\left(\frac{m}{m_0}\right).
$$


## Bayesian inference pipeline
In this work we aim at infering the parameters of the DC2 redMaPPer mass-richness relation by drawing the posterior distribution

$$
\mathcal{P}(\theta|\mathrm{data}) = \frac{\mathcal{L}(\mathrm{data}|\theta)|\pi(\theta)}{\mathcal{L}(\mathrm{data})}
$$

- `CL_COUNT_class_likelihood`: provides binned Gaussian, Poissonian and unbinned Poissonian likelihoods for cluster count cosmology. We use a Gaussian likelihood for either the cluster counts, the stacked lensing profiles of stacked lensign masses.
- `STAT_forecast`: module for Fisher forecast.

# Data extraction
in the /extract_data directory
The python files in this directory are dedicated to the background source extraction from cosmoDC2 using GCRCatalogs and Qserv, as well as computing cluster lensing individual lensing profiles. All the data (cluster/halo catalog, lensing profiles, ...) that are used in this repository are not publicly available, and stored in the private LSSTDESC repository [CLCosmo_Sim_database](https://github.com/LSSTDESC/CLCosmo_Sim_database). Data are only available for DESC members. If you have access and want to use the data, please clone the `CLCosmo_Sim_database` repository in the same directory as `CLCosmo_Sim`.

## Extraction of redMaPPer cluster catalog
- First, `python run_extract_cluster_catalog_redMaPPer.py` extracts the catalog of redMaPPer clusters (position, richness, redshift), and their member galaxies (position, redshifts).
## cosmoDC2: Extraction of galaxy catalogs behind redMaPPer clusters
- `_config_extract_sources_in_cosmoDC2.py` is a configuration file for the cosmoDC2 source selection. The current selection is based on the photometric redshift PDFs, such as a galaxy is chosen as a source if

$$
\langle z\rangle_s > z_l + 0.2\hspace{0.5cm} \mathrm{and}\hspace{0.5cm} P(z_s > z_l) = \int_{z_{\rm cl}}^{\rm +\infty} dz_{\rm gal}p_{\rm photoz}(z_{\rm gal})> 0.8
$$

- `_utils_extract_sources_in_cosmoDC2.py` : functions to extract galaxy data using GCRcatalogs in the cosmoDC2 photoz-addons, and "truth" quantities using Qserv queries.
- `_utils_photometric_redshifts.py` : functions. From photometric probability density functions, compute mean redshift, mean critical surface mass density for each galaxy.
- `run_extract_sources_in_cosmoDC2.py` : module to extract background galaxy catalogs for each redMaPPer cluster, with photozs, each one saved in a pickle file. Individual catalogs can also not be saved, and lensing profiles directly computed.


# Estimation of stacked excess surface mass density profiles arround redMaPPer clusters
- `_config_lensing_profiles.py`: configuration file for estimation of lensing profiles
- `_utils_lensing_profiles.py`: set of functions for computing individual lensing profiles with the DESC [Cluster Lensing Mass Modeling (CLMM)](https://github.com/LSSTDESC/CLMM) code. This module is called at the level of the extraction of source sample catalog (above), to compute directy the individual lensing profiles. 
- `compute_stacked_lensing_profiles.py` : compute stacked lensing profiles in bins of redshift and richness.
