# CLCosmo_Sim repository

Authors: C. Payerne, Z. Zhang, C. Combet, M. Aguena, T. Guillemin, M. Ricci, S. Vitenti

Cluster cosmology analysis using DESC tools and simulations. This repository presents the codes that are used in the DESC project "Impact of modeling and observational systematics from cluster lensing and abundance"

# Modeling
In the /modeling directory
## Modeling cluster abundance and cluster lensing
- `CL_COUNT_cluster_abundance`: provides prediction for cluster abundance in the proxy-redshift space (binned) with observable-mass relation, completeness and purity parametrisation.
$$
N_{ij}=  \Omega \int_{z_{i}}^{z_{i+1}} dz \int_{\lambda_{j}}^{\lambda_{j+1}} d\lambda\int_{m_{\rm min}}^{+\infty}dm\ \frac{dn(m, z)}{dm}\frac{d^2V(z)}{dzd\Omega} \Phi(\lambda, m, z) P(\lambda|m,z)
$$
- `CL_LENSING_cluster_lensing`: provides prediction for cluster lensing signal in the proxy-redshift space (binned) with observable-mass relation, completeness and purity parametrisation.
- `CL_MASS_cluster_lensing`: provides prediction for cluster mean mass in the proxy-redshift space (binned) with observable-mass relation, completeness and purity parametrisation.
- `CL_COUNT_modeling_halo_mass_function`: modeling the halo mass function and cosmoving volume (CCL as backend).
- `CL_COUNT_modeling_purity`: modeling purity.
- `CL_COUNT_modeling_completeness`: modeling completeness.
- `CL_COUNT_modeling_richness_mass_relation`: modeling the observable-mass relation (mean, disperion, probability density function).

## Bayesian inference pipeline
- `CL_COUNT_class_likelihood`: provides binned Gaussian, Poissonian and unbinned Poissonian likelihoods for cluster count cosmology.
- `CL_COUNT_DATAOPS_cluster_abundance_covariance`: provides estimator of cluster abundance covariance based on data (still under construction).
- `STAT_forecast`: module for Fisher forecast.

# Data extraction and lensing signal estimation
in the /extract_data directory
The python files in this directory are dedicated to the background source extraction from cosmoDC2 using GCRCatalogs and Qserv, as well as computing cluster lensing individual lensing profiles. 
## cosmoDC2: Extraction of galaxy catalogs behind redMaPPer clusters
- `_config_extract_sources_in_cosmoDC2.py`: configuration file for source extraction in cosmoDC2
- `_utils_extract_sources_in_cosmoDC2.py` : functions to extract galaxy data using GCRcatalogs in the cosmoDC2 photoz-addons, and "truth" quantities using Qserv queries
- `_utils_photometric_redshifts.py` : functions. From photometric probability density functions, compute mean redshift, mean critical surface mass density for each galaxy
- `run_extract_sources_in_cosmoDC2.py` : module to extract background galaxy catalogs for each redMaPPer cluster, with photozs, each one saved in a pickle file. Individual catalogs can also not be saved, and lensing profiles directly computed. 
## Estimate cluster lensing profiles from background source galaxies
- `_config_lensing_profiles.py`: configuration file for estimation of lensing profiles
- `_utils_lensing_profiles.py`: set of functions for computing lensing profiles
- `run_lensing_profiles.py` : compute individual lensing profiles 

