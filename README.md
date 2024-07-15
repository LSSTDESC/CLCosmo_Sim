# CLCosmo_Sim repository

Authors: C. Payerne, Z. Zhang, C. Combet, M. Aguena, T. Guillemin, M. Ricci, S. Vitenti

Cluster cosmology analysis using DESC tools and simulations. This repository presents the codes that are used in the DESC project "Impact of modeling and observational systematics from cluster lensing and abundance"

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
P(\lambda|m,z) \propto \frac{1}{\lambda}\exp\left(-\frac{[\ln\lambda - \langle \ln \lambda|m, z\rangle]^2}{\sigma_{\ln\lambda|m,z}^2}\right)
$$

## Bayesian inference pipeline
In this work we aim at infering the parameters of the DC2 redMaPPer mass-richness relation by drawing the posterior distribution

$$
\mathcal{P}(\theta|\mathrm{data}) = \frac{\mathcal{L}(\mathrm{data}|\theta)|\pi(\theta)}{\mathcal{L}(\mathrm{data})}
$$

- `CL_COUNT_class_likelihood`: provides binned Gaussian, Poissonian and unbinned Poissonian likelihoods for cluster count cosmology. We use the count likelihood

$$
\mathcal{L}_{\rm N} \propto \mathrm{det}(\Sigma_{\rm N})^{-1}\exp \left( -\frac{1}{2}[N-\widehat{N}]^T\Sigma^{-1}_{\rm N}[N-\widehat{N}] \right)
$$

and we use a Gaussian likelihood for either the stacked lensing profiles of stacked lensign masses.
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

