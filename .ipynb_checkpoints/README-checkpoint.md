# CLCosmo_Sim repository

Cluster cosmology analysis using DESC tools and simulations

# Modules
## Modeling cluster abundance and cluster lensing
- `CL_COUNT_cluster_abundance`: provides prediction for cluster abundance in the proxy-redshift space (binned) with observable-mass relation, completeness and purity parametrisation.
- `CL_LENSING_cluster_lensing`: provides prediction for cluster lensing signal in the proxy-redshift space (binned) with observable-mass relation, completeness and purity parametrisation.
- `CL_COUNT_modeling_halo_mass_functio`: modeling the halo mass function and cosmoving volume (CCL as backend).
- `CL_COUNT_modeling_purity`: modeling purity.
- `CL_COUNT_modeling_completenes`: modeling completeness.
- `CL_COUNT_modeling_richness_mass_relation`: modeling the observable-mass relation (mean, disperion, probability density function).

## Bayesian inference pipeline
- `CL_COUNT_class_likelihood`: provides binned Gaussian, Poissonian and unbinned Poissonian likelihoods for cluster count cosmology.
- `CL_COUNT_DATAOPS_cluster_abundance_covariance`: provides estimator of cluster abundance covariance based on data.
- `STAT_forecast`: module for Fisher forecast.
