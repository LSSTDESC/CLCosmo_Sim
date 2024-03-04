# Modules
## cosmoDC2: Extraction of galaxy catalogs behind redMaPPer clusters
- `extract_in_cosmodc2_utils.py` : functions to extract photometric probability distributions using GCRcatalogs in the cosmoDC2 photoz-addons, and "truth" quantities using Qserv queries
- `photoz_utils.py` : From photometric probability density functions, compute mean redshift, mean critical surface mass density for each galaxy
- `cosmodc2_weak_lensing_catalog.py` : module to extract background galaxy catalogs for each redMaPPer cluster, with photozs, each one saved in a pickle file
## Estimate cluster lensing profiles from background source galaxies
- `compute_individual_profile_cosmodc2_redmapper_clusters.py` : compute individual lensing profiles 
