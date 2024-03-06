# Modules
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
