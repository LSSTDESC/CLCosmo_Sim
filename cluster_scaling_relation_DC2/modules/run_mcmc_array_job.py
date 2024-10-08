import sys
import _mcmc_cluster_scaling_relation_from_lensing_profiles
import _analysis_scaling_relation as analysis_list

code, config_name, index_analysis = sys.argv
analysis_metadata = analysis_list.config[config_name][int(index_analysis)]
_mcmc_cluster_scaling_relation_from_lensing_profiles.mcmc(analysis_metadata)
