import argparse
import _mcmc_cluster_scaling_relation_from_lensing_profiles

def collect_argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--type", type=str, required=False, default='N')
    parser.add_argument("--fit_cosmo", type=str, required=False, default='False')
    parser.add_argument("--hmf", type=str, required=False, default='Bocquet16')
    parser.add_argument("--radius_max", type=float, default=1.0)
    parser.add_argument("--radius_min", type=float, default=5.0 )
    parser.add_argument("--cM_relation", type=str, default='Duffy08')
    parser.add_argument("--two_halo",type=str, required=False, default='False')
    parser.add_argument("--photoz", type=str, required=False, default='Truez')
    parser.add_argument("--lensing_data", type=str, required=False, 
                        default='../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl')
    parser.add_argument("--mass_file", type=str, required=False, 
                        default='../../cluster_mass_measurement/cluster-masses_1-halo=nfw+c-M=Duffy08_+2-halo_rmin=1.0-rmax=5.5_photoz=Truez.pkl')
    return parser.parse_args()

_config_extract_sources_in_cosmoDC2 = collect_argparser()
analysis_metadata = vars(collect_argparser())
_mcmc_cluster_scaling_relation_from_lensing_profiles.mcmc(analysis_metadata)
