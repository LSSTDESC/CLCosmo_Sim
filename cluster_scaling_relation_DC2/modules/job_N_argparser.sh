#!/usr/bin/bash
# SLURM options:
#SBATCH --job-name=N    # Job name
#SBATCH --output=log/%x-%j.log
#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=3                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=7000                    # Memory in MB per default
#SBATCH --time=0-10:00:00             # 7 days by default on htc partition
source /pbs/home/c/cpayerne/setup_mydesc.sh



#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Duffy08 --two_halo False --hmf Despali16 --radius_max 3.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --redshift_range Full --richness_range Full --redshift_bin_index [0 1 2 3 4 5 6] --richness_bin_index [0 1 2 3] --redshift_corner_index [0 1 2 3 4 5 6 7] --richness_corner_index [0 1 2 3 4] --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --add_bias_lensing False --name_plot baseline --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Duffy08 --two_halo False --hmf Despali16 --radius_max 3.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --redshift_range Partial --richness_range Full --redshift_bin_index [0 1 2 3 4] --richness_bin_index [0 1 2 3] --redshift_corner_index [0 1 2 3 4 5] --richness_corner_index [0 1 2 3 4] --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --add_bias_lensing False --name_plot baseline - low z --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Duffy08 --two_halo False --hmf Despali16 --radius_max 3.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --redshift_range Full --richness_range Partial --redshift_bin_index [0 1 2 3 4 5 6] --richness_bin_index [0 1] --redshift_corner_index [0 1 2 3 4 5 6 7] --richness_corner_index [0 1 2] --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --add_bias_lensing False --name_plot baseline - low richness --name MCMC_fit_m-r_N__hmf=Despali16 

