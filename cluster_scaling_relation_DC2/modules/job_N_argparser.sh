#!/usr/bin/bash
# SLURM options:
#SBATCH --job-name=N    # Job name
#SBATCH --output=log/%x-%j.log
#SBATCH --partition=htc               # Partition choice
#SBATCH --n=9                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=7000                    # Memory in MB per default
#SBATCH --time=0-10:00:00             # 7 days by default on htc partition
source /pbs/home/c/cpayerne/setup_mydesc.sh



#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Duffy08 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Duffy08 --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Diemer15 --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Prada12 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Prada12 --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Bhattacharya13 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Bhattacharya13 --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo True --hmf Despali16 --radius_max 15 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Diemer15 (2h) --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Bocquet16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Bocquet16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Bocquet16 --name MCMC_fit_m-r_N__hmf=Bocquet16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz BPZ --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_BPZ_full_coverage.pkl --name_plot BPZ --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz flex --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_flex_full_coverage.pkl --name_plot FleXZBoost --name MCMC_fit_m-r_N__hmf=Despali16 

#MCMC_fit_m-r_N__hmf=Despali16
#srun -n 1 python run_mcmc_argparser.py--type N --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov True --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw + Cov($\Delta\Sigma,\lambda$) --name MCMC_fit_m-r_N__hmf=Despali16 

