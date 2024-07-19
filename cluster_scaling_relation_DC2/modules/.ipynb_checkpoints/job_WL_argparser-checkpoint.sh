#!/usr/bin/bash
# SLURM options:
#SBATCH --job-name=WL    # Job name
#SBATCH --output=log/%x-%j.log
#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=5                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=7000                    # Memory in MB per default
#SBATCH --time=0-6:00:00             # 7 days by default on htc partition
#SBATCH --array=0-8
ID=$SLURM_ARRAY_TASK_ID
source /pbs/home/c/cpayerne/setup_mydesc.sh
#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Duffy08 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Duffy08 --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16 

#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Diemer15 --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16 

#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Prada12_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Prada12 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Prada12 --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Prada12_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16 

#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Bhattacharya13_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Bhattacharya13 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Bhattacharya13 --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Bhattacharya13_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16 

#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_+2-halo_rmin=1.0-rmax=15_photoz=Truez_hmf=Despali16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo True --hmf Despali16 --radius_max 15 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Diemer15 (2h) --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_+2-halo_rmin=1.0-rmax=15_photoz=Truez_hmf=Despali16 

#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Bocquet16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Bocquet16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw - Bocquet16 --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Bocquet16 

#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=BPZ_hmf=Despali16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz BPZ --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_BPZ_full_coverage.pkl --name_plot BPZ --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=BPZ_hmf=Despali16 

#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=flex_hmf=Despali16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz flex --shear_richness_cov False --lensing_data ../../data/stacked_esd_profiles_redmapper_flex_full_coverage.pkl --name_plot FleXZBoost --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=flex_hmf=Despali16 

#MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5+cov(g,richness)_photoz=Truez_hmf=Despali16
#python run_mcmc_argparser.py--type WL --fit_cosmo False --density_profile nfw --cM_relation Diemer15 --two_halo False --hmf Despali16 --radius_max 5.5 --radius_min 1.0 --photoz Truez --shear_richness_cov True --lensing_data ../../data/stacked_esd_profiles_redmapper_true_full_coverage.pkl --name_plot nfw + Cov($\Delta\Sigma,\lambda$) --name MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5+cov(g,richness)_photoz=Truez_hmf=Despali16 

