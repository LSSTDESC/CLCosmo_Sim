#!/usr/bin/bash
# SLURM options:
#SBATCH --job-name=M    # Job name
#SBATCH --output=log/%x-%j.log
#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=5                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=7000                    # Memory in MB per default
#SBATCH --time=0-10:00:00             # 7 days by default on htc partition
#SBATCH --array=0-4
ID=$SLURM_ARRAY_TASK_ID
source /pbs/home/c/cpayerne/setup_mydesc.sh

python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/cluster_scaling_relation_DC2/modules/run_mcmc_array_job.py M $ID



#likelihood = WL
#0: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=Truez_hmf=Despali16
#1: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=BPZ_hmf=Despali16
#2: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=flex_hmf=Despali16
#3: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Duffy08+lensing_bias_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=BPZ_hmf=Despali16
#4: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Duffy08+lensing_bias_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=flex_hmf=Despali16


#likelihood = WLxN
#0: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=Truez_hmf=Despali16
#1: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=BPZ_hmf=Despali16
#2: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=flex_hmf=Despali16
#3: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Duffy08+lensing_bias_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=BPZ_hmf=Despali16
#4: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Duffy08+lensing_bias_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=flex_hmf=Despali16


#likelihood = M
#0: MCMC_fit_m-r_M_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=Truez_hmf=Despali16
#1: MCMC_fit_m-r_M_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=BPZ_hmf=Despali16
#2: MCMC_fit_m-r_M_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=flex_hmf=Despali16
#3: MCMC_fit_m-r_M_1-halo=nfw+c-M=Duffy08+lensing_bias_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=BPZ_hmf=Despali16
#4: MCMC_fit_m-r_M_1-halo=nfw+c-M=Duffy08+lensing_bias_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=flex_hmf=Despali16


#likelihood = MxN
#0: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=Truez_hmf=Despali16
#1: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=BPZ_hmf=Despali16
#2: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=flex_hmf=Despali16
#3: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Duffy08+lensing_bias_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=BPZ_hmf=Despali16
#4: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Duffy08+lensing_bias_rmin=1.0-rmax=3.5+cov(g,richness)_photoz=flex_hmf=Despali16


#likelihood = N
#0: MCMC_fit_m-r_N_hmf=Despali16
#1: MCMC_fit_m-r_N_hmf=Despali16
#2: MCMC_fit_m-r_N_hmf=Despali16
#3: MCMC_fit_m-r_N_hmf=Despali16
#4: MCMC_fit_m-r_N_hmf=Despali16


