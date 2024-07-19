#!/usr/bin/bash
# SLURM options:
#SBATCH --job-name=M    # Job name
#SBATCH --output=log/%x-%j.log
#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=5                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=7000                    # Memory in MB per default
#SBATCH --time=0-6:00:00             # 7 days by default on htc partition
#SBATCH --array=0-8
ID=$SLURM_ARRAY_TASK_ID
source /pbs/home/c/cpayerne/setup_mydesc.sh

python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/cluster_scaling_relation_DC2/modules/run_mcmc_array_job.py M $ID



#likelihood = WL
#0: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#1: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#2: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Prada12_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#3: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Bhattacharya13_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#4: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_+2-halo_rmin=1.0-rmax=15_photoz=Truez_hmf=Despali16
#5: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Bocquet16
#6: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=BPZ_hmf=Despali16
#7: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=flex_hmf=Despali16
#8: MCMC_fit_m-r_WL_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5+cov(g,richness)_photoz=Truez_hmf=Despali16


#likelihood = WLxN
#0: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#1: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#2: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Prada12_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#3: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Bhattacharya13_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#4: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Diemer15_+2-halo_rmin=1.0-rmax=15_photoz=Truez_hmf=Despali16
#5: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Bocquet16
#6: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=BPZ_hmf=Despali16
#7: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=flex_hmf=Despali16
#8: MCMC_fit_m-r_WLxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5+cov(g,richness)_photoz=Truez_hmf=Despali16


#likelihood = M
#0: MCMC_fit_m-r_M_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#1: MCMC_fit_m-r_M_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#2: MCMC_fit_m-r_M_1-halo=nfw+c-M=Prada12_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#3: MCMC_fit_m-r_M_1-halo=nfw+c-M=Bhattacharya13_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#4: MCMC_fit_m-r_M_1-halo=nfw+c-M=Diemer15_+2-halo_rmin=1.0-rmax=15_photoz=Truez_hmf=Despali16
#5: MCMC_fit_m-r_M_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Bocquet16
#6: MCMC_fit_m-r_M_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=BPZ_hmf=Despali16
#7: MCMC_fit_m-r_M_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=flex_hmf=Despali16
#8: MCMC_fit_m-r_M_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5+cov(g,richness)_photoz=Truez_hmf=Despali16


#likelihood = MxN
#0: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Duffy08_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#1: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#2: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Prada12_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#3: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Bhattacharya13_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Despali16
#4: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Diemer15_+2-halo_rmin=1.0-rmax=15_photoz=Truez_hmf=Despali16
#5: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=Truez_hmf=Bocquet16
#6: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=BPZ_hmf=Despali16
#7: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5_photoz=flex_hmf=Despali16
#8: MCMC_fit_m-r_MxN_1-halo=nfw+c-M=Diemer15_rmin=1.0-rmax=5.5+cov(g,richness)_photoz=Truez_hmf=Despali16


#likelihood = N
#0: MCMC_fit_m-r_N__hmf=Despali16
#1: MCMC_fit_m-r_N__hmf=Despali16
#2: MCMC_fit_m-r_N__hmf=Despali16
#3: MCMC_fit_m-r_N__hmf=Despali16
#4: MCMC_fit_m-r_N__hmf=Despali16
#5: MCMC_fit_m-r_N__hmf=Bocquet16
#6: MCMC_fit_m-r_N__hmf=Despali16
#7: MCMC_fit_m-r_N__hmf=Despali16
#8: MCMC_fit_m-r_N__hmf=Despali16


