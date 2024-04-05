#!/usr/bin/bash
# SLURM options:
#SBATCH --job-name=WLxN    # Job name
#SBATCH --output=log/%x-%j.log
#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=5                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=7000                    # Memory in MB per default
#SBATCH --time=0-6:00:00             # 7 days by default on htc partition
#SBATCH --array=0-8
ID=$SLURM_ARRAY_TASK_ID
source /pbs/home/c/cpayerne/setup_mydesc.sh

python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/cluster_scaling_relation_DC2/modules/run_mcmc_cluster_scaling_relation_from_lensing_profiles.py WLxN $ID
