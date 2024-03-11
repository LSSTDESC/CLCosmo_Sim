#!/usr/bin/bash

# SLURM options:

#SBATCH --job-name=cm    # Job name
#SBATCH --output=log/%x-%j.log
#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=5                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=7000                    # Memory in MB per default
#SBATCH --time=0-5:00:00             # 7 days by default on htc partition
#SBATCH --array=0-9
ID=$SLURM_ARRAY_TASK_ID

source /pbs/home/c/cpayerne/setup_mydesc.sh
#python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/richness_mass_from_DC2/posterior_modules/mcmc.py c_M $ID
#python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/richness_mass_from_DC2/posterior_modules/mcmc.py GammaLambda $ID
python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/richness_mass_from_DC2/posterior_modules/mcmc.py $ID

