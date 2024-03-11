#!/usr/bin/bash

# SLURM options:

#SBATCH --job-name=cm    # Job name
#SBATCH --output=log.log
#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=5                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=7000                    # Memory in MB per default
#SBATCH --time=0-03:00:00             # 7 days by default on htc partition
# SBATCH --array=1-4
ID=0
echo startof $ID
# Print the task and run range

source /pbs/home/c/cpayerne/setup_mydesc.sh
python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/richness_mass_from_DC2/mcmc.py


