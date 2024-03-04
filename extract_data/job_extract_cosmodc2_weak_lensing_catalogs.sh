#!/bin/sh

# SBATCH options:

# SBATCH --job-name=extract_cosmodc2_redmapper    # Job name
# SBATCH --output=log.log
# SBATCH --partition=htc               # Partition choice
# SBATCH --ntasks=1                    # Run a single task (by default tasks == CPU)
# SBATCH --mem=8000                    # Memory in MB per default
# SBATCH --time=0-5:00:00             # 7 days by default on htc partition
# SBATCH --array=1-50

N_JOB=$SLURM_ARRAY_TASK_MAX
N_POINTS=4240
N_PER_JOB=$(($N_POINTS / $N_JOB))

START=0
START_NUM=$(($START+($SLURM_ARRAY_TASK_ID-1) * $N_PER_JOB ))
END_NUM=$(($START+($SLURM_ARRAY_TASK_ID) * $N_PER_JOB-1 + 1))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs forthe $START_NUM to $END_NUM

source /pbs/home/c/cpayerne/setup_mydesc.sh
python /pbs/throng/lsst/users/cpayerne/CLMassDC2/data/data_extraction/job/cosmodc2_weak_lensing_catalogs.py $START_NUM $END_NUM 

# python /pbs/throng/lsst/users/cpayerne/CLMassDC2/data/data_extraction/job/cosmodc2_weak_lensing_catalogs.py $DOWN $UP
