import os, sys
import _analysis_scaling_relation

for config in _analysis_scaling_relation.config.keys():
    n = len(_analysis_scaling_relation.config[config])
    lines_base = [
        '#!/usr/bin/bash',
        '# SLURM options:',
       f'#SBATCH --job-name={config}    # Job name',
        '#SBATCH --output=log/%x-%j.log',
        '#SBATCH --partition=htc               # Partition choice',
        '#SBATCH --ntasks=5                    # Run a single task (by default tasks == CPU)',
        '#SBATCH --mem=7000                    # Memory in MB per default',
        '#SBATCH --time=0-6:00:00             # 7 days by default on htc partition',
        f'#SBATCH --array=0-{n-1:.0f}',
        'ID=$SLURM_ARRAY_TASK_ID',
        'source /pbs/home/c/cpayerne/setup_mydesc.sh']
    cmd = ['',f'python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/cluster_scaling_relation_DC2/modules/run_mcmc_cluster_scaling_relation_from_lensing_profiles.py {config} $ID']
    lines = lines_base + cmd
    name_job = f'job_{config}.sh'
    with open(name_job, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')
    os.system(f'sbatch {name_job}')
    #os.remove(name_job)
