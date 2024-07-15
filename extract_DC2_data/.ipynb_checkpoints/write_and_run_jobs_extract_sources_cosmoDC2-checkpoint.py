import os, sys
import numpy as np
n_splits = 60
catalog = '../data/random_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl'
for index_split in np.arange(n_splits):
    lines_base = [
        '#!/bin/sh',
        '# SBATCH options:',
        '#SBATCH --job-name=extract_cosmodc2_redmapper    # Job name',
        f'#SBATCH --output=Logs/{index_split}_nsplits_{n_splits}.log',
        '#SBATCH --partition=htc               # Partition choice',
        '#SBATCH --ntasks=1                    # Run a single task (by default tasks == CPU)',
        '#SBATCH --mem=8000                    # Memory in MB per default',
        '#SBATCH --time=0-6:00:00             # 7 days by default on htc partition',
        'source /pbs/home/c/cpayerne/setup_mydesc.sh']
    cmd =   f'python run_extract_sources_in_cosmoDC2.py --index_split {index_split} --n_splits {n_splits} '
    cmd +=  f'--lens_catalog_name {catalog} '
    cmd +=  f'--compute_individual_lensing_profile True --save_catalog False '
    lines = lines_base + [cmd]
    name_job = f'job_index_split={index_split}_n_splits={n_splits}.job'
    with open(name_job, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')
    os.system(f'sbatch {name_job}')
    os.remove(name_job)