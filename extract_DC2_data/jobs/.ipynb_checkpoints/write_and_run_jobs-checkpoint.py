import os, sys
import numpy as np
sys.path.append('..')
name_lens_cat = '/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/data/lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl'
ra_name, dec_name, z_name = 'ra', 'dec', 'redshift'
obs_name = 'richness'
lens_cat = np.load(name_lens_cat, allow_pickle = True)
lens_cat_to_extract = lens_cat[(lens_cat['richness'] > 20)*(lens_cat['redshift'] > .2)]
n_cl = len(lens_cat_to_extract)
index_cl = np.arange(n_cl)
split_lists = np.array_split(index_cl, 50)
#split_lists = [[0, 5]]
for split_list in split_lists:
    start, end = min(split_list), max(split_list)
    lines_base = [
        '#!/bin/sh',
        '# SBATCH options:',
        '#SBATCH --job-name=extract_cosmodc2_redmapper    # Job name',
        f'#SBATCH --output=Logs/log_cluster_index_from_{start}_to_{end}.log',
        '#SBATCH --partition=htc               # Partition choice',
        '#SBATCH --ntasks=1                    # Run a single task (by default tasks == CPU)',
        '#SBATCH --mem=8000                    # Memory in MB per default',
        '#SBATCH --time=0-5:00:00             # 7 days by default on htc partition',
        'source /pbs/home/c/cpayerne/setup_mydesc.sh']
    cmd = [F'python ../run_extract_sources_in_cosmoDC2.py {start} {end}'] 
    lines = lines_base + cmd
    name_job = f'job_source_gal_extract_clusters_from_{start}_to_{end}.job'
    with open(name_job, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')
    os.system(f'sbatch {name_job}')