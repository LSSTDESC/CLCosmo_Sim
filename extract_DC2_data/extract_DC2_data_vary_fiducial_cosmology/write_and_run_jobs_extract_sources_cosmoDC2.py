import os, sys
import numpy as np
n_splits = 60
catalog = '../../../CLCosmo_Sim_database/data/lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl'
for index_split in np.arange(n_splits):
    lines_base = [
        '#!/bin/sh',
        '# SBATCH options:',
        '#SBATCH --job-name=extract_cosmodc2_redmapper    # Job name',
        f'#SBATCH --output=Logs/{index_split}_nsplits_{n_splits}.log',
        '#SBATCH --partition=htc               # Partition choice',
        '#SBATCH --ntasks=1                    # Run a single task (by default tasks == CPU)',
        '#SBATCH --mem=8000                    # Memory in MB per default',
        '#SBATCH --time=0-5:00:00             # 7 days by default on htc partition',
        'source /pbs/home/c/cpayerne/setup_mydesc.sh']
#python run_extract_sources_in_cosmoDC2_vary_cosmology_gammat.py --index_split 0 --n_splits 3 --lens_catalog_name ../../../CLCosmo_Sim_database/data/lens_catalog_cosmoDC2_v1.1.4_redmapper_v0.8.1.pkl --compute_individual_lensing_profile True --save_catalog False
    cmd =   f'python run_extract_sources_in_cosmoDC2_vary_cosmology_gammat.py --index_split {index_split} --n_splits {n_splits} '
    cmd +=  f'--lens_catalog_name {catalog} '
    cmd +=  f'--compute_individual_lensing_profile True --save_catalog False '
    lines = lines_base + [cmd]
    name_job = f'job_index_split={index_split}_n_splits={n_splits}.job'

    #file1 = f'/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim_database/data_vary_fuducial_cosmology/ind_gammat_profile_redmapper_per_cluster_index/ind_profile_redmapper_split={index_split}_nsplits={n_splits}_ncl=71.pkl'
   # file2 = f'/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim_database/data_vary_fuducial_cosmology/ind_profile_redmapper_per_cluster_index/ind_profile_redmapper_split={index_split}_nsplits={n_splits}_ncl=70.pkl'

   # if os.path.isfile(file1):
   #     continue
   # if os.path.isfile(file2):
   #     continue
    with open(name_job, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')
    #break
    os.system(f'sbatch {name_job}')
    os.remove(name_job)
    #break
