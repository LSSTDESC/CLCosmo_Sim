import os, sys
import _analysis_scaling_relation

import os, sys
import _analysis_scaling_relation
line= []
for config in _analysis_scaling_relation.config.keys():
    line+=['likelihood = '+ config]
    for i in range(len(_analysis_scaling_relation.config[config])):
        line+=[str(i) + ': ' + _analysis_scaling_relation.config[config][i]['name']]
    line+=['================================================', '']

with open('inventory.txt', 'w') as f:
    for line_ in line:
        f.write(line_)
        f.write('\n')


# for config in _analysis_scaling_relation.config.keys():
#     n = len(_analysis_scaling_relation.config[config])
#     lines_base = [
#         '#!/usr/bin/bash',
#         '# SLURM options:',
#        f'#SBATCH --job-name={config}    # Job name',
#         '#SBATCH --output=log/%x-%j.log',
#         '#SBATCH --partition=htc               # Partition choice',
#         F'#SBATCH --ntasks={n}                    # Run a single task (by default tasks == CPU)',
#         '#SBATCH --mem=7000                    # Memory in MB per default',
#         '#SBATCH --time=0-10:00:00             # 7 days by default on htc partition',
#         'source /pbs/home/c/cpayerne/setup_mydesc.sh','','','']

#     cmd = []
#     cmd += lines_base 
#     python_cmd = f'#srun -n 1 python run_mcmc_argparser.py'
#     for case in range(n):
#         keys_ = _analysis_scaling_relation.config[config][case].keys()
#         cmd_analysis_k = '' + python_cmd
#         for k in keys_:
#             cmd_analysis_k += f'--{k} {_analysis_scaling_relation.config[config][case][k]} '
#             analysis_name = _analysis_scaling_relation.config[config][case]['name']
#         cmd += ['#'+analysis_name, cmd_analysis_k, '']
#     name_job = f'job_{config}_argparser.sh'
#     #with open(name_job, 'w') as f:
#     #    for line in cmd:
#      #       f.write(line)
#      #       f.write('\n')
#     #print('[write job file]: '+config+' --> saved (argparser) !')
#     #os.system(f'sbatch {name_job}')
#     #os.remove(name_job)


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
        '#SBATCH --time=0-10:00:00             # 7 days by default on htc partition',
        f'#SBATCH --array=0-{n-1:.0f}',
        'ID=$SLURM_ARRAY_TASK_ID',
        'source /pbs/home/c/cpayerne/setup_mydesc.sh']
    cmd = ['',f'python /pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/cluster_scaling_relation_DC2/modules/run_mcmc_array_job.py {config} $ID']
    lines = lines_base + cmd

    print('[write job file]: '+config+' --> saved (array job) !')
    line= ['','','']
    for config_ in _analysis_scaling_relation.config.keys():
        line+=['#likelihood = '+ config_]
        for i in range(len(_analysis_scaling_relation.config[config_])):
            line+=['#'+str(i) + ': ' + _analysis_scaling_relation.config[config_][i]['name']]
        line+=['', '']
    lines+=line
    name_job = f'job_{config}_array_job.sh'
    with open(name_job, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')

    #os.system(f'sbatch {name_job}')
    #os.remove(name_job)
