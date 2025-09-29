## How to run the pipeline ?
# Fit cluster masses
- Add a configuration in `_analysis_cluster_mass_measurement.py`
- Do `python write_run_cluster-mass_measurement.py`, that will fit the stacked masses according the configurations
# Fit the cluster scaling relation
- Add a configuration you want to test in `_analysis_scaling_relation.py`.
- `python write_run_cluster-scaling-relation.py` writes and saves jobs as `job_{config}_{mode}*`, where `config = {WL, M, N, WLxN, MxN}`. In `write_run_cluster-scaling-relation.py`. Here, `mode={argparser, job_array}` is the method of submission. 
- sbatch `job_{config}_{mode}.sh` to submit the job at [CC-IN2P3](https://cc.in2p3.fr/en/).
