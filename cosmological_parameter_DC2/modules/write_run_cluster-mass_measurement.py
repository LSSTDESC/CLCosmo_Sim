import os, sys
import _analysis_cluster_mass_measurement

for analysis_name in  _analysis_cluster_mass_measurement.analysis.keys():
    analysis = _analysis_cluster_mass_measurement.analysis[analysis_name]
    for i in range(len(analysis)):
        os.system(f'python run_mcmc_cluster-mass_measurement.py {analysis_name} {i}')
