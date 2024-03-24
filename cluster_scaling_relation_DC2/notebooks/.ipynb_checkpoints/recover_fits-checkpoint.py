import numpy as np
import glob
def best_fit(config, fit):
    
    names = [config[fit][i]['name'] for i in range(len(config[fit]))]
    print(names)
    where = '../chains/'
    name_plot = []
    p = []
    p_err = []
    for i in range(len(names)):
        name_file = names[i] + '.pkl'
        f = glob.glob(f'/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/cluster_scaling_relation_DC2/chains/{fit}/*')
        w = f'/pbs/throng/lsst/users/cpayerne/CLCosmo_Sim/cluster_scaling_relation_DC2/chains/{fit}/'
        if w+names[i]+'.pkl' not in f: continue
        data = np.load(w+names[i]+'.pkl', allow_pickle=True)['flat_chains']
        param = np.mean(data[15000:], axis=0)
        param_err = np.std(data[15000:], axis=0)
        p.append(param)
        p_err.append(param_err)
        name_plot.append(config[fit][i]['name_plot'])
    return np.array(p), np.array(p_err), name_plot