import clmm

#compute lensing profiles
ra_name, dec_name, z_name, obs_name = 'ra', 'dec', 'redshift', 'richness'
bin_edges = clmm.dataops.make_bins(0.5, 10, 15, method='evenlog10width')
label_pz = ['true']
label_prf = ['DSt', 'DSx', 'W_l', 'radius']
names_cl=['id', ra_name, dec_name, z_name, obs_name]
label_prf_full = [label_prf_ + '_' + label_pz_ for label_pz_ in label_pz for label_prf_ in label_prf]
names = names_cl + label_prf_full
suff_lensing = None