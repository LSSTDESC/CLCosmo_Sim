import clmm
#source selection

mag_i_max =  24.5
mag_r_max =  28
def source_selection_magnitude(table):
    mask = table['mag_i'] < mag_i_max
    mask *= table['mag_r'] < mag_r_max
    return table[mask]
p_background_min = 0.8

#compute lensing profiles
ra_name, dec_name, z_name = 'ra', 'dec', 'redshift'
obs_name = 'richness'
bin_edges = clmm.dataops.make_bins(0.5, 10, 15, method='evenlog10width')
label_pz = ['true', 'flex', 'bpz']
label_prf = ['DSt', 'DSx', 'W_l', 'radius']
names_cl=['id', ra_name, dec_name, z_name, obs_name]
label_prf_full = [label_prf_ + '_' + label_pz_ for label_pz_ in label_pz for label_prf_ in label_prf]
print(label_prf_full)
names = names_cl + label_prf_full
suff_lensing = None