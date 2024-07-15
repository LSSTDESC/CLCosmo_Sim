#source selection

mag_i_max =  24.25
mag_r_max =  28
def source_selection_magnitude(table):
    mask = table['mag_i'] < mag_i_max
    mask *= table['mag_r'] < mag_r_max
    return table[mask]
p_background_min = 0.9
Dz_background = 0.2