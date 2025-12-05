#cosmology

H0_true = 71
h = H0_true/100
Omega_b_true = 0.02258 / (h**2)
Omega_c_true = 0.1109 / (h**2)
Omega_m_true = Omega_b_true + Omega_c_true
sigma8_true = 0.8
ns_true = 0.963
cosmo_clmm = Cosmology(H0 = H0_true, Omega_dm0 = Omega_c_true, Omega_b0 = Omega_b_true, Omega_k0 = 0.0)
cosmo = ccl.Cosmology(Omega_c = Omega_c_true, Omega_b = Omega_b_true, h = H0_true/100, sigma8 = sigma8_true, n_s=ns_true)

halo_profile = 'nfw', 'einasto'
massdef = ccl.halos.massdef.MassDef(200, 'critical',)
concDiemer15 = ccl.halos.concentration.ConcentrationDiemer15(mass_def=massdef)
concDuffy08 = ccl.halos.concentration.ConcentrationDuffy08(mass_def=massdef)
concPrada12 = ccl.halos.concentration.ConcentrationPrada12(mass_def=massdef)
concBhattacharya13 = ccl.halos.concentration.ConcentrationBhattacharya13(mass_def=massdef)

Omega_m_list = np.linspace(0.1, 0.6, 30)

radial_range = [1, 3.5] #in Mpc at fiducial cosmology

file_format ? 
Astropy_table with colnames ['z_mean', 'obs_mean', 'log10M200c_WL', 'err_log10M200c_WL'] (nz * nrich size table) for each Omega_m value