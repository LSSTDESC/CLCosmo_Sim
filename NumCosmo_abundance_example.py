try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')

except:
  pass
import math
import matplotlib.pyplot as plt
from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm 
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")

#
#  New homogeneous and isotropic reionization object.
#
reion = Nc.HIReionCamb.new () 

#
#  New homogeneous and isotropic primordial object.
#
prim = Nc.HIPrimPowerLaw.new () 

#
# Adding submodels to the main cosmological model.
#
cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

#
#  New cosmological distance objects optimizied to perform calculations
#  up to redshift 2.0.
#
dist = Nc.Distance.new (2.0)

#
# New transfer function 'NcTransferFuncEH' using the Einsenstein, Hu
# fitting formula.
#
tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")

#
# New linear matter power spectrum object based of the EH transfer function.
# 
psml = Nc.PowspecMLTransfer.new (tf)
psml.require_kmin (1.0e-6)
psml.require_kmax (1.0e3)

#
# Apply a tophat filter to the psml object, set best output interval.
#
psf = Ncm.PowspecFilter.new (psml, Ncm.PowspecFilterType.TOPHAT)
psf.set_best_lnr0 ()

#
# New multiplicity function 'NcMultiplicityFuncTinkerMean'
#
mulf = Nc.MultiplicityFuncTinker.new ()
mulf.set_mdef (Nc.MultiplicityFuncMassDef.MEAN)
mulf.set_Delta (200.0)

#
# New mass function object using the objects defined above.
#
mf = Nc.HaloMassFunction.new (dist, psf, mulf)

#
# New Cluster Mass object using Log normal distribution
#
lnMobs_min = math.log (1.0e14)
lnMobs_max = math.log (1.0e15)
cluster_m = Nc.ClusterMass.new_from_name ("NcClusterMassLnnormal{'lnMobs-min':<%20.15e>, 'lnMobs-max':<%20.15e>}" % (lnMobs_min, lnMobs_max))

#
# New Cluster Redshift object using a global gaussian distribution
#
z_min = 0.3
z_max = 0.7
cluster_z = Nc.ClusterRedshift.new_from_name ("NcClusterPhotozGaussGlobal{'pz-min':<%20.15e>, 'pz-max':<%20.15e>, 'z-bias':<0.0>, 'sigma0':<0.03>}" % (z_min, z_max))

#
# New Cluster abundance object that uses all objects above
#
cad = Nc.ClusterAbundance.new (mf, None)

#
# New NcmData object for number count calculations
#
ncdata = Nc.DataClusterNCount.new (cad)

#
#  Creating a new Modelset and set cosmo as the HICosmo model to be used
#  and cluster_m as the distribution of the mass-observable relation
#
mset = Ncm.MSet.new_array ([cosmo, cluster_z, cluster_m])

#
#  Setting values for the cosmological model, those not set stay in the
#  default values. Remember to use the _orig_ version to set the original
#  parameters when a reparametrization is used.
#
cosmo.props.H0      = 70.0
cosmo.props.Omegab  = 0.05
cosmo.props.Omegac  = 0.25
cosmo.props.Omegax  = 0.70
cosmo.props.Tgamma0 = 2.72
cosmo.props.w       = -1.0

cosmo.param_set_lower_bound (Nc.HICosmoDESParams.OMEGA_C, 0.12)

#
#  Setting values for the mass distribution model
#
cluster_m.props.bias       = 0.0
cluster_m.props.sigma      = 0.2

#
#  Printing the parameters used.
#
mset.pretty_log ()

#
# Creates a new random number generator from a pool named "example_ca_sampling"
# it implicitly creates this pool.
#
rng = Ncm.RNG.pool_get ("example_ca_sampling");

#
# Since ncdata is currently empty, run init_from_sampling
# using the objects above and an survey area of 300deg^2
#
ncdata.init_from_sampling (mset, 300.0 * (math.pi / 180.0)**2, rng)
ncdata.true_data (True)

#
# Save to a fits file
#
ncdata.catalog_save ("ca_data.fits", True)

#
# Generate another sample by resampling from mset
#
ncdata.resample (mset, rng)

cosmo.props.Omegac        = 0.2
cosmo.props.Omegac_fit    = True
cosmo.props.Omegab_fit    = False
#cosmo.props.w_fit         = True
cosmo.props.w_fit         = False
prim.props.ln10e10ASA_fit = True

#cluster_m.props.sigma_fit  = True
#snia = Nc.DataDistMu.new_from_id (dist, Nc.DataSNIAId.SIMPLE_UNION2_1)

dset = Ncm.Dataset.new ()
dset.append_data (ncdata)
#dset.append_data (snia)

#
# New likelihood object using dset
#
lh = Ncm.Likelihood.new (dset)

#
#  Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
#  fit the Modelset mset using the Likelihood lh and using a numerical differentiation
#  algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
#
fit = Ncm.Fit.new (Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD)
fit.run (Ncm.FitRunMsgs.FULL)
#
# Printing fitting informations.
#
fit.log_info ()
#
# Setting single thread calculation.
#
Ncm.func_eval_set_max_threads (5)
Ncm.func_eval_log_pool_stats ()
#
# Additional functions
#
mfunc_oa = Ncm.ObjArray.new ()
mfunc_sigma8 = Ncm.MSetFuncList.new ("NcHICosmo:sigma8", psf)
mfunc_Omegam = Ncm.MSetFuncList.new ("NcHICosmo:Omega_m0", None)
mfunc_oa.add (mfunc_sigma8)
mfunc_oa.add (mfunc_Omegam)
print (mfunc_sigma8.eval0 (mset))
print (mfunc_Omegam.eval0 (mset))

#
# New Gaussian prior to provide the initial points for the chain.
# It was created with size 0 (number of parameters), but once 
# initialized with mset the correct size is assigned. 
#
# The initial sampler will use a diagonal covariance with the
# diagonal terms being the parameters scale set by each model.
#
init_sampler = Ncm.MSetTransKernGauss.new (0)
init_sampler.set_mset (mset)
init_sampler.set_prior_from_mset ()
init_sampler.set_cov_from_rescale (1.0) #1

#
# Creates the ESMCMC walker object, this object is responsible
# for moving the walkers in each interation, the stretch move
# is affine invariant and therefore gives good results even for
# very correlated parametric space.
# 
sampler = 'apes'
#sampler  = 'stretch'
nwalkers = int (math.ceil (800)) #500
ssize    = 1000000 #1000000

if sampler == 'apes':
  walker = Ncm.FitESMCMCWalkerAPES.new (nwalkers, mset.fparams_len ())
elif sampler == "stretch":
  walker = Ncm.FitESMCMCWalkerStretch.new (nwalkers, mset.fparams_len ())
#
# The methods below set the walk scale, which controls the size of the
# step done between two walkers and circumscribe the walkers inside
# the box defined by the parameters inside the mset object.
#
#walker.set_scale (3.0)
#walker.set_box_mset (mset)
#
# Initialize the ESMCMC object using the objects above. It will
# use 50 walkers, i.e., each point in the MCMC chain contains
# 50 points in the parametric space. Each step uses the last point
# in the chain (the last 50 parametric points) to calculate the
# proposal points.
#
esmcmc  = Ncm.FitESMCMC.new_funcs_array (fit, nwalkers, init_sampler, walker, Ncm.FitRunMsgs.SIMPLE, mfunc_oa)

#
# These methods enable the auto-trim options on ESMCMC. This option 
# makes the sampler check the chains' health and trim any unnecessary 
# burn-in part. We set the number of divisions to 100 so we test the
# chains in blocks of n/100. The last method asserts that each 2min
# the catalog will be checked.
#
#esmcmc.set_auto_trim (True)
#esmcmc.set_auto_trim_div (100)
#esmcmc.set_max_runs_time (2.0 * 60.0)
#esmcmc.set_nthreads (4)
esmcmc.set_data_file ("Teste_NC_%s_st_%d_true.fits" % (sampler, nwalkers))

#
# Running the esmcmc, it will first calculate 1000 points, after that
# it will estimate the error in the parameters mean. Using the current
# errors the algorithm tries to calculated how many extra steps are 
# necessary to obtain the required error `10^-3' in every parameters,
# and it will run such extra steps. It will repeat this procedure
# until it attains the required error in every parameter.
# 
#
esmcmc.start_run ()
#esmcmc.run (ssize / nwalkers)
esmcmc.run_lre (30, 1.e-3)
esmcmc.end_run ()
#
# Calculates the parameter means and covariance and set it into 
# the fit object and then print.
# 
esmcmc.mean_covar ()
fit.log_covar ()