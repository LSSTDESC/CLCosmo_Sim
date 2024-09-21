import numpy as np
from itertools import combinations, chain
import healpy
import PySSC

def binning(corner): return [[corner[i],corner[i+1]] for i in range(len(corner)-1)]

class Covariance_matrix():
    r"""
    Class for the computation of covariance matrices for cluster abundance:
    a. Bootstrap covariance matrix
    b. Jackknife covariance matrix
    c. Sample covariance matrix
    """
    def __init__(self):
        self.name = 'name'
        return None
    
    def compute_theoretical_Sij(self, Z_bin, cosmo, f_sky):
        default_cosmo_params = {'omega_b':cosmo['Omega_b']*cosmo['h']**2, 
                                'omega_cdm':cosmo['Omega_c']*cosmo['h']**2, 
                                'H0':cosmo['h']*100, 
                                'n_s':cosmo['n_s'], 
                                'sigma8': cosmo['sigma8'],
                                'output' : 'mPk'}
        z_arr = np.linspace(0.2,1.,1000)
        nbins_T   = len(Z_bin)
        windows_T = np.zeros((nbins_T,len(z_arr)))
        for i, z_bin in enumerate(Z_bin):
            Dz = z_bin[1]-z_bin[0]
            z_arr_cut = z_arr[(z_arr > z_bin[0])*(z_arr < z_bin[1])]
            for k, z in enumerate(z_arr):
                if ((z>z_bin[0]) and (z<=z_bin[1])):
                    windows_T[i,k] = 1  
        Sij_fullsky = PySSC.Sij_alt_fullsky(z_arr, windows_T, order=1, cosmo_params=default_cosmo_params, cosmo_Class=None, convention=0)
        f_sky = (439.78987/(360**2/np.pi))
        Sij_partialsky = Sij_fullsky/f_sky
        return Sij_partialsky 
    
    def sample_covariance_full_sky(self, Z_bin, Proxy_bin, NBinned_halo_bias, Sij):
        r"""
        returns the sample covariance matrix for cluster count
        Attributes:
        -----------
         Redshift_bin : list of lists
            list of redshift bins
        Proxy_bin : list of lists
            list of mass bins
        Binned_Abundance: array
            predicted abundance
        Binned_halo_bias: array
            predicted binned halo bias
            Sij: array
        matter fluctuation amplitude per redshift bin
        Returns:
        --------
        sample_covariance: array
            sample covariance for cluster abundance
            #uses the calculation of the fluctuation apmplitude Sij
        """
        index_LogM, index_Z =  np.meshgrid(np.arange(len(Proxy_bin)), np.arange(len(Z_bin)))
        index_Z_flatten = index_Z.flatten()
        len_mat = len(Z_bin) * len(Proxy_bin)
        cov_SSC = np.zeros([len_mat, len_mat])
        Nb = NBinned_halo_bias#np.multiply(Binned_Abundance, Binned_halo_bias)
        Nbij = np.tensordot(Nb, Nb)
        for i, Nbi in enumerate(Nb.flatten()):
           # if i%100==0: print(i)
            for j, Nbj in enumerate(Nb.flatten()):
                if i >= j:
                    index_z_i, index_z_j = index_Z_flatten.flatten()[i], index_Z_flatten.flatten()[j]
                    cov_SSC[i,j] = Nbi * Nbj * Sij[index_z_i, index_z_j]
                    cov_SSC[j,i] = cov_SSC[i,j]
        return cov_SSC

    def compute_boostrap_covariance(self, catalog = None, 
                                    proxy_colname = 'mass', redshift_colname = 'redshift', 
                                    proxy_corner = None, z_corner = None, 
                                    n_boot = 100):
        r"""
        Attributes:
        -----------
        catalog: Table
            catalog of clusters (ra, dec, z, proxy, etc...)
        proxy_colname: str
            name of the proxy column
        redshift_colname: str
            name of the redshift column
        proxy_corner: str
            values of proxues to be binned
        z_corner: str
            values of redshifts to be binned
        n_boot: int
            number of bootstrap resampling
        fct_modify: fct
            function to add optional modifications
        Returns:
        --------
        cov_N: array
            bootstrap covariance matrix
        """
        proxy, redshift = np.array(catalog[proxy_colname]), np.array(catalog[redshift_colname])
        index = np.arange(len(proxy))
        data_boot = []
        for i in range(n_boot):
            index_bootstrap = np.random.choice(index, len(index))
            data, proxy_edges, z_edges = np.histogram2d(redshift[index_bootstrap], 
                                                        proxy[index_bootstrap],
                                                   bins=[z_corner, proxy_corner])
            data_boot.append(data.flatten())
        data_boot = np.array(data_boot)
        N = np.stack((data_boot.astype(float)), axis = 1)
        mean = np.mean(data_boot, axis = 0)
        cov_N = np.cov(N, bias = False)
        self.Bootstrap_covariance_matrix = cov_N
        
        return cov_N

    def compute_jackknife_covariance_healpy(self, catalog = None, 
                                            proxy_colname = 'mass', redshift_colname = 'redshift',
                                            z_corner = None, proxy_corner = None, 
                                            ra_colname = 'ra', dec_colname = 'dec', 
                                            n_power = 8, N_delete = 1):
        r"""
        Attributes:
        -----------
        proxy_colname: str
            name of the proxy column
        redshift_colname: str
            name of the redshift column
        proxy_corner: str
            values of proxues to be binned
        ra_colname: str
            name of the ra column
        dec_colname: str
            name of the dec column
        z_corner: str
            values of redshifts to be binned
        n_power: int
            defines the number of healpix pixels
        N_delete: int
            number of jackknife region to delete each repetition
        Returns:
        --------
        cov_N: array
            Jackknife covariance matrix
        """
        proxy, redshift = np.array(catalog[proxy_colname]), np.array(catalog[redshift_colname])
        ra, dec =  np.array(catalog[ra_colname]), np.array(catalog[dec_colname])
        index = np.arange(len(proxy))
        healpix = healpy.ang2pix(2**n_power, ra, dec, nest=True, lonlat=True)
        healpix_list_unique = np.unique(healpix)
        print(f'Number of JK regions: {len(healpix_list_unique)}')
        healpix_combination_delete = list(combinations(healpix_list_unique, N_delete))
        data_jack = []
        for i, hp_list_delete in enumerate(healpix_combination_delete):
            if (i/2000.).is_integer():
                print(i)
            mask_in_area = np.isin(healpix, hp_list_delete)
            mask_out_area = np.invert(mask_in_area)
            #data, mass_edges, z_edges = np.histogram2d(redshift[mask_out_area], 
            #                                           proxy[mask_out_area],
            #                                           bins=[z_corner, proxy_corner]) 
            data, mass_edges, z_edges = np.histogram2d(
                                                       proxy[mask_out_area],redshift[mask_out_area],
                                                       bins=[ proxy_corner, z_corner])
            data_jack.append(data.flatten())
        data_jack = np.array(data_jack)
        N = np.stack((data_jack.astype(float)), axis = 1)
        n_jack = len(healpix_combination_delete)
        cov_N = (n_jack - 1) * np.cov(N, bias = False,ddof=0)
        coeff = (n_jack-N_delete)/(N_delete*n_jack)
        Jackknife_covariance_matrix = cov_N * coeff
        
        return Jackknife_covariance_matrix