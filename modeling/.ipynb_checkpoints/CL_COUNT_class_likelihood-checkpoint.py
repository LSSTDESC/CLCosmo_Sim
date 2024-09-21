import numpy as np
import scipy
from scipy.special import erfc
from scipy.stats import poisson
from scipy.stats import multivariate_normal
from scipy.integrate import quad,simps, dblquad

def binning(corner): return [[corner[i],corner[i+1]] for i in range(len(corner)-1)]

class Likelihood():
    r"""
        compute likelihood :
            a. for the binned gaussian case
            b. for the binned poissonian case
            c. for the un-binned poissonian case
    """
    def ___init___(self):
        self.name = 'Likelihood for cluster count Cosmology'
        
    def lnLikelihood_Binned_Poissonian(self, N_th_matrix, N_obs_matrix):
        r"""
        returns the value of the log-likelihood for Poissonian binned approach
        Attributes:
        -----------
        N_th_matrix: array
            cosmological prediction for binned cluster abundance
        N_obs_matrix:
            observed binned cluster abundance
        Returns:
        --------
        add attributes with total log-likelihood for Poissonian binned approach
        """
        lnL_Poissonian = np.sum(N_obs_matrix.flatten() * np.log(N_th_matrix.flatten()) - N_th_matrix.flatten())
        self.lnL_Binned_Poissonian = lnL_Poissonian

    def lnLikelihood_Binned_Gaussian(self, N_th_matrix, N_obs_matrix, covariance_matrix):
        r"""
        Attributes:
        -----------
        N_th_matrix: array
            cosmological prediction for binned cluster abundance
        N_obs_matrix:
            observed binned cluster abundance
        covariance_matrix: array
            full covariance matrix for binned cluster abundance
        Returns:
        --------
        add attributes with total log-likelihood for Gaussian binned approach
        """
        delta = (N_obs_matrix - N_th_matrix).flatten()
        inv_covariance_matrix = np.linalg.inv((covariance_matrix))
        self.lnL_Binned_Gaussian = -0.5*np.sum(delta*inv_covariance_matrix.dot(delta)) - np.log(np.linalg.det(covariance_matrix)**.5)

        
    def lnLikelihood_UnBinned_Poissonian(self, dN_dzdlogMdOmega, N_tot):
        r"""
        Attributes:
        -----------
       dN_dzdlogMdOmega: array
            cosmological prediction for multiplicu-ity function
        N_tot: float
            cosmological prediction for total number of cluster
        Returns:
        --------
        add attributes with total log-likelihood for Poissonian unbinned approach
        """
        self.lnL_UnBinned_Poissonian = np.sum(np.log(dN_dzdlogMdOmega)) - N_tot
