"""
Module mkde provides methods to compute multivariate kernel density
estimations.
(see https://yketa.github.io/UBC_2018_Wiki/#Multivariate%20kernel%20density%20estimation)
"""

import numpy as np

import random

import scipy.optimize as scop

import ctypes

import os

import subprocess

from KDEpy import FFTKDE

# C EXTENSION

_dir_path = os.path.dirname(os.path.realpath(__file__)) # script directory path

_c_ext_so_path = os.path.join(_dir_path, 'cmkde.so')    # path to C extension shared object .so file
_c_ext_c_path = os.path.join(_dir_path, 'cmkde.c')      # path to C extension source .c file
_c_ext_h_path = os.path.join(_dir_path, 'cmkde.h')      # path to C extension source .h file

if (not(os.path.isfile(_c_ext_so_path))                                     # C extention shared object does not exist
    or os.path.getctime(_c_ext_c_path) < os.path.getctime(_c_ext_so_path)   # C extension source .c file is more recent than shared object
    or os.path.getctime(_c_ext_h_path) < os.path.getctime(_c_ext_so_path)): # C extension source .h file is more recent than shared object
    subprocess.run(['gcc', '-shared', '-o', _c_ext_so_path,
        '-std=c99', '-fPIC', _c_ext_c_path],
        cwd=_dir_path)                                                      # compile C extension

class _MKDECExt:
    """
    Wrapper of multivariate kerndel density estimation (MKDE) C extension.
    """

    def __init__(self, data):
        """
        Initialise structure containing all data relevant to computation in the
        C extension.

        Parameters
        ----------
        data : 2D float Numpy array
            Data points.
        """

        # C EXTENSION

        self.cmkde = ctypes.CDLL(_c_ext_so_path)    # C extension shared library

        # DATA FOR COMMUNICATION WITH C EXTENSION

        self.data = np.array(data, dtype=float)
        self.n, self.d = self.data.shape    # number of data points and dimensions

        self.sd = _SampleData(                                              # data structure for communication between C and Python
            ctypes.c_int(self.n),                                           # int n
            ctypes.c_int(self.d),                                           # int d
            (ctypes.POINTER(ctypes.c_double)*self.n)(
                *[np.ctypeslib.as_ctypes(_) for _ in self.data]),           # double **data
            ctypes.POINTER(ctypes.c_double)(ctypes.c_double()),             # double *AMISE
            ctypes.POINTER(ctypes.c_double)((ctypes.c_double*self.d)()),    # double *gradAMISE
            (ctypes.POINTER(ctypes.c_double)*self.d)(*[
                ctypes.POINTER(ctypes.c_double)((ctypes.c_double*self.d)())
                for _ in range(self.d)]),                                   # double **hessAMISE
            ctypes.POINTER(ctypes.c_double)((ctypes.c_double*self.d)()))    # double *h

    def _update_h(self, h):
        """
        Update bandwidths in data structure.

        Parameters
        ----------
        h : 1D array-like
            Array of bandwidths.
        """

        for i in range(self.d):
            self.sd.h[i] = h[i]

    def AMISE(self, h):
        """
        Compute asymptotic mean integrated squared error (AMISE).

        Parameters
        ----------
        h : 1D array-like
            Array of bandwidths.

        Returns
        -------
        AMISE : float
            Asymptotic mean integrated squared error (AMISE).
        """

        self._update_h(h)
        self.cmkde.AMISE(ctypes.byref(self.sd))

        return self.sd.AMISE[0]

    def gradAMISE(self, h):
        """
        Compute gradient of AMISE with respect to the bandwidths.

        Parameters
        ----------
        h : 1D array-like
            Array of bandwidths.

        Returns
        -------
        gradAMISE : 1D Numpy array
            Gradient of AMISE with respect to the bandwidths.
        """

        self._update_h(h)
        self.cmkde.gradAMISE(ctypes.byref(self.sd))

        return np.array([self.sd.gradAMISE[i]
            for i in range(self.d)])

    def hessAMISE(self, h):
        """
        Compute Hessian matrix of AMISE with respect to the bandwidths.

        Parameters
        ----------
        h : 1D array-like
            Array of bandwidths.

        Returns
        -------
        gradAMISE : 2D Numpy array
            Hessian matrix of AMISE with respect to the bandwidths.
        """

        self._update_h(h)
        self.cmkde.hessAMISE(ctypes.byref(self.sd))

        return np.array([[self.sd.hessAMISE[i][j]
            for j in range(self.d)]
            for i in range(self.d)])

class _SampleData(ctypes.Structure):
	"""
	Structure containing all data relevant to computation in the C extension.
	"""

	_fields_ = [
        # data sample
        ('n', ctypes.c_int),                                            # number of data points
        ('d', ctypes.c_int),                                            # number of dimensions
        ('data', ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),      # array of data points
        ('AMISE', ctypes.POINTER(ctypes.c_double)),                     # AMISE
        ('gradAMISE', ctypes.POINTER(ctypes.c_double)),                 # gradient of AMISE
        ('hessAMISE', ctypes.POINTER(ctypes.POINTER(ctypes.c_double))), # Hessian matrix of AMISE
        # optimised bandwidths
        ('h', ctypes.POINTER(ctypes.c_double))]                         # bandwidths

# FUNCTIONS AND CLASSES

class MKDE:
    """
    Perform multivariate kernel density estimation with Gaussian kernel
    functions and diagonal bandwidth matrix. (see
    https://yketa.github.io/UBC_2018_Wiki/#Multivariate%20kernel%20density%20estimation)
    """

    def __init__(self, data):
        """
        Sample data input.

        Parameters
        ----------
        data : array-like
            Data sample which we want to estimate the probability distribution
            function.
        """

        self.data = np.array(data)              # data sample
        if len(self.data.shape) == 1:           # unidimensional sample
            self.data = np.reshape(self.data, (len(self.data), 1))
        self.points, self.d = self.data.shape   # size and dimension of sample data

    def bw(self, n_max=1000, lower_bound=1e-10):
        """
        Compute optimal bandwidths from the minimisation of the asymptotic mean
        integrated squared error (AMISE).

        Parameters
        ----------
        n_max : int
            Maximum number of points considered in the computation of AMISE.
            NOTE: Computation time of AMISE and its derivatives is quadratic in
                  this number of points.
            (default: 1000)
        lower_bound : float
            Lower bound for bandwidths.
            (default: 1e-10)

        Returns
        -------
        self : active_particles.mkde.MKDE
            MKDE object.
        """

        # RESTRICTION OF DATA

        n_max = int(n_max)
        if n_max >= self.points:
            self.res_data = self.data
        else:
            self.res_data = self.data[random.sample(range(self.points), n_max)] # resticted data sample
        self.res_points, _ = self.res_data.shape                                # size of restricted data sample

        # C EXTENSION

        self.cext = _MKDECExt(self.res_data)    # wrapper of C extension

        # MINIMISAITON ALGORITHM

        self.lower_bound = lower_bound
        self.h0 = (((4/(self.res_points*(self.d + 2)))**(1/(self.d + 4)))
            *np.array([np.std(self.res_data[:, i]) for i in range(self.d)]))    # initial guess for the bandwidths based on Silverman's rule of thumb
        self.minimisation_res = scop.minimize(
            self.cext.AMISE, self.h0, method='trust-constr',                    # minimise self._AMISE with respect to the bandwidths with the Trust-Region Constrained Algorithm
            jac=self.cext.gradAMISE, hess=self.cext.hessAMISE,
            bounds=scop.Bounds([self.lower_bound]*self.d, [np.inf]*self.d)
            , callback=lambda *x: x
            , options={'verbose': 1}
            )
        self.h = self.minimisation_res.x                                        # optimised bandwidths

        return self
