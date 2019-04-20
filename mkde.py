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

from statsmodels.nonparametric.kernel_density import KDEMultivariate

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

    def sm_bw(self, n_max=None, method='cv_ml'):
        """
        Compute optimal bandwidths with the statsmodels package.

        Parameters
        ----------
        n_max : int
            Maximum number of points considered in the computation of AMISE.
            NOTE: Computation time of AMISE and its derivatives is quadratic in
                  this number of points.
            NOTE: if n_max=None then n_max=self.points.
            (default: None)
        method : str
            Type of solver. (see
            https://www.statsmodels.org/stable/generated/statsmodels.nonparametric.kernel_density.KDEMultivariate.html)
            (default: 'cv_ml')

        Returns
        -------
        self : active_particles.mkde.MKDE
            MKDE object.
        """

        # RESTRICTION OF DATA

        self._res_data(n_max)

        # MINIMISAITON ALGORITHM

        self.min_method = ('sm', method)
        self.sm_minimisation_res = KDEMultivariate(self.res_data, 'c'*self.d,
            bw=self.min_method[1])
        self.h = self.sm_minimisation_res.bw    # optimised bandwidths

        return self

    def c_bw(self, n_max=None, method='trust-exact', callback=False):
        """
        Compute optimal bandwidths from the minimisation of the asymptotic mean
        integrated squared error (AMISE) given by the C library.

        Parameters
        ----------
        n_max : int
            Maximum number of points considered in the computation of AMISE.
            NOTE: Computation time of AMISE and its derivatives is quadratic in
                  this number of points.
            NOTE: if n_max=None then n_max=self.points.
            (default: None)
        method : str
            Type of solver. (see
            https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html &
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)
            NOTE: see 'trust-constr' for the constrained version of
                  'trust-exact'.
            NOTE: unconstrained methods (such as 'trust-exact') may result in
                  negative bandwidths.
            (default: 'trust-exact')
        callback : bool
            Print minimisation algorithm state at each iteration.
            (default: False)

        Returns
        -------
        self : active_particles.mkde.MKDE
            MKDE object.
        """

        # RESTRICTION OF DATA

        self._res_data(n_max)

        # C EXTENSION

        self.cext = _MKDECExt(self.res_data)    # wrapper of C extension

        # MINIMISAITON ALGORITHM

        self.min_method = ('c', method)
        self.c_h0 = (((4/(self.res_points*(self.d + 2)))**(1/(self.d + 4)))
            *np.array([np.std(self.res_data[:, i]) for i in range(self.d)]))    # initial guess for the bandwidths based on Silverman's rule of thumb
        self.c_minimisation_res = scop.minimize(
            self.cext.AMISE, self.c_h0, method=self.min_method[1],              # minimise self._AMISE with respect to the bandwidths
            jac=self.cext.gradAMISE, hess=self.cext.hessAMISE,
            bounds=scop.Bounds([0]*self.d, [np.inf]*self.d)
            , callback=lambda *x: (print(x) if callback else None)
            # , options={'verbose': 1}
            )
        self.h = self.c_minimisation_res.x                                      # optimised bandwidths

        return self

    def pdf(self, pdf_points, bw=None):
        """
        Compute probability density function at points pdf_points.

        Parameters
        ----------
        pdf_points : 2D array-like
            Points at which to compute the probability density function.
        bw : 1D array-like
            Bandwidths.
            NOTE: if bw=None then bw=self.h.

        Returns
        -------
        pdf : 1D array-like
            Probability density function at points pdf_points.
        """

        if bw == None: bw = self.h
        return KDEMultivariate(self.data, 'c'*self.d, bw=bw).pdf(pdf_points)

    def _res_data(self, n_max):
        """
        Compute a randomly restricted data sample and store it in sef.res_data
        and its length in self.res_points.

        Parameters
        ----------
        n_max : int
            Maximum number of points considered in the computation of AMISE.
            NOTE: if n_max=None then n_max=self.points.

        Returns
        -------
        self : active_particles.mkde.MKDE
            MKDE object.
        """

        if n_max == None or int(n_max) >= self.points:
            self.res_data = self.data
        else:
            n_max = int(n_max)
            self.res_data = self.data[random.sample(range(self.points), n_max)] # resticted data sample
        self.res_points, _ = self.res_data.shape                                # size of restricted data sample

        return self
