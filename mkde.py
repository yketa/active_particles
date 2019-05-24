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

# from multiprocessing import Pool

from KDEpy import FFTKDE

from statsmodels.nonparametric.kernel_density import KDEMultivariate

# C EXTENSION

_dir_path = os.path.dirname(os.path.realpath(__file__)) # script directory path

_c_ext_so_path = os.path.join(_dir_path, 'cmkde.so')    # path to C extension shared object .so file
_c_ext_c_path = os.path.join(_dir_path, 'cmkde.c')      # path to C extension source .c file
_c_ext_h_path = os.path.join(_dir_path, 'cmkde.h')      # path to C extension source .h file

def compile():
    """
    Compile C extension.
    """

    subprocess.run(['gcc', '-shared', '-o', _c_ext_so_path,
        '-std=c99', '-fPIC', _c_ext_c_path],
        cwd=_dir_path)

if (not(os.path.isfile(_c_ext_so_path))                                     # C extention shared object does not exist
    or os.path.getctime(_c_ext_c_path) > os.path.getctime(_c_ext_so_path)   # C extension source .c file is more recent than shared object
    or os.path.getctime(_c_ext_h_path) > os.path.getctime(_c_ext_so_path)): # C extension source .h file is more recent than shared object
    compile()                                                               # compile C extension

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

    def grid_pdf(self, bw=None):
        """
        Compute grid of probability density function with KDEpy.FFTKDE and a
        Gaussian kernel.

        NOTE: The probability density function returned by KDEpy.FFTKDE seems
              not to be normalised for multivariate variables, this function
              should thus be considered with extra care.

        Parameters
        ----------
        bw : 1D array-like
            Bandwidths.
            NOTE: if bw=None then bw=self.h.

        Returns
        -------
        x : 2D array-like
            Coordinates at which the probability density function is evaluated.
        y : array-like
            Values of the probability density function.
        """

        if bw == None: bw = self.h
        bw = np.array(bw, ndmin=2)
        data = self.data/bw # rescaled data
        x, y = FFTKDE(kernel='gaussian', bw=1).fit(data).evaluate()
        return x*bw, y

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

    # DEPRECATED
    # AMISE IS NOW COMPUTED FROM THE C LIBRARY

    def _AMISE(self, h):
        """
        Compute asymptotic mean integrated squared error (AMISE).

        Parameters
        ----------
        h : 1D Numpy array
            Array of bandwidths.

        Returns
        -------
        AMISE : float
            Asymptotic mean integrated squared error (AMISE).
        """

        rdelta = self.delta/h    # rescaled different between data points

        # with Pool(processes=self.processes) as pool:    # pool of worker processes for multiprocessing
        #
        #     b = [pool.apply_async(self._B, args=(rd,))
        #         for rd in rdelta]
        #     c = [pool.apply_async(self._C, args=(rd,))
        #         for rd in rdelta]
        #
        #     B = np.array([_.get() for _ in b])
        #     C = np.array([_.get() for _ in c])

        B = np.array(list(map(
            self._B,
            rdelta)))
        C = np.array(list(map(
            self._C,
            rdelta)))

        sumBC = np.sum(B*C)

        AMISE = self._A(h)*(
            1/(((2*np.sqrt(np.pi))**self.d)*self.res_points)
            + sumBC/(2*self.res_points*(self.res_points - 1)))
        return AMISE

    def _AMISE_grad(self, h):
        """
        Compute the gradient of AMISE with respect to bandwidths.

        Parameters
        ----------
        h : 1D Numpy array
            Array of bandwidths.

        Returns
        -------
        AMISE_grad : 1D float Numpy array of length self.d
            Gradient of AMISE.
        """

        AMISE_grad = np.array([self._pAMISE(h, p)
            for p in range(self.d)])
        return AMISE_grad

    def _pAMISE(self, h, p):
        """
        Compute first derivative of AMISE with respect to the p-th bandwidth
        h_p.

        Parameters
        ----------
        h : 1D Numpy array
            Array of bandwidths.
        p : int
            Index of bandwidth for which to compute the derivative.

        Returns
        -------
        pAMISE : float
            \partial AMISE/\partial h_p (h).
        """

        rdelta = self.delta/h    # rescaled different between data points

        # with Pool(processes=self.processes) as pool:    # pool of worker processes for multiprocessing
        #
        #     b = [pool.apply_async(self._B, args=(rd,))
        #         for rd in rdelta]
        #     c = [pool.apply_async(self._C, args=(rd,))
        #         for rd in rdelta]
        #     pb = [pool.apply_async(self._pB, args=(h, p, rd))
        #         for rd in rdelta]
        #     pc = [pool.apply_async(self._pC, args=(h, p, rd))
        #         for rd in rdelta]
        #
        #     B = np.array([_.get() for _ in b])
        #     C = np.array([_.get() for _ in c])
        #     pB = np.array([_.get() for _ in pb])
        #     pC = np.array([_.get() for _ in pc])

        B = np.array(list(map(
            self._B,
            rdelta)))
        C = np.array(list(map(
            self._C,
            rdelta)))
        pB = np.array(list(map(
            lambda rd: self._pB(h, p, rd),
            rdelta)))
        pC = np.array(list(map(
            lambda rd: self._pC(h, p, rd),
            rdelta)))

        sumBC = np.sum(B*C)
        sumpBCBpC = np.sum(pB*C + B*pC)

        pAMISE = (
            self._pA(h, p)*(
                1/(((2*np.sqrt(np.pi))**self.d)*self.res_points)
                + sumBC/(2*self.res_points*(self.res_points - 1)))
            + self._A(h)*sumpBCBpC/(2*self.res_points*(self.res_points - 1)))
        return pAMISE

    def _AMISE_hess(self, h):
        """
        Compute the Hessian matrix of AMISE with respect to bandwidths.
        NOTE: the Hessian matrix is symmetric.

        Parameters
        ----------
        h : 1D Numpy array
            Array of bandwidths.

        Returns
        -------
        AMISE_hess : 2D float Numpy array of length self.d x self.d
            Hessian matrix of AMISE.
        """

        AMISE_hess = np.zeros((self.d, self.d))
        for p in range(self.d):
            for q in range(p, self.d):
                pqAMISE = self._pqAMISE(h, p, q)
                AMISE_hess[p, q] = pqAMISE
                AMISE_hess[q, p] = pqAMISE
        return AMISE_hess

    def _pqAMISE(self, h, p, q):
        """
        Compute second cross-derivative of AMISE with respect to the p-th and
        q-th bandwidths h_p and h_q.

        Parameters
        ----------
        h : 1D Numpy array
            Array of bandwidths.
        p : int
            Index of the first bandwidth for which to compute the derivative.
        q : int
            Index of the second bandwidth for which to compute the derivative.

        Returns
        -------
        pqAMISE : float
            \partial^2 AMISE/\partial h_p \partial h_q (h).
        """

        rdelta = self.delta/h    # rescaled different between data points

        # with Pool(processes=self.processes) as pool:    # pool of worker processes for multiprocessing
        #
        #     b = [pool.apply_async(self._B, args=(rd,))
        #         for rd in rdelta]
        #     c = [pool.apply_async(self._C, args=(rd,))
        #         for rd in rdelta]
        #     pb = [pool.apply_async(self._pB, args=(h, p, rd))
        #         for rd in rdelta]
        #     pc = [pool.apply_async(self._pC, args=(h, p, rd))
        #         for rd in rdelta]
        #     qb = [pool.apply_async(self._pB, args=(h, q, rd))
        #         for rd in rdelta]
        #     qc = [pool.apply_async(self._pC, args=(h, q, rd))
        #         for rd in rdelta]
        #     pqb = [pool.apply_async(self._pqB, args=(h, p, q, rd))
        #         for rd in rdelta]
        #     pqc = [pool.apply_async(self._pqC, args=(h, p, q, rd))
        #         for rd in rdelta]
        #
        #     B = np.array([_.get() for _ in b])
        #     C = np.array([_.get() for _ in c])
        #     pB = np.array([_.get() for _ in pb])
        #     pC = np.array([_.get() for _ in pc])
        #     qB = np.array([_.get() for _ in qb])
        #     qC = np.array([_.get() for _ in qc])
        #     pqB = np.array([_.get() for _ in pqb])
        #     pqC = np.array([_.get() for _ in pqc])

        B = np.array(list(map(
            self._B,
            rdelta)))
        C = np.array(list(map(
            self._C,
            rdelta)))
        pB = np.array(list(map(
            lambda rd: self._pB(h, p, rd),
            rdelta)))
        pC = np.array(list(map(
            lambda rd: self._pC(h, p, rd),
            rdelta)))
        qB = (pB if p == q else
            np.array(list(map(
                lambda rd: self._pB(h, q, rd),
                rdelta))))
        qC = (pC if p == q else
            np.array(list(map(
                lambda rd: self._pC(h, q, rd),
                rdelta))))
        pqB = np.array(list(map(
            lambda rd: self._pqB(h, p, q, rd),
            rdelta)))
        pqC = np.array(list(map(
            lambda rd: self._pqC(h, p, q, rd),
            rdelta)))

        sumBC = np.sum(B*C)
        sumpBCBpC = np.sum(pB*C + B*pC)
        sumqBCBqC = np.sum(qB*C + B*qC)
        sumpqBCpBqCqBpCBpqC = np.sum(pqB*C + pB*qC + qB*pC + B*pqC)

        pqAMISE = (
            self._pqA(h, p, q)*(
                1/(((2*np.sqrt(np.pi))**self.d)*self.res_points)
                + sumBC/(2*self.res_points*(self.res_points - 1)))
            + (self._pA(h, p)*sumqBCBqC + self._pA(h, q)*sumpBCBpC
                + self._A(h)*sumpqBCpBqCqBpCBpqC)
                /(2*self.res_points*(self.res_points - 1)))
        return pqAMISE

    def _A(self, h):
        """
        Parameters
        ----------
        h : 1D Numpy array
            Array of bandwidths.
        """

        return 1/np.prod(h)

    def _pA(self, h, p):
        """
        Parameters
        ----------
        h : 1D Numpy array
            Array of bandwidths.
        p : int
            Index of the first bandwidth for which to compute the derivative.
        """

        return -self._A(h)/h[p]

    def _pqA(self, h, p, q):
        """
        Parameters
        ----------
        h : 1D Numpy array
            Array of bandwidths.
        p : int
            Index of the first bandwidth for which to compute the derivative.
        q : int
            Index of the second bandwidth for which to compute the derivative.
        """

        return (1 + (p==q))*self._A(h)/(h[p]*h[q])

    def _B(self, rdelta):
        """
        Parameters
        ----------
        rdelta : 1D Numpy array
            Difference between a couple of data points, rescaled by the
            corresponding bandwidths.
        """

        sumsq = np.sum(rdelta**2)

        return sumsq**2 - (2*self.d + 4)*sumsq + (self.d**2 + 2*self.d)

    def _pB(self, h, p, rdelta):
        """
        Parameters
        ----------
        rdelta : 1D Numpy array
            Difference between a couple of data points, rescaled by the
            corresponding bandwidths.
        h : 1D Numpy array
            Array of bandwidths.
        p : int
            Index of the first bandwidth for which to compute the derivative.
        """

        sumsq = np.sum(rdelta**2)

        return (rdelta[p]**2)*(2*(2*self.d + 4) - 4*sumsq)/h[p]

    def _pqB(self, h, p, q, rdelta):
        """
        Parameters
        ----------
        rdelta : 1D Numpy array
            Difference between a couple of data points, rescaled by the
            corresponding bandwidths.
        h : 1D Numpy array
            Array of bandwidths.
        p : int
            Index of the first bandwidth for which to compute the derivative.
        q : int
            Index of the second bandwidth for which to compute the derivative.
        """

        sumsq = np.sum(rdelta**2)

        return (8*(rdelta[p]**2)*(rdelta[q]**2)/(h[p]*h[q])
            + (p==q)*(rdelta[p]**2)*(12*sumsq - 6*(2*self.d + 4))/(h[p]**2))

    def _C(self, rdelta):
        """
        Parameters
        ----------
        rdelta : 1D Numpy array
            Difference between a couple of data points, rescaled by the
            corresponding bandwidths.
        """

        return np.prod(np.exp(-(rdelta**2)/2)/np.sqrt(2*np.pi))

    def _pC(self, h, p, rdelta):
        """
        Parameters
        ----------
        rdelta : 1D Numpy array
            Difference between a couple of data points, rescaled by the
            corresponding bandwidths.
        h : 1D Numpy array
            Array of bandwidths.
        p : int
            Index of the first bandwidth for which to compute the derivative.
        """

        return (rdelta[p]**2)*self._C(rdelta)/h[p]

    def _pqC(self, h, p, q, rdelta):
        """
        Parameters
        ----------
        rdelta : 1D Numpy array
            Difference between a couple of data points, rescaled by the
            corresponding bandwidths.
        h : 1D Numpy array
            Array of bandwidths.
        p : int
            Index of the first bandwidth for which to compute the derivative.
        q : int
            Index of the second bandwidth for which to compute the derivative.
        """

        return ((rdelta[q]**2)/h[q] - 3*(p==q)/h[p])*self._pC(h, p, rdelta)
