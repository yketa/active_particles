"""
Module scde provides methods to compute self-consistent density estimations.
(see
https://yketa.github.io/UBC_2018_Wiki/#Self-consistent%20density%20estimation)
"""

import numpy as np

from scipy.interpolate import griddata

from fastkde import fastKDE

from multiprocessing import Pool

from itertools import product
from functools import partial

from operator import itemgetter

from active_particles.maths import Wrap

class PDF:
    """
    Compute and manipulate probability density functions computed via
    self-consistent density estimation.
    """

    def __init__(self, *vars, renormalise=True,
        wrap_period=None, wrap_method='linear', wrap_fill_value=0,
        wrap_processes=None, **fastKDE_kwargs):
        """
        Compute probability density function.
        (see fastkde.fastKDE.pdf)

        NOTE: Coordinates in self.axes are in the same order as the input
              variables while it is in reversed order in self.pdf (see
              fastkde.fastKDE.pdf).
              See active_particles.scde.PDF.evaluate for probability density
              function evaluation.

        Positional arguments
        --------------------
        vars : array-like
            Input variables.

        Parameters
        ----------
        renormalise : bool
            Rescale probability density function values by the integral over
            the computed volume.
            DEFAULT: True
        wrap_period : float
            Period over which to wrap the computed probability density function.
            NOTE: If wrap_period == None, the computed probability density
                  function remains unwrapped.
            DEFAULT: None
        wrap_method : string
            Method of interpolation. (see scipy.interpolate.griddata)
            DEFAULT: linear
        wrap_fill_value : float
            Value used to fill in for requested points outside of the convex
            hull of the input points. (see scipy.interpolate.griddata)
            DEFAULT: 0
        wrap_processes : int
            Number of worker processes to use. (see multiprocessing.Pool)
            NOTE: If processes == None then processes = os.cpu_count().
            DEFAULT: None

        Optional keyword arguments
        --------------------------
        (see fastkde.fastKDE.pdf)
        """

        self.vars = vars
        self.n = len(self.vars)
        self.fastKDE_kwargs = fastKDE_kwargs

        self.pdf, self.axes = fastKDE.pdf(*self.vars, **self.fastKDE_kwargs)
        if self.n == 1: self.axes = [np.array(self.axes)]
        self._extended_axes()

        if wrap_period != None: self.wrap(wrap_period,
            method=wrap_method, fill_value=wrap_fill_value,
            processes=wrap_processes)

        if renormalise: self.renormalise()

    def wrap(self, p, method='linear', fill_value=0, processes=None):
        """
        Wrap self.pdf so that it is p-periodic and evaluates it at
        self.pdf.shape linearly spaced points.

        NOTE: self.pdf is interpolated.

        NOTE: As a matter of efficiency, wrapped function is evaluated with a
              pool of worker processes. (see  multiprocessing.Pool and
              multiprocessing.Pool.starmap)

        Parameters
        ----------
        p : float
            Period of the function in each direction.
        method : string
            Method of interpolation. (see scipy.interpolate.griddata)
            DEFAULT: linear
        fill_value : float
            Value used to fill in for requested points outside of the convex
            hull of the input points. (see scipy.interpolate.griddata)
            DEFAULT: 0
        processes : int
            Number of worker processes to use. (see multiprocessing.Pool)
            NOTE: If processes == None then processes = os.cpu_count().
            DEFAULT: None
        """

        old_pdf_flat, old_extended_axes_flat = self._flat()
        wrap = Wrap(old_extended_axes_flat, old_pdf_flat, [p]*self.n)

        self.axes = [np.linspace(-p/2, p/2, len(self.axes[i]))
            for i in range(self.n)]
        self._extended_axes()
        self.pdf = np.transpose(
            np.reshape(
                wrap.evaluate(
                    *itemgetter(*np.ndindex(self.pdf.shape[::-1]))
                        (self.extended_axes),
                    method=method, fill_value=fill_value, processes=processes),
                self.pdf.shape[::-1]))

    def evaluate(self, *coordinates, method='linear', fill_value=0,
        processes=None):
        """
        Evaluate interpolated probability density function from evaluated
        points.

        NOTE: As a matter of efficiency, probability density function at
              coordinates is evaluated with a pool of worker processes. (see
              multiprocessing.Pool and multiprocessing.Pool.starmap)

        Parameters
        ----------
        method : string
            Method of interpolation. (see scipy.interpolate.griddata)
            DEFAULT: linear
        fill_value : float
            Value used to fill in for requested points outside of the convex
            hull of the input points. (see scipy.interpolate.griddata)
            DEFAULT: 0
        processes : int
            Number of worker processes to use. (see multiprocessing.Pool)
            NOTE: If processes == None then processes = os.cpu_count().
            DEFAULT: None

        Positional arguments
        --------------------
        coordinates : (self.n,) array-like
            Coordinates at which to evaluate the interpolated probability
            density function.

        Returns
        -------
        pdfs : (len(coordinates),) array-like
            Interpolated probability density function at coordinates.
        """

        pdf_flat, extended_axes_flat = self._flat()

        with Pool(processes=processes) as pool: # pool of worker processes
            pdfs = pool.starmap(
                partial(griddata, method=method, fill_value=fill_value),
                product([extended_axes_flat], [pdf_flat],
                    np.array(coordinates)[:, ::-1]))
        return np.array(pdfs).flatten()

    def integrate(self, apply_func=None):
        """
        Integrate apply_func(pdf) over the whole computed volume.

        Parameters
        ----------
        apply_func : function
            Function to apply to the probability density function.
            NOTE: If apply_func == None, then apply_func = lambda x: x.
            DEFAULT: None

        Returns
        -------
        integral : float
            Integrated probability density function.
        """

        if apply_func == None: apply_func = lambda x: x

        integral = apply_func(self.pdf)
        for ax in range(self.n):
            integral = np.trapz(integral, self.axes[-ax-1], axis=0)
        return integral

    def entropy(self):
        """
        Returns the Shannon entropy of the distribution.

        Returns
        -------
        entropy : float
            Shannon entropy of the distribution.
        """

        return self.integrate(apply_func=
            lambda x: - x * np.ma.log(x).filled(0))

    def renormalise(self):
        """
        Rescale probability density function values self.pdf by the integral
        over the the computed volume.
        """

        self.pdf /= self.integrate(apply_func=lambda x: x)

    def _extended_axes(self):
        """
        Sets self.extended_axes as self.axes coordinates in extended form.
        """

        self.extended_axes = np.empty(self.pdf.shape + (self.n,))   # extended points coordinates at which the probability density function is evaluated
        for ind in np.ndindex(self.pdf.shape):
            self.extended_axes[ind] = [self.axes[-i-1][ind[i]]
                for i in range(self.n)]

    def _flat(self):
        """
        Returns flat self.pdf and self.extended_axes.

        Returns
        -------
        pdf_flat : (np.prod(self.pdf.shape),) Numpy array
            Flat self.pdf.
        extended_axes_flat : (np.prod(self.pdf.shape), self.n) Numpy array
            Flat self.extended_axes.
        """

        pdf_flat = self.pdf.flatten()
        extended_axes_flat = self.extended_axes.reshape((-1, len(self.axes)))

        return pdf_flat, extended_axes_flat
