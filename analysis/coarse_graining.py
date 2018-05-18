"""
Module coarse_graining implements a Gaussian coarse-graining adapted from
Illing et al., Phys. Rev. Lett. 117, 208002 (2016) following Goldhirsch and
Goldenberg, Eur. Phys. J. E 9, 245–251 (2002).
"""

import numpy as np

class GaussianCG:
    """
    Gaussian coarse-graining.
    """

    def __init__(self, sigma, r_cut):
        """
        Parameters
        ----------
        sigma : float
            Length scale of Gaussian function.
        r_cut : float
            Coarse-graining cut-off radius.
        """

        self.sigma = sigma  # length scale of Gaussian function
        self.r_cut = r_cut  # coarse-graining cut-off radius

    def function(self, r):
        """
        Parameters
        ----------
        r : float
            Radius.

        Returns
        -------
        phi : float
            Coarse-graining factor at radius r.
        """

        if r > self.r_cut: return 0 # coarse-graining function is zero after cut-off

        Dg = 2*np.pi*(self.sigma**2)*(1 -
            np.exp(-0.5*((self.r_cut/self.sigma)**2)))  # normalisation factor
        return np.exp(-0.5*((r/self.sigma)**2))/Dg      # coarse-graining factor

    def factors(self, positions):
        """
        Parameters
        ----------
        positions : float array
            Coordinates at which coarse-graining is desired.

        Returns
        -------
        CGfactors : Numpy float array
            Coarse-graining factors at positions.
        """

        return np.array(list(map(
            lambda r: self.function(r),
            np.sqrt(np.sum(positions**2, axis=-1))
            ))) # coarse graining factors at positions

class SquareUniformCG:
    """
    Square uniform coarse-graining.
    """

    def __init__(self, dL):
        """
        Parameters
        ----------
        dL : float
            Length of square box on which to average.
        """

        self.dL = dL    # averaging square length

    def function(self, position):
        """
        Parameters
        ----------
        position : float array
            Coordinates.

        Returns
        -------
        phi : float
            Coarse-graining factor at position position.
        """

        if (np.abs(np.array(position)) > dL/2).any(): return 0  # coarse-graining function is zero outside square
        return 1                                                # is one in

    def factors(self, positions):
        """
        Parameters
        ----------
        positions : float array
            Coordinates at which coarse-graining is desired.

        Returns
        -------
        CGfactors : Numpy float array
            Coarse-graining factors at positions.
        """

        CGfactors = np.array(list(map(
            lambda position:
            self.function(position),
            positions
            )))
        sumCGfactors = np.sum(CGfactors)
        if np.sum(CGfactors) == 0: return 0
        return CGfactors/sumCGfactors   # coarse graining factors at positions

class CoarseGraining:
    """
    Enables unique calculation of coarse-graining factors and then calculation
    of coarse-graining avergages.
    """

    def __init__(self, factors_function, positions):
        """
        Parameters
        ----------
        factors_function : function
            Function of array of coordinates which returns coarse-graining
            factors at these coordinates.
        positions : float array
            Coordinates at which coarse-graining is desired.
        """

        self.CGfactors = np.array(factors_function(positions))  # coarse-graining factors at positions

    def average(self, var):
        """
        Coarse-graining averaging.

        Parameters
        ----------
        var : float array
            Values of variable to coarse-grain at different positions from
            point at which coarse-graining is desired.

        Returns
        -------
        average : float
            Coarse-grained variable.
        """

        return np.sum(
            np.transpose(np.array(self.CGfactors,
            ndmin=len(np.array(var).shape)))
            *np.array(var), axis=0) # coarse-grained variable
