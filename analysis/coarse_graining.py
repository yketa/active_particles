"""
Module coarse_graining implements a Gaussian coarse-graining adapted from
Illing et al., Phys. Rev. Lett. 117, 208002 (2016) following Goldhirsch and
Goldenberg, Eur. Phys. J. E 9, 245–251 (2002).
"""

import numpy as np

def gaussianCG_function(r, sigma, r_cut):
    """
    Gaussian coarse-graining function.

    Parameters
    ----------
    r : float
        Radius.
    sigma : float
        Length scale of Gaussian function.
    r_cut : float
        Coarse-graining cut-off radius.

    Returns
    -------
    phi : float
        Coarse-graining factor at radius r.
    """

    if r > r_cut: return 0  # coarse-graining function is zero after cut-off

    Dg = 2*np.pi*(sigma**2)*(1 - np.exp(-0.5*((r_cut/sigma)**2)))   # normalisation factor
    return np.exp(-0.5*((r/sigma)**2))/Dg                           # coarse-graining factor

def gaussianCG_factors(radii, sigma, r_cut):
    """
    Gaussian coarse-graining factors.

    Parameters
    ----------
    var : unidimensional iterable of length N containing floats
        Values of variable to coarse-grain at N different radii from point
        at which coarse-graining is desired.
    radii : unidimensional iterable of length N containing floats
        Radii of points from the one at which coarse-graining is desired.
    sigma : float
        Length scale of Gaussian function.
    r_cut : float
        Coarse-graining cut-off radius.

    Returns
    -------
    CGfactors : Numpy array of length N containing floats
        Coarse-graining factors at radii.
    """

    CGfactors = np.array(list(map(
        lambda r: gaussianCG_function(r, sigma, r_cut),
        radii
        )))                         # coarse graining factors at radii
    return CGfactors

def gaussianCG_average(var, CGfactors):
    """
    Gaussian coarse-graining.

    Parameters
    ----------
    var : unidimensional iterable of length N containing floats
        Values of variable to coarse-grain at N different radii from point
        at which coarse-graining is desired.
    CGfactors : unidimensional iterable of length N containing floats
        Coarse-graining factors at radii.

    Returns
    -------
    average : float
        Coarse-grained variable.
    """

    return np.sum(np.array(var)*np.array(CGfactors))    # coarse-grained variable

class CoarseGraining:
    """
    Enables unique calculation of coarse-graining factors and then calculation
    of coarse-graining avergages.
    """

    def __init__(self, radii, sigma, r_cut):
        """
        Parameters
        ----------
        radii : unidimensional iterable of length N containing floats
            Radii of points from the one at which coarse-graining is desired.
        sigma : float
            Length scale of Gaussian function.
        r_cut : float
            Coarse-graining cut-off radius.
        """

        self.CGfactors = gaussianCG_factors(radii, sigma, r_cut)    # coarse graining factors at radii

    def average(self, var):
        """
        Gaussian coarse-graining.

        Parameters
        ----------
        var : unidimensional iterable of length N containing floats
            Values of variable to coarse-grain at N different radii from point
            at which coarse-graining is desired.

        Returns
        -------
        average : float
            Coarse-grained variable.
        """

        return gaussianCG_average(var, self.CGfactors)  # coarse-grained variable
