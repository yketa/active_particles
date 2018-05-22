"""
Module maths provides useful mathematic tools.
"""

import numpy as np

def relative_positions(positions, point, box_size):
    """
    Returns relative positions to point in box of extent
    (-box_size/2, box_size) in both dimensions of space.

    Parameters
    ----------
    positions : float array
        Position of single point or array of positions.
    point : float array
        Position of the new centre.
    box_size : float or array
        Length of the box in one dimension or all dimensions.

    Returns
    -------
    rel_positions : float array
        Relative positions.
    """

    return (np.array(positions) - np.array(point)
        + np.array(box_size)/2)%np.array(box_size) - np.array(box_size)/2

def wo_mean(arr):
    """
    Returns deviation of values in array with respect to mean of array.

    Parameters
    ----------
    arr : array like
        Array of values.

    Returns
    -------
    dev_arr : array like
        Deviations from mean of array.
    """

    return np.array(arr) - np.mean(arr, axis=0)

class DictList(dict):
    """
    Custom hash table class to give value [] to uninitialised keys.
    """
    def __init__(self):
        super().__init__()
    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            return []

def g2Dto1D(g2D, L):
    """
    Returns cylindrical average of 2D grid.

    Parameters
    ----------
    g2D : 2D array
        2D grid.
        NOTE: g2D[0, 0] is considered the r=0 point on the grid, and we
        consider periodic boundaries.
    L : float or float array
        Length of the box represented by the grid in one dimension or all
        dimensions.

    Returns
    -------
    g1D : Numpy array
        Array of (r, g1D(r)) with g1D(r) the averaged 2D grid at radius r.
    """

    g2D = np.array(g2D)
    dL = np.array(L)/np.array(g2D.shape)    # boxes separation in each direction
    r_max = np.min(L)/2                     # maximum radius to be calculated in number of boxes

    g1D_dic = DictList()    # hash table of radii and values at radii

    for i in range(g2D.shape[0]):
        for j in range(g2D.shape[1]):
            radius = np.sqrt(np.sum((np.array((i, j))*dL)**2))  # radius corresponding to coordinates [i, j], [-i, j], [i, -j], [-i, -j]
            if radius <= r_max:
                g1D_dic[radius] += [g2D[i, j], g2D[-i, j], g2D[i, -j],
                    g2D[-i, -j]]

    return np.array(list(map(
        lambda radius: [radius, np.mean(g1D_dic[radius])],
        sorted(g1D_dic))))

def g2Dto1Dsquare(g2D, L):
    """
    Returns cylindrical average of square 2D grid.

    Parameters
    ----------
    g2D : 2D array
        Square 2D grid.
        NOTE: g2D[0, 0] is considered the r=0 point on the grid, and we
        consider periodid boundaries.
    L : float
        Length of the box represented by the grid in one dimension.

    Returns
    -------
    g1D : Numpy array
        Array of (r, g1D(r)) with g1D(r) the averaged 2D grid at radius r.
    """

    g2D = np.array(g2D)
    dL = L/g2D.shape[0]     # boxes separation in each direction
    r_max = g2D.shape[0]/2  # maximum radius to be calculated in number of boxes

    g1D_dic = DictList()    # hash table of radii and values at radii

    for i in range(g2D.shape[0]):
        for j in range(g2D.shape[1]):
            sqradius = i**2 + j**2  # radius corresponding to coordinates [i, j], [-i, j], [i, -j], [-i, -j]
            if sqradius <= r_max**2:
                g1D_dic[sqradius] += [g2D[i, j], g2D[-i, j], g2D[i, -j],
                    g2D[-i, -j]]

    return np.array(list(map(
        lambda sqradius: [dL*np.sqrt(sqradius), np.mean(g1D_dic[sqradius])],
        sorted(g1D_dic))))

def normalise1D(vector):
    """
    Returs 1D vector of unitary norm with same direction.

    Parameters
    ----------
    vector : 1D array-like
        Vector to normalise.

    Returns
    -------
    u_vector : 1D Numpy array
        Unitary vector with same direction.
    """

    norm = np.linalg.norm(vector)           # vector norm
    if norm == 0: return np.array(vector)   # vector is 0
    return np.array(vector)/norm
