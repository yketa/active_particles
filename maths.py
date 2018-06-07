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

def g2Dto1Dgrid(g2D, grid, average_grid=False):
    """
    Returns cylindrical average of square 2D grid with values of radius given
    by other parameter grid.

    Parameters
    ----------
    g2D : 2D array
        Square 2D grid.
    grid : 2D array
        Array of radii.
    average_grid : bool
        Return g2D grid with cylindrically averaged values.

    Returns
    -------
    g1D : Numpy array
        Array of (r, g1D(r)) with g1D(r) the averaged 2D grid at radius r.
    g2D_cylindrical [average_grid] : Numpy array
        Cylindrically averaged g2D.
    """

    g2D = np.array(g2D)
    grid = np.array(grid)

    g1D_dic = DictList()    # hash table of radii and values at radii

    for i in range(g2D.shape[0]):
        for j in range(g2D.shape[1]):
            g1D_dic[grid[i, j]] += [g2D[i, j]]

    g1D = np.array(list(map(
        lambda radius: [radius, np.mean(g1D_dic[radius])],
        sorted(g1D_dic))))

    if not(average_grid): return g1D

    g2D_cylindrical = np.zeros(grid.shape)
    for radius, mean_g in zip(*np.transpose(g1D)):
        for i, j in zip(*np.where(grid == radius)):
            g2D_cylindrical[i, j] = mean_g

    return g1D, g2D_cylindrical

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

def amplogwidth(arr, factor=2):
    """
    Calculates the amplitudes of elements in array arr and, excluding the
    zeros, returns the mean of the logarithms of these amplitudes plus and
    minus factor times their standard deviation.

    Parameters
    ----------
    arr : array like
        Array.
    factor : float
        Width factor. (default: 2)

    Returns
    -------
    min : float
        E(log(||arr||)) - factor*V(log(||arr||))
    max : float
        E(log(||arr||)) + factor*V(log(||arr||))
    """

    log = np.ma.log10(np.sqrt(np.sum(arr**2, axis=-1))) # logarithms of amplitudes
    mean = log.mean()                                   # means of logarithms of amplitudes
    std = log.std()                                     # standard deviation of logarithms of amplitudes

    return mean - factor*std, mean + factor*std

def mean_sterr(values):
    """
    Returns mean and standard error of values.

    Parameters
    ----------
    values : float array
        Values.

    Returns
    -------
    mean : float
        Mean of values.
    sterr : float
        Standard error of values.
    """

    values = np.array(values)
    if values.size == 0: return None, None

    return np.mean(values), np.std(values)/np.sqrt(np.prod(values.shape))

class Grid:
    """
    Manipulate 2D grids, in which we consider the values to correspond to a
    variable at uniformly distributed positions in space.
    """

    def __init__(self, grid, extent=(-1, 1, -1, 1)):
        """
        Sets the grid.

        Parameters
        ----------
        grid : array-like
            2D grid.
        extent : scalars (left, right, bottom, top)
            Values of space variables at corners. (default: (-1, 1, -1, 1))
        """

        self.grid = np.array(grid)
        self.shape = self.grid.shape    # shape of the grid

        self.extent = extent
        self.sep_boxes_x = (self.extent[1] - self.extent[0])/self.grid.shape[0] # distance between consecutive boxes in x (first) direction
        self.sep_boxes_y = (self.extent[3] - self.extent[2])/self.grid.shape[1] # distance between consecutive boxes in y (second) direction

    def __getitem__(self, *key):
        """
        Associates Grid[key] to Grid.grid[key].

        Parameters
        ----------
        key : *
            Key to access.

        Returns
        -------
        value : *
            Grid.grid[key]
        """

        return self.grid.__getitem__(*key)

    def in_grid(self, x, y):
        """
        Indicates if point (x, y) in cartesian coordinates is in grid.

        Parameters
        ----------
        x : float
            x-coordinate
        y : float
            y-coordinate

        Returns
        -------
        is_in_grid : bool
            (x, y) in grid.
        """

        return (x >= self.extent[0] and x <= self.extent[1]
            and y >= self.extent[2] and y <= self.extent[3])

    def get_value_cartesian(self, x, y):
        """
        Get value of grid at position in cartesian coordinates.

        Parameters
        ----------
        x : float
            x-coordinate
        y : float
            y-coordinate

        Returns
        -------
        value : *
            Value at (x, y).
        """

        if not(self.in_grid(x, y)): return None # point not on grid

        index_x = int((x - self.extent[0])//self.sep_boxes_x)%self.shape[0] # index correponding to x
        index_y = int((y - self.extent[2])//self.sep_boxes_y)%self.shape[1] # index correponding to y

        return self.grid[index_x, index_y]

    def get_value_polar(self, r, angle, centre=(0, 0)):
        """
        Get value of grid at position in polar coordinates.

        Parameters
        ----------
        r : float
            Radius from centre.
        angle : float
            Angle from x-direction.
        centre : float tuple
            Origin for calculation. (default: (0, 0))

        Returns
        -------
        value : *
            Value at (r, angle) from centre.
        """

        x = centre[0] - r*np.sin(angle) # corresponding cartesian x-coordinate
        y = centre[1] + r*np.cos(angle) # corresponding cartesian y-coordinate

        return self.get_value_cartesian(x, y)

def vector_vector_grid(vector1, vector2):
    """
    From vector1 = (v1_i)_i and vector2 = (v2_i)_i, returns matrix
    M = (M_{i, j})_{i, j} = ((v1_i, v2_j))_{i, j}.

    Parameters
    ----------
    vector1 : 1D array-like
        Vector 1.
    vector2 : 1D array-like
        Vector 2.

    Returns
    -------
    M : 2D array-like
        Matrix M.
    """

    M = np.zeros((len(vector1), len(vector2), 2))
    M[:, :, 0] = vector1
    M = np.transpose(M, (1, 0, 2))
    M[:, :, 1] = vector2

    return M

def kFFTgrid(grid, d=1):
    """
    Calculates the Fast Fourier Transform (FFT) of 2D grid and returns its dot
    and cross product with corresponding wave vector as well as grid of wave
    vector.

    Parameters
    ----------
    grid : array-like
        2D grid of 2D vectors (i.e., (_, _, 2) grid).
    d : float
        Distance between two consecutive grid boxes. (default: 1)

    Returns
    -------
    wave_vectors : (*grid.shape, 2) Numpy array
        Grid of wave vectors.
    k_cross_dot_grid : (*grid.shape, 2) Numpy array
        Concatenanated grids of cross and dot products between wave vectors and
        grid Fourier transform.
        NOTE: Concatenation is for dimension reasons.
    """

    FFTgrid = np.fft.fft2(grid, axes=(0, 1))    # Fourier transform of grid
    wave_vectors = 2*np.pi*vector_vector_grid(
        np.fft.fftfreq(grid.shape[0], d=d),
        np.fft.fftfreq(grid.shape[1], d=d))     # grid of wave vectors

    k_cross_grid = np.cross(wave_vectors, FFTgrid)  # k cross FFTgrid

    k_dot_grid = np.zeros(FFTgrid.shape[:2], dtype=np.complex128)
    for i in range(FFTgrid.shape[0]):
        for j in range(FFTgrid.shape[1]):
            k_dot_grid[i, j] = np.dot(wave_vectors[i, j], FFTgrid[i, j])    # k dot FFTgrid

    return wave_vectors, np.concatenate(
        (np.reshape(k_cross_grid, k_cross_grid.shape + (1,)),
        np.reshape(k_dot_grid, k_dot_grid.shape + (1,))),
        axis=-1)

def divide_arrays(array1, array2):
    """
    Divide array1 by array2, and outputs 0 values where array2 is equal to 0.
    NOTE: array1, array2 and out must have the same shapes.

    Parameters
    ----------
    array1 : array-like
        Numerator array.
    array2 : array-like
        Denominator array.

    Returns
    -------
    array : array-like
        Quotient array.
    """

    array2 = np.array(array2)

    return np.divide(array1, array2,
        out=np.zeros(array2.shape), where=array2!=0)

def grid_from_function(grid_values, function, dimension=None):
    """
    Returns grid of dimension dimension from function evaluated at
    grid_values.

    NOTE: in 1D, use list(map(function, grid_values)).

    Parameters
    ----------
    grid_values : array-like
        Grid of values at which to evaluate function.
    function : function
        Function of grid values variables.
    dimension : int or None
        Dimension of the grid to return. (default: None)
        NOTE: None is considered as dimension equal to dimension of
        grid_values.
        NOTE : dimension can be lesser than grid_values dimension, in this case
        values in grid_values remaining dimensions are passed as positional
        parameters to function.

    Returns
    -------
    grid : Numpy array
        Grid created from function evaluated at grid_values.
    """

    grid_values = np.array(grid_values)
    grid_values_dimension = len(grid_values.shape)  # dimension of grid_values

    if dimension == None or dimension > grid_values_dimension:
        dimension = grid_values_dimension

    grid_shape = grid_values.shape[:dimension]  # shape of grid
    values_length = np.prod(grid_shape)         # number of elements in grid

    grid = list(map(
        lambda value: function(*np.array(value).flatten()),
        np.reshape(grid_values,
        (values_length, *grid_values.shape[dimension:]))))

    return np.reshape(grid, grid_shape)

def step_function(X, Y):
    """
    Returns step function f from array-likes X and Y, such that
           | Y[0]      if x <= X[0]
    f(x) = | Y[i + 1]  if X[i] < x <= X[i + 1]
           | Y[-1]     if x >= X[-1]

    NOTE: X and Y must have the same shape.

    Parameters
    ----------
    X : 1D array-like
        x-coordinates.
    Y : 1D array-like
        y-coordinates.

    Returns
    -------
    f : lambda function
        f function.
    """

    return lambda x: Y[next((i for i in range(len(X)) if X[i] >= x), -1)]
