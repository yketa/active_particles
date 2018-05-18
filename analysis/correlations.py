"""
Module correlations provides functions to calculate 2D fields correlations.
"""

import numpy as np

def corField2D_scalar(field):
    """
    2D correlation field of a scalar field. Correlations are calculated with
    use of Fast Fourier Transform.

    Parameters
    ----------
    field : 2D array like
        Scalar field to extract correlations from.
        Points are supposed to be uniformly distributed.

    Returns
    -------
    C : 2D numpy array
        Unnormalised correlation field.
        C[0, 0] is the origin, points are uniformly distributed.
    Norm : float
        Norm of correlation field.
    """

    FFT = np.fft.fft2(field)                       # FFT of scalar field
    C = np.real(np.fft.ifft2(np.conj(FFT)*FFT))    # Unnormalised correlation field
    Norm = np.sum(field**2)                        # Norm of correlation field

    return C, Norm

def corField2D_scalar_average(field_list):
    """
    2D correlation field, averaged from a list of scalar fields. Correlations
    are calculated with use of Fast Fourier Transform.

    Parameters
    ----------
    field_list : list of 2D array like
        List of scalar fields to extract correlations from.
        Points are supposed to be uniformly distributed.

    Returns
    -------
    C : 2D numpy array
        Normalised averaged correlation field.
        C[0, 0] is the origin, points are uniformly distributed.
    """

    C = (lambda c, Norm: c/Norm)(*tuple(
        np.sum(list(map(corField2D_scalar, field_list)), axis=0)
        ))  # normalised averaged correlation field
    return C

def corField2D_vector(field):
    """
    2D correlation field of a vector field. Correlations are calculated with
    use of Fast Fourier Transform.

    Parameters
    ----------
    field : (n, n, 2) shaped array like
        Vector field to extract correlations from.
        Points are supposed to be uniformly distributed.

    Returns
    -------
    C : 2D numpy array
        Unnormalised correlation field.
        C[0, 0] is the origin, points are uniformly distributed.
    xCL : float
        Unnormalised longitudinal correlation of field projected on the first
        direction of space at distance equal to field grid spacing.
    yCL : float
        Unnormalised longitudinal correlation of field projected on the second
        direction of space at distance equal to field grid spacing.
    xCT : float
        Unnormalised transversal correlation of field projected on the first
        direction of space at distance equal to field grid spacing.
    yCT : float
        Unnormalised transversal correlation of field projected on the second
        direction of space at distance equal to field grid spacing.
    Norm : float
        Norm of correlation field.
    """

    xfield = field[:, :, 0]                 # projection of field on the first direction of space
    xC, xNorm = corField2D_scalar(xfield)   # unnormalised correlation field and its norm associated to field projection on the first direction of space

    yfield = field[:, :, 1]                 # projection of field on the second direction of space
    yC, yNorm = corField2D_scalar(yfield)   # unnormalised correlation field and its norm associated to field projection on the second direction of space

    C = xC + yC                    # correlation field of field
    xCL, yCL = xC[0, 1], yC[1, 0]  # longitudinal correlations in first and second directions of space
    xCT, yCT = xC[1, 0], yC[0, 1]  # transversal correlations in first and second directions of space
    Norm = xNorm + yNorm           # norm of correlation field

    return C, xCL, yCL, xCT, yCT, Norm

def corField2D_vector_average(field_list):
    """
    2D correlation field, averaged from a list of vector fields. Correlations
    are calculated with use of Fast Fourier Transform.

    Parameters
    ----------
    field_list : list of (n, n, 2) shaped array like
        List of vector fields to extract correlations from.
        Points are supposed to be uniformly distributed.

    Returns
    -------
    C : 2D numpy array
        Normalised averaged correlation field.
        C[0, 0] is the origin, points are uniformly distributed.
    CL : float
        Normalised averaged longitudinal correlation of field at distance equal
        to field grid spacing.
    CT : float
        Normalised averaged transversal correlation of field at distance equal
        to field grid spacing.
    """

    C, CL, CT = (lambda c, xCL, yCL, xCT, yCT, Norm: (
        c/Norm,
        (xCL + yCL)/(2*Norm),
        (xCT + yCT)/(2*Norm)
        ))(*tuple(np.sum(list(map(corField2D_vector, field_list)), axis=0)))    # normalised averaged correlation field, longitudinal and transversal correlations

    return C, CL, CT

def corField2D_vector_average_Cnn(field_list, Cnn):
    """
    2D correlation field, averaged from a list of vector fields. Correlations
    are calculated with use of Fast Fourier Transform.

    Compared to correlations.corField2D_vector_average, this function also
    divides values of longitudial and transversal correlations, in each
    direction of space, by the value of the corresponding normalised
    density correlation.
    WARNING: This correction with the density correlation is not applied to the
    correlation field.

    Parameters
    ----------
    field_list : list of (n, n, 2) shaped array like
        List of vector fields to extract correlations from.
        Points are supposed to be uniformly distributed.
    Cnn : (n, n) shaped array like
        Normalised density correlation.

    Returns
    -------
    C : 2D numpy array
        Normalised averaged correlation field.
        C[0, 0] is the origin, points are uniformly distributed.
    CL : float
        Normalised averaged longitudinal correlation of field at distance equal
        to field grid spacing, corrected with density correlation.
    CT : float
        Normalised averaged transversal correlation of field at distance equal
        to field grid spacing, corrected with density correlation.
    """

    C, CL, CT = (lambda c, xCL, yCL, xCT, yCT, Norm: (
        c/Norm,
        (xCL/Cnn[0, 1] + yCL/Cnn[1, 0])/(2*Norm),
        (xCT/Cnn[1, 0] + yCT/Cnn[0, 1])/(2*Norm)
        ))(*tuple(np.sum(list(map(corField2D_vector, field_list)), axis=0)))    # normalised averaged correlation field, longitudinal and transversal correlations
