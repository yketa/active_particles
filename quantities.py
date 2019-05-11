"""
Module quantities provides functions to calculate relevant physical quantities.
"""

def nD0_active(N, vzero, dr, L):
    """
    Returns product of particle density n = N/L^2 and active diffusion constant
    D_0 = vzero^2/2*dr.

    Parameters
    ----------
    N : int or float
        Number of particles.
    vzero : float
        Self-propelling velocity.
    dr : float
        Rotation diffusion constant.
    L : float
        Characteristic system length.

    Returns
    -------
    product : float
        n D_0
    """

    return (N*(vzero**2))/(2*dr*(L**2))

def nD0_thermal(N, kT, gamma, L):
    """
    Returns product of particle density n = N/L^2 and diffusion constant
    D_0 = 2*kT*N/gamma*L^2.

    Parameters
    ----------
    N : int or float
        Number of particles.
    kT : float
        Dimensionless temperature.
    gamma : float
        Brownian dumping coefficient.
    L : float
        Characteristic system length.

    Returns
    -------
    product : float
        n D_0
    """

    return (2*kT*N)/(gamma*(L**2))
