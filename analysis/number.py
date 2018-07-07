"""
Module number provides functions to count particles in a system.
"""

import numpy as np

from active_particles.maths import count

def count_particles(w_traj, *frames, box_size=None, centre=(0, 0)):
    """
    Returns the number of particles at frames frames contained in the square
    box of length box_size centred on centre.

    Parameters
    ----------
    w_traj : active_particles.dat.Gsd
		Wrapped trajectory object.
    frames : int
        Frame indexes.
    box_size : float
        Box length. (default: None)
        NOTE: if box_size is None, then box_size is considered to be the
        system box size.
    centre : 2-uple of float
        Box centre. (default: (0, 0))

    Returns
    -------
    N : list of int
        Number of particles in each frames.
    """

    return list(map(
        lambda frame: count(w_traj.position(frame, centre=centre),
            box_size if box_size != None else w_traj.box_size(time=frame)),
        frames))
