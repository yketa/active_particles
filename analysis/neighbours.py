"""
Module neighbours defines the object NeighboursGrid, which is initiated by
calculating the neighbours grid and allows one to then access neighbouring
particles of any point.
"""

import numpy as np

class NeighboursGrid:
    """
    Considering a 2D square box, a grid is built by dividing the box in smaller
    boxes of length at least equal to r_cut. We call neighbours grid the hash
    table which keys are 2-uples of indexes of the created grid and which
    values are lists of the indexes of particles contained in the corresponding
    boxes.

    Neighbours grids are used to speed up calculations when looking for
    neighbouring particles within a cut-off radius. Once a neighbours grid is
    calculated by initiating this object, neighbouring particles within r_cut
    of any point can be found with self.get_neighbours(...).
    """

    def __init__(self, positions, box_size, r_cut):
        """
        Calculates neighbours grid from particle positions, box size and
        cut-off radius.

        Parameters
        ----------
        positions : (N, 2) shaped array
            Positions of the particles.
        box_size : float
            Length of the 2D square box.
        r_cut : float
            Cut-off radius.
        """

        self.cases_rcut = int(box_size/r_cut)   # number of boxes in one direction for the neighbours grid
        self.spacing = box_size/self.cases_rcut # spacing between two consecutive boxes
        self.neighbours_grid = {(x, y): [] for x in range(self.cases_rcut)\
            for y in range(self.cases_rcut)}    # neighbours grid

        for index in range(len(positions)):     # for all particles
            if (np.abs(positions[index]) > box_size).any(): continue    # do not consider particles outside the box
            self.neighbours_grid[tuple(
                (np.array(positions[index])//self.spacing + self.cases_rcut)
                %self.cases_rcut)] += [index]                           # add index to corresponding list

    def get_neighbours(self, point):
        """
        Returns the list of indexes of particles contained in boxes
        neighbouring the box containing a given point as well as the particles'
        indexes in this box.

        WARNING: Not all particles which indexes are returned by this function
        are within a distance r_cut of point, however all particles which are
        within r_cut will be returned.

        Parameters
        ----------
        point : array of length 2
            Position of the point of which we want neighbouring particles
            within r_cut.

        Returns
        -------
        neighbours : list
            List of the indexes of particles neighbouring point within r_cut.
        """

        index =\
            (np.array(point)//self.spacing + self.cases_rcut)%self.cases_rcut   # index of the box containing point
        sum_index = lambda ind, inc:\
            tuple((ind + inc + self.cases_rcut)%self.cases_rcut)                # function summing box indexes

        neighbours = []                 # neighbours list
        for inc_x in [-1, 0, 1]:        # increment in x index
            for inc_y in [-1, 0, 1]:    # increment in y index
                neighbours += self.neighbours_grid[
                    sum_index(index, np.array([inc_x, inc_y]))
                    ]                   # adds indexes of neighbouring box to list of neighbours

        return neighbours
