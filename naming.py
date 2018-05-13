"""
Module naming provides functions and objects to get and read names of all
directories and files used throughout this project, as well as a glossary of
commonly used variables.
"""

from active_particles.exponents import float_to_letters, letters_to_float,\
    significant_figures
from active_particles.init import get_env

from collections import OrderedDict
from itertools import chain
from copy import deepcopy

# GLOSSARY

var_glossary = {    # definitions of variables commonly used
    'density': 'surface packing fraction',
    'vzero': 'self-propelling velocity',
    'dr': 'rotational diffusion constant',
    'N': 'number of particles',
    'init_frame': 'initial frame',
    'dt': 'lag time',
    'int_max': 'maximum number of intervals',
    'Ncases': 'number of boxes in one direction of a grid',
    'r_cut': 'cut-off radius',
    'sigma': 'length scale of Gaussian function'
}

# DEFAULT NAMES

sim_directory = get_env('HOME') + '/active_particles_data'  # simulation data directory

parameters_file = 'param.p'                     # simulation parameters file
log_file = 'log-output.log'                     # simulation log output file
wrapped_trajectory_file = 'trajectory.gsd'      # wrapped trajectory file (with periodic boundary conditions)
unwrapped_trajectory_file = 'trajectory.dat'    # unwrapped trajectory file (without periodic boundary conditions)

# FILES NAMING

_image_extension = '.eps'  # default image extension

class _File:
    """
    Naming files.

    This is the superclass which provides methods to get and read names for
    all subclasses defined specifically per files which we will call
    'standards'.

    All subclasses must be initiated with
        > self.name : string
            Generic name.
        > self.parameters : ordered dictionary
            Hash table of parameters and their abbreviations.
        > self.extension : string
            File extension.
    """

    def filename(self, **definitions):
        """
        Returns a list of name parts which the hash table parameters enables to
        define.

        Optional keyword arguments
        --------------------------
        definitions : float, int or bool
            Defined parameters.

        Returns
        -------
        name_parts : list of strings
            List of defined name parts.
        """

        name_parts = []                 # list of defined name parts
        buffer = self.name              # consecutive defined attributes
        for param in self.parameters:   # looping through file name parameters

            if param in definitions:                        # if parameter is defined
                buffer += str(
                    self.parameters[param] +                # parameter abbreviation
                    float_to_letters(definitions[param]))   # defined value of parameter

            else:
                if buffer != '': name_parts += [buffer] # add buffer to defined name parts
                buffer = ''

        buffer += self.extension   # adding extension to name parts
        name_parts += [buffer]

        return name_parts

    def get_data(self, file, *parameter):
        """
        Returns values of parameters which can be read from the file name.

        Parameters
        ----------
        file : string
            File name.

        Optional positional arguments
        -----------------------------
        parameter : string
            Name of parameter in self.parameters.

        Returns
        -------
        par_values : list of floats
            List of requested parameters values.
        """

        par_values = [] # list of requested parameters values.

        for param in parameter:             # for all requested parameters
            if param in self.parameters:    # requested parameter exists
                par_values += [letters_to_float(
                    file.split(self.parameters[param])[1][
                    :significant_figures + 1]
                )]                          # parameter value translated to float

        return par_values

    def add_ext(self, ext_parameters, ext_extension):
        """
        From a default standard, this function returns an extended standard
        with additional parameters and different extension.

        Parameters
        ----------
        ext_parameters : ordered dictionary
            Hash table of additional parameters and their abbreviations.
        ext_extension : string
            Different file extension.

        Returns
        -------
        ext_self : active_particles.naming standard
            Extended standard object.
        """

        ext_self = deepcopy(self)                               # creates a deep copy of the standard
        ext_self.parameters = OrderedDict(chain(
            self.parameters.items(), ext_parameters.items()))   # extended ordered dictionary of parameters
        ext_self.extension = ext_extension                      # extended standard extension

        return ext_self

    def image(self):
        """
        This function is the default image name generator, which only changes
        file extension with __image_extension.
        """

        return self.add_ext(OrderedDict(), _image_extension)

class Css(_File):
    """
    Naming shear strain maps and shear strain correlation files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__()  # initialise with superclass

        self.name = 'Css' + endpoint()  # generic name
        self.parameters = OrderedDict([
            ('density', '_D'), ('vzero', '_V'), ('dr', '_R'), ('N', '_N'),
            ('init_frame', '_I'), ('dt', '_T'), ('int_max', '_M'),
            ('Ncases', '_C'), ('r_cut', '_RCUT'), ('sigma', '_SIGM')
        ])                              # parameters and corresponding abbreviations (in order)
        self.extension = '.pickle'      # file extension

class Ccc(Css):
    """
    Naming displacement vorticity maps and displacement vorticity correlation
    files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__()  # initialise with superclass

        self.name = 'Ccc' + endpoint()

def endpoint():
    """
    Consider a variable A which depends on space, time and a lag time (e.g.,
    displacement). It is possible to measure
    (i) A between times t and t + dt for a particle which is at position r at
    time t: A(r, t ; t, t + dt)
    (ii) A between times t and t + dt for a particle which is at position r at
    time t + dt: A(r, t + dt ; t, t + dt)
    and calculate correlations over space and time in both cases.

    We will add a 'b' following the name of the correlation for files which
    have produced considering case (i).
        e.g., for displacement correlations, filenames will begin with 'Cuub'
        if correlations were calculated for displacements between times t and
        t + dt between particles at position r at time t.

    Case (i) corresponds to environment variable 'ENDPOINT' set as False or not
    set, while case (ii) corresponds to environment variable 'ENDPOINT' set as
    True.
    """

    if get_env('ENDPOINT', default=False, vartype=bool): return ''
    return 'b'
