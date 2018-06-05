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

from os import getcwd
from os import environ as envvar
from os import listdir as ls
from os.path import join as joinpath

# DEFAULT NAMES

sim_directory = joinpath(get_env('HOME'), 'active_particles_data')  # simulation data directory
out_directory = joinpath(sim_directory, 'out')                      # launch output directory

parameters_file = 'param.p'                     # simulation parameters file
log_file = 'log-output.log'                     # simulation log output file
wrapped_trajectory_file = 'trajectory.gsd'      # wrapped trajectory file (with periodic boundary conditions)
unwrapped_trajectory_file = 'trajectory.dat'    # unwrapped trajectory file (without periodic boundary conditions)

# GLOSSARY

class Glossary:
    """
    This class references VarInfo objects.
    """

    def __init__(self, *entries):
        """
        Positional arguments
        --------------------
        entries : active_particles.naming.VarInfo object
            Information about a variable.
        """

        self.entries = entries

    def __getitem__(self, varname):
        """
        Called by self[varname].

        Parameters
        ----------
        varname : string
            Generic name of variable.

        Returns
        -------
        entry : active_particles.naming.VarInfo
            Informations about variable 'varname' in glossary.
        """

        return next((entry for entry in self.entries
            if entry.name == varname), None)

    def __add__(self, entry):
        """
        Called by self + ...

        Parameters
        ----------
        entry : active_particles.naming.VarInfo object
            New entry to add to glossary.

        Returns
        -------
        gloassry : active_particles.naming.Glossary object
            Augmented glossary.
        """

        self.entries += (entry,)
        return self

class VarInfo:
    """
    This class save informations about a given variable.
    """

    def __init__(self, name, definition, symbol, format):
        """
        Parameters
        ----------
        name : string
            Generic name.
        definition : string
            Brief description of the variable.
        symbol : string
            LaTeX symbol to use in legends.
        format : formatting string
            Format to use in legends.
        """

        self.name = name
        self.definition = definition
        self.symbol = symbol
        self.format = format

    def __str__(self):
        """
        Prints description of variable.
        """

        return self.definition

    def legend(self, value):
        """
        Parameters
        ----------
        value : *
            Value of variable.

        Returns
        -------
        legend : string
            Formatted string '{variable} = {value}'.
        """

        return self.symbol + '=' + self.format.format(value)

glossary = Glossary(*map(lambda entry: VarInfo(*entry), (
        # | NAME | DEFINITION | SYMBOL | FORMAT |
        ('density', 'surface packing fraction', r'$\phi', '{:1.2f}'),
        ('vzero', 'self-propelling velocity', r'$\tilde{v}$', '{:.2e}'),
        ('dr', 'rotational diffusion constant', r'$\tilde{\nu}_r$', '{:.2e}'),
        ('N', 'number of particles', r'$N$', '{:.2e}'),
        ('init_frame', 'initial frame', r'$S_{init}$', '{:.2e}'),
        ('frame', 'time frame', r'$T$', '{:.2e}'),
        ('dt', 'lag time', r'$\Delta t$', '{:.2e}'),
        ('int_max', 'maximum number of intervals', r'$S_{max}$', '{:.2e}'),
        ('int_period', 'period of intervals', r'$S_{period}$', '{:.2e}'),
        ('Ncases', 'number of boxes in one direction of a grid',
            r'$N_{cases}$', '{:.2e}'),
        ('r_cut', 'cut-off radius', r'$r_{cut}$', '{:.2e}'),
        ('sigma', 'length scale of Gaussian function', r'$\sigma$', '{:.2e}'),
        ('box_size', 'length of the box in one dimension', r'$L$', '{:.3e}')
    )))

# FILES NAMING

_image_extension = '.eps'   # default image extension
_movie_extension = '.mp4'   # default movie extension

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
        Returns a list of name parts which the keyword arguments enables to
        define.

        Optional keyword arguments
        --------------------------
        definitions : float, int or bool
            Defined parameters (e.g., density=0.1).

        Returns
        -------
        name_parts : list of strings
            List of defined name parts.
        """

        name_parts = []                 # list of defined name parts
        buffer = self.name              # consecutive defined attributes
        for param in self.parameters:   # looping through file name parameters

            try:
                if definitions[param] == None: raise KeyError
                buffer += str(
                    self.parameters[param] +                # parameter abbreviation
                    float_to_letters(definitions[param]))   # defined value of parameter

            except KeyError:                            # if parameter is not defined
                if buffer != '': name_parts += [buffer] # add buffer to defined name parts
                buffer = ''

        buffer += self.extension   # adding extension to name parts
        name_parts += [buffer]

        return name_parts

    def filename_length(self):
        """
        Returns the length of the filename according to its name, parameters
        and extensions.

        Returns
        -------
        length : int
            Length of filename.
        """

        length = len(self.name)
        for parameter in self.parameters.values():
            length += len(parameter) + 1 + significant_figures
        length += len(self.extension)

        return length

    def get_files(self, directory=getcwd(), **definitions):
        """
        Returns list of files which correpond to the keyword arguments.

        Parameters
        ----------
        directory : string
            Files directory. (default: current working directory)

        Optional keyword arguments
        --------------------------
        definitions : float, int or bool
            Defined parameters (e.g., density=0.1).

        Returns
        -------
        files : list of strings
            List of files correponding to parameters.
        """

        length = self.filename_length()             # length of filename
        name_parts = self.filename(**definitions)   # list of defined name parts

        return [file for file in ls(directory) if
            all(part in file for part in name_parts)
            and len(file) == length]

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
        file extension with _image_extension.
        """

        return self.add_ext(OrderedDict(), _image_extension)

    def out(self):
        """
        This function is the default output file name generator, which only
        changes file extension with '.out'.
        """

        return self.add_ext(OrderedDict(), '.out')

class _CorFile(_File):
    """
    Naming correlation files.
    """

    def __init__(self, name, **kwargs):
        """
        Architecture of file name.

        Parameters
        ----------
        name : string
            Generic name of correlation file.

        Optional keyword arguments
        --------------------------
        ext_parameters : ordered dictionary
            Hash table of additional parameters and their abbreviations.
        """

        self.name = name + endpoint()   # generic name
        self.parameters = OrderedDict([
            ('density', '_D'), ('vzero', '_V'), ('dr', '_R'), ('N', '_N'),
            ('init_frame', '_I'), ('dt', '_T'), ('int_max', '_M'),
            ('Ncases', '_C')
        ])                              # parameters and corresponding abbreviations (in order)
        self.extension = '.pickle'      # file extension

        if 'ext_parameters' in kwargs:
            self.parameters = self.add_ext(kwargs['ext_parameters'],
                self.extension).parameters  # add additional parameters

        if 'BOX_SIZE' in envvar:                        # modified box size
            self.parameters = OrderedDict(chain(self.parameters.items(),
                {'box_size': '_B'}.items()))
        if 'X_ZERO' in envvar or 'Y_ZERO' in envvar:    # modified centre of the box
            self.parameters = OrderedDict(chain(self.parameters.items(),
                OrderedDict([('x_zero', '_X'), ('y_zero', '_Y')]).items()))

class Css(_CorFile):
    """
    Naming shear strain maps and shear strain correlation files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Css', ext_parameters=OrderedDict([('r_cut', '_RCUT'),
            ('sigma', '_SIGM')]))   # initialise with superclass

class Ccc(_CorFile):
    """
    Naming displacement vorticity maps and displacement vorticity correlation
    files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Ccc', ext_parameters=OrderedDict([('r_cut', '_RCUT'),
            ('sigma', '_SIGM')]))   # initialise with superclass

class Cuu(_CorFile):
    """
    Naming displacmeent correlation files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Cuu')  # initialise with superclass

class Cnn(_CorFile):
    """
    Naming density correlation files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Cnn')  # initialise with superclass

class Cww(_CorFile):
    """
    Naming displacement relative to centre of mass displacement correlation
    files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Cww')  # initialise with superclass

class Cdd(_CorFile):
    """
    Naming displacement norm correlation files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Cdd')  # initialise with superclass

class Cee(_CorFile):
    """
    Naming displacement direction correlation files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Cee')  # initialise with superclass

class Ctt(_CorFile):
    """
    Naming mean square norm of cross product of wave vector and displacement
    Fourier transform.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Ctt')  # initialise with superclass

class Cll(_CorFile):
    """
    Naming mean square norm of dot product of wave vector and displacement
    Fourier transform.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('Cll')  # initialise with superclass

class _FrameFile(_File):
    """
    Naming system image files.
    """

    def __init__(self, name, **kwargs):
        """
        Architecture of file name.

        Parameters
        ----------
        name : string
            Generic name of image file.

        Optional keyword arguments
        --------------------------
        ext_parameters : ordered dictionary
            Hash table of additional parameters and their abbreviations.
        """

        self.name = name        # generic name
        self.parameters = OrderedDict([
            ('density', '_D'), ('vzero', '_V'), ('dr', '_R'), ('N', '_N'),
            ('init_frame', '_I')
        ])                      # parameters and corresponding abbreviations (in order)
        self.extension = '.eps' # file extension

        if 'ext_parameters' in kwargs:
            self.parameters = self.add_ext(kwargs['ext_parameters'],
                self.extension).parameters  # add additional parameters

        if 'BOX_SIZE' in envvar:                        # modified box size
            self.parameters = OrderedDict(chain(self.parameters.items(),
                {'box_size': '_B'}.items()))
        if 'X_ZERO' in envvar or 'Y_ZERO' in envvar:    # modified centre of the box
            self.parameters = OrderedDict(chain(self.parameters.items(),
                OrderedDict([('x_zero', '_X'), ('y_zero', '_Y')]).items()))

    def movie(self, folder=False):
        """
        Movie name generator.

        Parameters
        ----------
        folder : bool
            Generate folder name.
        """

        if folder:
            ext_extension = ''
        else: ext_extension = _movie_extension

        ext_parameters = OrderedDict()
        if 'FINAL_FRAME' in envvar:
            ext_parameters = OrderedDict(chain(ext_parameters.items(),
                OrderedDict([('frame_fin', '_F')]).items()))
        if 'FRAME_PERIOD' in envvar:
            ext_parameters = OrderedDict(chain(ext_parameters.items(),
                OrderedDict([('frame_per', '_P')]).items()))
        if 'FRAME_MAXIMUM' in envvar:
            ext_parameters = OrderedDict(chain(ext_parameters.items(),
                OrderedDict([('frame_max', '_M')]).items()))
        return self.add_ext(ext_parameters, ext_extension)

class Velocity(_FrameFile):
    """
    Naming system images with arrows along velocity directions files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('v')   # initialise with superclass

class Trajectory(_FrameFile):
    """
    Naming system images with displacements shown as arrows files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('t', ext_parameters=OrderedDict([('dt', '_T')]))   # initialise with superclass

class Displacement(_FrameFile):
    """
    Naming system images with arrows along displacement directions files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        super().__init__('u', ext_parameters=OrderedDict([('dt', '_T')]))   # initialise with superclass

class Msd(_File):
    """
    Naming mean square displacement files.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        self.name = 'msd_sterr' # generic name
        self.parameters = OrderedDict([
            ('density', '_D'), ('vzero', '_V'), ('dr', '_R'), ('N', '_N'),
            ('init_frame', '_I'), ('int_max', '_M'), ('int_period', '_P')
        ])                      # parameters and corresponding abbreviations (in order)
        self.extension = '.csv' # file extension

class VarN(_File):
    """
    Naming local densities densities.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        self.name = 'varN'          # generic name
        self.parameters = OrderedDict([
            ('density', '_D'), ('vzero', '_V'), ('dr', '_R'), ('N', '_N'),
            ('init_frame', '_I'), ('int_max', '_M'), ('Ncases', '_C'),
            ('box_size', '_B')
        ])                          # parameters and corresponding abbreviations (in order)
        self.extension = '.pickle'  # file extension

class AHB2D(_File):
    """
    Naming simulations of 2D active brownian uniformly polydisperse particles
    directories.
    """

    def __init__(self):
        """
        Architecture of file name.
        """

        self.name = ''      # generic name
        self.parameters = OrderedDict([
            ('density', 'D'), ('vzero', '_V'), ('dr', '_R'), ('N', '_N'),
            ('launch', '_L')
        ])                  # parameters and corresponding abbreviations (in order)
        self.extension = '' # file extension

    def filename(self, directory=getcwd(), **definitions):
        """
        Returns a list of name parts which the keyword arguments enables to
        define.

        Parameters
        ----------
        directory : string
            Directory in which the folder will be created.

        Optional keyword arguments
        --------------------------
        definitions : float, int or bool
            Defined parameters (e.g., density=0.1).

        Returns
        -------
        name_parts : list of strings
            List of defined name parts.
        """

        if not('launch' in definitions):    # launch number not defined
            definitions['launch'] = len(super().get_files(**definitions,
                launch=None))
        return super().filename(**definitions)


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
