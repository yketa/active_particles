"""
Module init provides functions useful when initialising simulations or
analysis.
"""

from os.path import join as joinpath
from os.path import exists as pathexists
from os import makedirs
from shutil import rmtree as rmr
from os import environ as envvar
import sys
import atexit

def to_vartype(input, default=None, vartype=str):
    """
    Returns input converted to vartype or default if the conversion fails.

    Parameters
    ----------
    input : *
        Input value (of any type).
    default : *
        Default value to return if conversion fails.
    vartype : data type
        Desired data type. (default: string)

    Returns
    -------
    output : vartype or type(default)
        Input converted to vartype data type or default.
    """

    if vartype == bool and input == 'False': return False   # special custom case

    try:
        try:
            return vartype(input)
        except ValueError: return vartype(eval(input))
    except: return default

def get_env(var_name, default=None, vartype=str):
    """
    Returns environment variable with desired data type.

    WARNING: get_env function uses eval function to evaluate environment
    variable strings if necessary, therefore extra cautious is recommended when
    using it.

    Parameters
    ----------
    var_name : string
        Name of environment variable.
    default : *
        Default value to return if environment variable does not exist or if
        conversion fails. (default: None)
    vartype : data type
        Desired data type. (default: string)

    Returns
    -------
    var : vartype or type(default)
        Environment variable converted to vartype data type of default.
    """

    return to_vartype(envvar[var_name], default=default, vartype=vartype)

def get_env_list(var_name, delimiter=';', default=None, vartype=str):
    """
    Returns list from environment variable containing values delimited with
    delimiter to be converted to vartype data type or taken to be default if
    the conversion fails.
    NOTE: Returns empty list if the environment variable does not exist or is
    an empty string.

    Parameters
    ----------
    var_name : string
        Name of environment variable.
    delimiter : string
        Pattern which delimits values to be evaluated in environment variable.
    default : *
        Default value to return if individual value in environment variable
        does not exist or if conversion fails. (default: None)
    vartype : data type
        Desired data type. (default: string)

    Returns
    -------
    var_list : list of vartype of type(default)
        List of individual environment variables values converted to vartype
        data type or default.
    """

    if not(var_name in envvar) or envvar[var_name] == '': return []
    return list(map(
        lambda var: to_vartype(var, default=default, vartype=vartype),
        envvar[var_name].split(delimiter)
        ))

class StdOut:
    """
    Enables to set output stream to file and revert this setting.
    """

    def __init__(self):
        """
        Saves original standard output as attribute.
        """

        self.stdout = sys.stdout    # original standard output

    def set(self, output_file):
        """
        Sets output to file.

        Parameters
        ----------
        output_file : file object
            Output file.
        """

        try:
            self.output_file.close()    # if output file already set, close it
        except AttributeError: pass

        self.output_file = output_file  # output file
        sys.stdout = self.output_file   # new output stream

        atexit.register(self.revert)    # close file when exiting script

    def revert(self):
        """
        Revers to original standard output.
        """

        try:
            self.output_file.close()
            sys.stdout = self.stdout    # revert to original standart output
        except AttributeError: pass     # no custom output was set

def mkdir(directory, replace=False):
    """
    Creates directory if not existing, erases and recreates it if replace is
    set to be True.

    Parameters
    ----------
    directory : string
        Name of directory.
    """

    if pathexists(directory):
        if not(replace): return
        rmr(directory)
    makedirs(directory)

def slurm_output(output_dir, naming_standard, attributes):
    """
    Sets standard output to file when launching from Slurm job scheduler.
    Writes job ID to output file.

    Parameters
    ----------
    output_dir : string
        Output file directory.
    naming_standard : active_particles.naming standard
        Naming standard to name output file.
    attributes : hash table
        Attributes which define ENTIRELY output file name.
    """

    mkdir(output_dir)   # create output directory if not existing

    output_filename, = naming_standard.out().filename(**attributes) # output file name
    output_file = open(joinpath(output_dir, output_filename), 'w')  # output file
    output_file.write('Job ID: %i\n\n'
        % get_env('SLURM_JOB_ID', vartype=int))					    # write job ID to output file

    stdout = StdOut()
    stdout.set(output_file)	# set output file as standard output
