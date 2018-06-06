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
        Default value to return if environment variable does not exist or does
        not evaluate. (default: None)
    vartype : data type
        Desired data type. (default: string)
    """

    try:
        if vartype == bool and envvar[var_name] == 'False': # special custom case for which we want to return False
                return False
        try:
            return vartype(envvar[var_name])
        except ValueError:
            try:
                return vartype(eval(envvar[var_name]))
            except ValueError: raise KeyError
    except KeyError:
        return default

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
