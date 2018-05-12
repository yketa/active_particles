"""
Module init provides functions useful when initialising simulations or
analysis.
"""

from os import environ as environment

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
        return vartype(environment[var_name])
    except ValueError:
        try:
            return vartype(eval(environment[var_name]))
        except:
            return default
    except:
        return default
