"""
Module plot provides useful functions for plots.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors

# DEFAULT VARIABLES

_markers = mpl.markers.MarkerStyle.filled_markers   # default markers list
_linestyles = ('-', '--', '-.', ':')                # default linestyles list

# FUNCTIONS AND CLASSES

def set_font_size(font_size):
    """
    Set matplotlib font size.

    Parameters
    ----------
    font_size : int
        Font size.
    """

    mpl.rcParams.update({'font.size': font_size})

def list_colormap(value_list, colormap='jet'):
    """
    Creates hash table of colors from colormap, defined according to value_list
    index, with value_list elements as keys.

    Parameters
    ----------
    value_list : list
        List of values.
    colormap : matplotlib colormap
        Colormap to use. (default: 'jet')

    Returns
    -------
    colors : hash table
        Hash table of colors.
    """

    cmap = plt.get_cmap(colormap)                               # colormap
    norm = colors.Normalize(vmin=0, vmax=len(value_list) - 1)   # normalise colormap according to list index
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)        # associates scalar to color

    return {value_list[index]: scalarMap.to_rgba(index)
        for index in range(len(value_list))}

def list_markers(value_list, marker_list=_markers):
    """
    Creates hash table of markers from markers_list, defined according to
    value_list index, with value_list elements as keys.

    Parameters
    ----------
    value_list : list
        List of values.
    marker_list : list of matplotlib markers
        List of markers to use. (default: active_particles.plot.plot._markers)

    Returns
    -------
    markers : hash table
        Hash table of markers.
    """

    return {value_list[index]: marker_list[index]
        for index in range(len(value_list))}
