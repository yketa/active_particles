"""
Module plot provides useful functions for plots.
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors

def list_colormap(value_list, colormap='jet'):
    """
    Creates hash table of colors from colormap, defined according to value_list
    index, with vaulue_list elements as keys.

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
