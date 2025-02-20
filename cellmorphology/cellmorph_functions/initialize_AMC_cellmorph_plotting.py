# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:27:56 2025

@author: nesseler
"""

# import packages for plotting
import matplotlib as mtl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sbn

# import custom functions
from functions.functions_plotting import save_figures, get_figure_size, apply_axis_settings, plot_half_violin, remove_spines_n_ticks, change_projection


# colors for neurites, dendrites, axons 
neurite_color_dict = {'soma' : '#D3D3D3',
                      'neurites'  : '#FFAD0A', 
                      'dendrites' : '#0424AC',
                          'glomerular_dendrites'    : '#0000FF',
                          'lateral_dendrites'       : '#55A0FB',
                          'LOTxing_dendrites'       : '#000080',
                          'nonglomerular_dendrites' : '#04A9CF',
                      'axons'     : '#B10318'}

#000080FF, #0000C0FF, #0000FFFF, #5757F9FF, #55A0FBFF, #90BFF9FF, #C8DEF9FF, #0000FFFF, #5757F9FF

# neurite_color_dict = {'soma' : '#D3D3D3',
#                       'neurites'  : '#FFAD0A', 
#                       'dendrites' : '#0424AC',
#                           'glomerular_dendrites'    : '#04A9CF',
#                           'lateral_dendrites'       : '#88A6E6',
#                           'LOTxing_dendrites'       : '#5103AB',
#                           'nonglomerular_dendrites' : '#03488D',
#                       'axons'     : '#B10318'}

# set colors
darkmode_bool = False

if darkmode_bool:
    plt.style.use('default')
    plt.style.use('dark_background')
    
elif not darkmode_bool:
    plt.style.use('default')

# set font size
mtl.rcParams.update({'font.size': 9, 'font.family' : 'Arial'})
