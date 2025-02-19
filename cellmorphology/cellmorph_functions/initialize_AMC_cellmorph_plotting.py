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
                          'glomerular_dendrites' : '#038CAB',
                          'basal_dendrites'      : '#3C51AB',
                          'LOT_dendrites'        : '#5103AB', 
                          'undefined_dendrites'  : '#0357AB',
                      'axons'     : '#B10318'}

# set colors
darkmode_bool = False

if darkmode_bool:
    plt.style.use('default')
    plt.style.use('dark_background')
    
elif not darkmode_bool:
    plt.style.use('default')

# set font size
mtl.rcParams.update({'font.size': 9, 'font.family' : 'Arial'})
