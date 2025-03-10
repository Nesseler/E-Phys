# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 11:33:49 2025

@author: nesseler
"""


# import packages for plotting
import matplotlib as mtl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import seaborn as sbn

# import custom functions
from functions.functions_plotting import save_figures, get_figure_size, apply_axis_settings, plot_half_violin, remove_spines_n_ticks, change_projection

# set colors for regions
region_colors = {'BAOT'     : '#7a66fc',
                 'MeA'      : '#ff8d00',
                 'BAOT/MeA' : 'gray'}

# colors for neurites, dendrites, axons 
neurite_color_dict = {'all' : {'neurites' : '#FFAD0AFF', 'dendrites' : '#0424ACFF', 'axons' : '#B10318FF', 'soma' : '#D3D3D3'},
                      'MeA' : {'neurites' : '#FDC067FF', 'dendrites' : '#6EC5ABFF', 'axons' : '#751C6DFF', 'soma' : '#D3D3D3'},
                      'BAOT': {'neurites' : '#FF8811FF', 'dendrites' : '#046E8FFF', 'axons' : '#D44D5CFF', 'soma' : '#D3D3D3'}}

# colors for spinyness categories
spines_color_dict = {'both' : {'high' : '#776B5D' , 'moderate' : '#B0A695', 'low': '#EBE3D5'},
                     'MeA'  : {'high' : '#A02334' , 'moderate' : '#EB5B00', 'low': '#FFB200'},
                     'BAOT' : {'high' : '#201E43' , 'moderate' : '#134B70', 'low': '#508C9B'}}

# set colors
darkmode_bool = False

if darkmode_bool:
    plt.style.use('default')
    plt.style.use('dark_background')
    primecolor = 'w'
    
elif not darkmode_bool:
    plt.style.use('default')
    primecolor = 'k'

# set font size
mtl.rcParams.update({'font.size': 9, 'font.family' : 'Arial'})
