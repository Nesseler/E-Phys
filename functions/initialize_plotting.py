# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 15:34:09 2025

@author: nesseler
"""


import matplotlib as mtl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

import seaborn as sbn

from functions.functions_plotting import save_figures, get_colors, get_figure_size, apply_axis_settings, plot_half_violin, get_blkr_colors, remove_spines_n_ticks, remove_spines_ticks_labels, change_projection

# load axis settings
from parameters.axis_settings import axis_dict

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)
blkr_colors = get_blkr_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9, 
                     'font.family' : 'Arial', 
                     'axes.titlesize' : 9})



