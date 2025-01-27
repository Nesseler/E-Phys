# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 15:34:09 2025

@author: nesseler
"""


import matplotlib as mtl
import matplotlib.pyplot as plt
import seaborn as sbn

from functions.functions_plotting import save_figures, get_colors, get_figure_size, apply_axis_settings, plot_half_violin, get_blkr_colors, remove_y_spines_n_ticks

# set colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)
blkr_colors = get_blkr_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})



