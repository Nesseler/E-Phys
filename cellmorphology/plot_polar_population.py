# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 09:58:17 2024

@author: nesseler
"""

import pandas as pd
from os import mkdir
from os.path import join, exists
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import math


from parameters.directories_win import cell_morph_descrip_dir, table_file, cell_morph_plots_dir

from functions.functions_plotting import save_figures, get_colors, get_figure_size

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# load histogram data
polar_plot_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_occurrances.xlsx'), index_col = 'orientation_rad')
polar_plot_dendrites_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_dendrites_occurrances.xlsx'), index_col = 'orientation_rad')
polar_plot_axons_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_axons_occurrances.xlsx'), index_col = 'orientation_rad')

orientation_labels = ['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv']

bins_angles = polar_plot_occurrances.index.to_list()

# get binsize
resul_binsize = np.diff(bins_angles)[-1]

# get colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

# define directory
cell_morph_plots_polar_pop_dir = join(cell_morph_plots_dir, 'polar_plots_population')

# %% all combined

# polar histogram 
fig_hist, ax_hist = plt.subplots(subplot_kw={'projection': 'polar'},
                                  layout = 'constrained',
                                  height_ratios= [1],
                                  width_ratios=[1],
                                  figsize = get_figure_size(width = 328.67/2))

# set title
fig_hist.suptitle('Terminal branch orientation\nall cells')

# get sum of all occurrances
polar_plot_occurrances_sum = polar_plot_occurrances.sum(axis = 1)

# plot histogram as barplot
ax_hist.bar(bins_angles, polar_plot_occurrances_sum,
            width = resul_binsize, 
            align = 'edge',
            edgecolor = 'none',
            color = 'gray')

# x axis
ax_hist.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
ax_hist.set_xticklabels(orientation_labels)

# yaxis
ymax = np.ceil(polar_plot_occurrances_sum.max()/10)*10
ax_hist.set_ylim([0, ymax])
ax_hist.set_yticks(np.arange(0, ymax + 1, 20))

# grid
ax_hist.grid(True, alpha = 0.5)

plt.show()

# save figures
save_figures(fig_hist, 'population_polar_plot-all_cells', cell_morph_plots_polar_pop_dir, darkmode_bool)

# %% separated regions

fig_hist_regions, ax_hist_regions = plt.subplots(nrows = 1, 
                                                 ncols = 2,
                                                 subplot_kw={'projection': 'polar'},
                                                 layout = 'constrained',
                                                 height_ratios= [1],
                                                 width_ratios=[1, 1],
                                                 figsize = get_figure_size(),
                                                 sharex = True,
                                                 sharey = True)

# set title
fig_hist_regions.suptitle('Terminal branch orientation')

regions = ['MeA', 'BAOT']

# ymax
ymax = 0

for i_region, region in enumerate(regions):
    
    # get cell IDs in region
    region_cell_IDs = MetaData[MetaData['Region'] == region].index.to_list()
    region_cell_IDs = [c for c in region_cell_IDs if c in polar_plot_occurrances.columns]
    
    # limit dataframe to cells within region
    region_polar_plot_occurrances = polar_plot_occurrances[region_cell_IDs]
    
    # get sum of all occurrances
    region_polar_plot_occurrances_sum = region_polar_plot_occurrances.sum(axis = 1)

    # define axis
    ax = ax_hist_regions[i_region]
    
    # set title
    ax.set_title(region)
    
    # plot histogram as bar plot
    ax.bar(bins_angles, region_polar_plot_occurrances_sum,
           width = resul_binsize, 
           align = 'edge',
           edgecolor = 'none',
           color = region_colors[region])
    
    ax.grid(True, alpha = 0.5)
    # ax.set_axisbelow(True)
    
    # get ymax
    if ymax < region_polar_plot_occurrances_sum.max():
        ymax = region_polar_plot_occurrances_sum.max()
    
# x axis
ax_hist_regions[-1].set_xticks(np.arange(0, np.pi*2, np.pi / 4))
ax_hist_regions[-1].set_xticklabels(orientation_labels)

# yaxis
ymax = np.ceil(ymax/10)*10
ax_hist_regions[-1].set_ylim([0, ymax])
ax_hist_regions[-1].set_yticks(np.arange(0, ymax + 1, 10))

plt.show()

# save figures
save_figures(fig_hist_regions, 'population_polar_plot-regions', cell_morph_plots_polar_pop_dir, darkmode_bool)


# %% dendrites & axons

# polar histogram 
fig_hist_cc_neurites, ax_hist_cc_neurites = plt.subplots(subplot_kw={'projection': 'polar'},
                                                         layout = 'constrained',
                                                         height_ratios= [1],
                                                         width_ratios=[1],
                                                         figsize = get_figure_size(width = 328.67/2))

# set title
fig_hist_cc_neurites.suptitle('Terminal branch orientation\nall cells')

# get sum of all occurrances
polar_plot_dendrites_occurrances_sum = polar_plot_dendrites_occurrances.sum(axis = 1)
polar_plot_axons_occurrances_sum = polar_plot_axons_occurrances.sum(axis = 1)

# plot histogram as barplot
ax_hist_cc_neurites.bar(bins_angles, polar_plot_dendrites_occurrances_sum,
                        width = resul_binsize, 
                        align = 'edge',
                        edgecolor = 'none',
                        color = 'gray')

ax_hist_cc_neurites.bar(bins_angles, polar_plot_axons_occurrances_sum,
                        width = resul_binsize,
                        bottom = polar_plot_dendrites_occurrances_sum,
                        align = 'edge',
                        edgecolor = 'none',
                        color = 'lightgray')

# x axis
ax_hist_cc_neurites.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
ax_hist_cc_neurites.set_xticklabels(orientation_labels)

# yaxis
ymax = np.ceil(polar_plot_occurrances_sum.max()/10)*10
ax_hist_cc_neurites.set_ylim([0, ymax])
ax_hist_cc_neurites.set_yticks(np.arange(0, ymax + 1, 20))

# grid
ax_hist_cc_neurites.grid(True, alpha = 0.5)

plt.show()

# save figures
save_figures(fig_hist_cc_neurites, 'population_polar_plot-all_cells-cc_neurites', cell_morph_plots_polar_pop_dir, darkmode_bool)


# %% dendrites & axons + regions

# polar histogram 
fig_hist_cc_neurites_regions, ax_hist_cc_neurites_regions = plt.subplots(nrows = 1, 
                                                                         ncols = 2,
                                                                         subplot_kw={'projection': 'polar'},
                                                                         layout = 'constrained',
                                                                         height_ratios= [1],
                                                                         width_ratios=[1, 1],
                                                                         figsize = get_figure_size(),
                                                                         sharex = True,
                                                                         sharey = True)

# set title
fig_hist_cc_neurites_regions.suptitle('Terminal branch orientation')


for i_region, region in enumerate(regions):
    
    # get cell IDs in region
    region_cell_IDs = MetaData[MetaData['Region'] == region].index.to_list()
    region_cell_IDs = [c for c in region_cell_IDs if c in polar_plot_occurrances.columns]
    
    # limit dataframe to cells within region
    region_polar_plot_dendrites_occurrances = polar_plot_dendrites_occurrances[region_cell_IDs]
    region_polar_plot_axons_occurrances = polar_plot_axons_occurrances[region_cell_IDs]
    
    # get sum of all occurrances
    region_polar_plot_dendrites_occurrances_sum = region_polar_plot_dendrites_occurrances.sum(axis = 1)
    region_polar_plot_axons_occurrances_sum = region_polar_plot_axons_occurrances.sum(axis = 1)

    # define axis
    ax = ax_hist_cc_neurites_regions[i_region]
    
    # set title
    ax.set_title(region)
    
    # plot histogram as bar plot
    ax.bar(bins_angles, region_polar_plot_dendrites_occurrances_sum,
           width = resul_binsize, 
           align = 'edge',
           edgecolor = 'none',
           color = region_colors[region])
    
    ax.bar(bins_angles, region_polar_plot_axons_occurrances_sum,
           bottom = region_polar_plot_dendrites_occurrances_sum,
           width = resul_binsize, 
           align = 'edge',
           edgecolor = 'none',
           color = 'lightgray')
    
    ax.grid(True, alpha = 0.5)
    # ax.set_axisbelow(True)


# x axis
ax_hist_cc_neurites_regions[-1].set_xticks(np.arange(0, np.pi*2, np.pi / 4))
ax_hist_cc_neurites_regions[-1].set_xticklabels(orientation_labels)

# yaxis
ymax = 120
ax_hist_cc_neurites_regions[-1].set_ylim([0, ymax])
ax_hist_cc_neurites_regions[-1].set_yticks(np.arange(0, ymax + 1, 10))

plt.show()

# save figures
save_figures(fig_hist_cc_neurites_regions, 'population_polar_plot-regions-cc_neurites', cell_morph_plots_polar_pop_dir, darkmode_bool)

