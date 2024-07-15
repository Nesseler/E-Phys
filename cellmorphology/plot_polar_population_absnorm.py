# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 19:06:08 2024

@author: nesseler
"""


import pandas as pd
from os import mkdir
from os.path import join, exists
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import math


from parameters.directories_win import cell_morph_descrip_dir, table_file, cell_morph_plots_dir, cell_morph_plots_dir

from functions.functions_plotting import save_figures, get_colors, get_figure_size, set_font_sizes

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
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# define directory
cell_morph_plots_polar_pop_dir = join(cell_morph_plots_dir, 'polar_plots_population')


# %% min max normalised version

### normalize polar occurrences
# all
polar_plot_occurrances_normed = polar_plot_occurrances.div(polar_plot_occurrances.sum(), axis = 1)

# dendrites
polar_plot_dendrites_occurrances_normed = polar_plot_dendrites_occurrances.div(polar_plot_dendrites_occurrances.sum(), axis = 1)

# axons
polar_plot_axons_occurrances_normed = polar_plot_axons_occurrances.div(polar_plot_axons_occurrances.sum(), axis = 1)


### create average per bin
# all
polar_plot_occurrances_normed_mean = polar_plot_occurrances_normed.mean(axis = 1)

# dendrites
polar_plot_dendrites_occurrances_normed_mean = polar_plot_dendrites_occurrances_normed.mean(axis = 1)

# axons
polar_plot_axons_occurrances_normed_mean = polar_plot_axons_occurrances_normed.mean(axis = 1)


### combine to one dataframe
ALL_pp_normed_mean_occs = pd.DataFrame({'all' : polar_plot_occurrances_normed_mean.to_list(),
                                        'dendrites' : polar_plot_dendrites_occurrances_normed_mean.to_list(), 
                                        'axons' : polar_plot_axons_occurrances_normed_mean.to_list()}, 
                                       index = orientation_labels)


### create normalized polar occurrences per region


def region_polar_plot_occurrences(pp_all, pp_dendrites, pp_axons, region, neurite_types, orientation_labels):
    
# pp_all = polar_plot_occurrances
# pp_dendrites = polar_plot_dendrites_occurrances
# pp_axons = polar_plot_axons_occurrances
# region = 'BAOT'
# neurite_types = ['all', 'dendrites', 'axons']
# orientation_labels = orientation_labels

    # define output dataframe
    pp_normed_mean_occs = pd.DataFrame(index = orientation_labels, columns = neurite_types)
    
    # get cell IDs in region
    region_cell_IDs = MetaData[MetaData['Region'] == region].index.to_list()
    region_cell_IDs = [c for c in region_cell_IDs if c in pp_all.columns]
    
    # loop through neurite types
    for occurrence_df, neurite_type in zip([pp_all, pp_dendrites, pp_axons], neurite_types):
        
        # limit dataframe to cells within region
        region_occ = occurrence_df[region_cell_IDs]
        
        # get normed occurrences
        region_occ_normed = region_occ.div(region_occ.sum(), axis = 1)
        
        # create average per bin
        region_occ_normed_mean = region_occ_normed.mean(axis = 1).to_list()
    
        # add to dataframe
        pp_normed_mean_occs[neurite_type] = region_occ_normed_mean
    
    # return dataframe
    return pp_normed_mean_occs


# get normed mean occurrences for MeA
MeA_pp_normed_mean_occs = region_polar_plot_occurrences(pp_all = polar_plot_occurrances,
                                                        pp_dendrites = polar_plot_dendrites_occurrances,
                                                        pp_axons = polar_plot_axons_occurrances,
                                                        region = 'MeA',
                                                        neurite_types = ['all', 'dendrites', 'axons'],
                                                        orientation_labels = orientation_labels)

# get normed mean occurrences for BAOT
BAOT_pp_normed_mean_occs = region_polar_plot_occurrences(pp_all = polar_plot_occurrances,
                                                          pp_dendrites = polar_plot_dendrites_occurrances,
                                                          pp_axons = polar_plot_axons_occurrances,
                                                          region = 'BAOT',
                                                          neurite_types = ['all', 'dendrites', 'axons'],
                                                          orientation_labels = orientation_labels)



# normed figure

# polar histogram 
fig_norm, axs_norm = plt.subplots(nrows = 3,
                                  ncols = 3,
                                  subplot_kw={'projection': 'polar'},
                                  layout = 'constrained',
                                  height_ratios= [1, 1, 1],
                                  width_ratios=[1, 1, 1],
                                  figsize = get_figure_size(width = 150, height = 150))

# flatten axis
axs_norm = axs_norm.flatten()

# # set title
# fig_hist.suptitle('Terminal branch orientation\nall cells')

## both regions

# plot histogram as barplot
## neurites (all)
axs_norm[0].bar(bins_angles, polar_plot_occurrances_normed_mean,
                width = resul_binsize, 
                align = 'edge',
                edgecolor = 'none',
                color = 'gray')

## dendrites
axs_norm[1].bar(bins_angles, polar_plot_dendrites_occurrances_normed_mean,
                width = resul_binsize, 
                align = 'edge',
                edgecolor = 'none',
                color = 'gray')

## axons
axs_norm[2].bar(bins_angles, polar_plot_axons_occurrances_normed_mean,
                width = resul_binsize, 
                align = 'edge',
                edgecolor = 'none',
                color = 'gray')



### regions
for plot_i, neurite_type in enumerate(['all', 'dendrites', 'axons']):
        
    ## MeA
    # plot histogram as barplot
    axs_norm[plot_i + 3].bar(bins_angles, MeA_pp_normed_mean_occs[neurite_type],
                             width = resul_binsize, 
                             align = 'edge',
                             edgecolor = 'none',
                             color = region_colors['MeA'])

    ## BAOT
    # plot histogram as barplot
    axs_norm[plot_i + 6].bar(bins_angles, BAOT_pp_normed_mean_occs[neurite_type],
                    width = resul_binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = region_colors['BAOT'])



# column header
for ax, neurite_type in zip(axs_norm[:3], ['Neurites', 'Dendrites', 'Axons']):
    ax.annotate(neurite_type, xy = (0.5, 1), xytext = (0, 30), 
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline', fontsize = 9)
    
# row header
for ax, plot_type in zip(axs_norm[::3], ['Combined', 'MeA', 'BAOT']):
    ax.annotate(plot_type, xy=(0, 0.5), xytext=(-20, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                ha='center', va='center', fontsize = 9, rotation = 90)


### y axis

# neurites
for ax in axs_norm[::3]:
    ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '10 %', '', '20 %', ''])
    ax.set_rlabel_position(80)

# dendrites
for ax in axs_norm[1::3]:
    ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '10 %', '', '20 %', ''])
    ax.set_rlabel_position(80)

# axons
for ax in axs_norm[2::3]:
    ax.set_yticks(ticks = np.arange(0, 0.35 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '', '30 %', ''], va = 'top')
    ax.set_rlabel_position(-80)


### x axis
for ax in axs_norm:
    # x axis
    ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
    ax.set_xticklabels(orientation_labels)


    # grid
    ax.grid(True, alpha = 0.5)
    ax.set_axisbelow(True)
    


fig_norm.align_labels()



set_font_sizes(9, 9)


save_figures(fig_norm, figure_name = 'population_polar_plots-normed_per_type-regions', save_dir = cell_morph_plots_dir,
             figure_format= 'both')


for occu_df, region_label in zip([ALL_pp_normed_mean_occs, MeA_pp_normed_mean_occs, BAOT_pp_normed_mean_occs], ['All', 'MeA', 'BAOT']):

    occu_df.to_excel(join(cell_morph_plots_dir, 
                          'polar_plots_population', 
                          'population_polar_plots-normed_per_type-regions-' + region_label + '.xlsx'), 
                     index = orientation_labels)
