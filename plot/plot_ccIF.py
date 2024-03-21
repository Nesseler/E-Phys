# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:53:21 2024

@author: nesseler
"""

import pandas as pd
from os.path import join
import seaborn as sbn
import matplotlib.pyplot as plt
import numpy as np

from functions.functions_plotting import set_font_sizes, get_colors, get_figure_size, save_figures

from parameters.directories_win import quant_data_dir, cell_descrip_dir, table_file, figure_dir


# load data

IF_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_col = 'i_input')
IF_inst_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_col = 'i_input')

active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')
passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col = 'cell_ID')

# get cell IDs
cell_IDs = IF_df.columns.to_list()

# %% import meta data sheet

MetaData = pd.read_excel(table_file,
                      sheet_name="MetaData",
                      index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

combined_df = pd.concat([active_properties_df, passiv_properties_df, MetaData], axis = 1)



# %% calc mean and std

IF_mean = IF_df.mean(axis = 1, skipna = True)
IF_std = IF_df.std(axis = 1, skipna = True)


# %% plotting


cell_ID = 'E-120'





darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)

# initialize figure
fig_IF, axs_IF = plt.subplots(nrows = 1,
                              ncols = 1,
                              dpi = 600,
                              layout = 'constrained',
                              figsize = get_figure_size())


# x axis at zero
axs_IF.hlines(y = 0, xmin = -100, xmax = 400,
                      lw = 0.5,
                      color = colors_dict['primecolor'],
                      linestyle = '--')


for cell_ID in cell_IDs:

    axs_IF.plot(IF_df[cell_ID], 
                color = 'gray', # color = region_colors[MetaData.at[cell_ID, 'Region']],  # color = colors_dict['primecolor'],#
                linewidth = 1.5)

axs_IF.fill_between(x = IF_mean.index.to_list(), 
                    y1 = IF_mean - IF_std,
                    y2 = IF_mean + IF_std, 
                    color = 'gray',
                    alpha = 0.5,
                    edgecolor = None)

axs_IF.plot(IF_mean.index.to_list(), IF_mean,
            color = colors_dict['primecolor'],
            linewidth = 2)


# x axis
xmax = 400
xmin = -100
axs_IF.set_xlabel('Input current [pA]')
axs_IF.set_xlim([xmin - 10, xmax + 10])
axs_IF.set_xticks(ticks = np.arange(xmin, xmax + 1, 100))
axs_IF.set_xticks(ticks = np.arange(xmin, xmax + 1, 25), minor = True)
axs_IF.spines['bottom'].set_bounds([xmin, xmax])



# y axis
ymax = 80
ymin = 0
axs_IF.set_ylabel('Frequency [Hz]')
axs_IF.set_ylim([ymin - 2, ymax + 2])
axs_IF.set_yticks(ticks = np.arange(-ymin, ymax + 1, 10))
axs_IF.set_yticks(ticks = np.arange(-ymin, ymax + 1, 5), minor = True)
axs_IF.spines['left'].set_bounds([ymin, ymax])

[axs_IF.spines[spine].set_visible(False) for spine in ['top', 'right']]

axs_IF.grid(False)

plt.show()

save_figures(fig_IF, 'ccIF-IF-all_gray+mean', figure_dir, darkmode_bool)


# %%


fig_IF_regions, axs_IF_regions = plt.subplots(nrows = 2,
                                              ncols = 2,
                                              layout = 'constrained',
                                              dpi = 600,
                                              sharey = True,
                                              sharex = True,
                                              figsize = get_figure_size())
axs_IF_regions[0][1].remove()

axs_IF_regions = axs_IF_regions.flatten()

axs_IF_regions = [axs_IF_regions[0], axs_IF_regions[2], axs_IF_regions[3]]

for ax in axs_IF_regions:
    # x axis at zero
    ax.hlines(y = 0, xmin = xmin, xmax = xmax,
              lw = 0.5,
              color = colors_dict['primecolor'],
              linestyle = '--')

regions = ['BAOT/MeA', 'MeA', 'BAOT']

for region_idx, region in enumerate(regions):
    
    for cell_ID in cell_IDs:
        
        cell_region = MetaData.at[cell_ID, 'Region']
        
        if cell_region == region:
            
            axs_IF_regions[region_idx].plot(IF_df[cell_ID], 
                                            color = region_colors[cell_region],
                                            linewidth = 1.5)



# x axis
xmax = 400
xmin = -100
axs_IF_regions[1].set_xlabel('Input current [pA]')
axs_IF_regions[2].set_xlabel('Input current [pA]')
axs_IF_regions[2].set_xlim([xmin - 5, xmax + 5])
axs_IF_regions[2].set_xticks(ticks = np.arange(xmin, xmax + 1, 100))
axs_IF_regions[2].set_xticks(ticks = np.arange(xmin, xmax + 1, 25), minor = True)

# y axis
ymax = 80
ymin = 0
axs_IF_regions[0].set_ylabel('Frequency [Hz]')
axs_IF_regions[1].set_ylabel('Frequency [Hz]')
axs_IF_regions[2].set_ylim([ymin - 2, ymax + 2])
axs_IF_regions[2].set_yticks(ticks = np.arange(-ymin, ymax + 1, 10))
axs_IF_regions[2].set_yticks(ticks = np.arange(-ymin, ymax + 1, 5), minor = True)


for ax in axs_IF_regions:
    # limit spines
    ax.spines['bottom'].set_bounds([xmin, xmax])
    ax.spines['left'].set_bounds([ymin, ymax])
    
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

    ax.grid(False)

plt.show()



save_figures(fig_IF_regions, 'ccIF-IF-sep_regions', figure_dir, darkmode_bool)


# %% sep regions + means

fig_IF_regions, axs_IF_regions = plt.subplots(nrows = 2,
                                              ncols = 2,
                                              layout = 'constrained',
                                              dpi = 600,
                                              sharey = True,
                                              sharex = True,
                                              figsize = get_figure_size())
axs_IF_regions[0][1].remove()

axs_IF_regions = axs_IF_regions.flatten()

axs_IF_regions = [axs_IF_regions[0], axs_IF_regions[2], axs_IF_regions[3]]


for ax in axs_IF_regions:
    # x axis at zero
    ax.hlines(y = 0, xmin = xmin, xmax = xmax,
              lw = 0.5,
              color = colors_dict['primecolor'],
              linestyle = '--')


regions = ['BAOT/MeA', 'MeA', 'BAOT']

for region_idx, region in enumerate(regions):
    
    cell_IDs_inregion = MetaData[MetaData['Region'] == region].index.to_list()
    
    IF_region_df = IF_df[cell_IDs_inregion]
    
    ax = axs_IF_regions[region_idx]
    
    for cell_ID in cell_IDs_inregion:
        ax.plot(IF_region_df[cell_ID], 
                color = region_colors[region],
                linewidth = 1.5)

    IF_region_mean = IF_region_df.mean(axis = 1, skipna = True)
    IF_region_std = IF_region_df.std(axis = 1, skipna = True)

    ax.fill_between(x = IF_region_mean.index.to_list(), 
                    y1 = IF_region_mean - IF_region_std,
                    y2 = IF_region_mean + IF_region_std, 
                    color = region_colors[region],
                    alpha = 0.5,
                    edgecolor = None)

    ax.plot(IF_region_mean.index.to_list(), IF_region_mean,
            color = colors_dict['primecolor'],
            linewidth = 2)



# x axis
xmax = 400
xmin = -100
axs_IF_regions[1].set_xlabel('Input current [pA]')
axs_IF_regions[2].set_xlabel('Input current [pA]')
axs_IF_regions[2].set_xlim([xmin - 5, xmax + 5])
axs_IF_regions[2].set_xticks(ticks = np.arange(xmin, xmax + 1, 100))
axs_IF_regions[2].set_xticks(ticks = np.arange(xmin, xmax + 1, 25), minor = True)

# y axis
ymax = 80
ymin = 0
axs_IF_regions[0].set_ylabel('Frequency [Hz]')
axs_IF_regions[1].set_ylabel('Frequency [Hz]')
axs_IF_regions[2].set_ylim([ymin - 2, ymax + 2])
axs_IF_regions[2].set_yticks(ticks = np.arange(-ymin, ymax + 1, 10))
axs_IF_regions[2].set_yticks(ticks = np.arange(-ymin, ymax + 1, 5), minor = True)


for ax in axs_IF_regions:
    # limit spines
    ax.spines['bottom'].set_bounds([xmin, xmax])
    ax.spines['left'].set_bounds([ymin, ymax])

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

    ax.grid(False)

plt.show()


save_figures(fig_IF_regions, 'ccIF-IF-sep_regions+means', figure_dir, darkmode_bool)








