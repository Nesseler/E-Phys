# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 18:22:40 2024

@author: nesseler
"""

import pandas as pd
from os.path import join
import seaborn as sbn
import matplotlib.pyplot as plt
import numpy as np

from functions.functions_plotting import set_font_sizes, get_colors, get_figure_size, save_figures

from parameters.directories_win import quant_data_dir, cell_descrip_dir, table_file, figure_dir
from parameters.PGFs import cc_IF_parameters

# load data

# IF_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_col = 'i_input')
# IF_inst_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_col = 'i_input')

active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')
passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col = 'cell_ID')

# # get parameters of first AP
# fstAP_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_col = 'cell_ID')

# # cc_rest
# activity_df = pd.read_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col = 'cell_ID')

# get cell_ID
cell_IDs = active_properties_df.index.to_list()

# import meta data sheet
MetaData = pd.read_excel(table_file, sheet_name="MetaData", index_col='cell_ID')
MetaData = MetaData.loc[cell_IDs, :]

plt_df = pd.concat([active_properties_df, passiv_properties_df, MetaData], axis = 1)

# parameters
params_toplot = ['max_freq', 'max_inst_freq', 'max_inst_initial_freq']


# %% figure

darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)

# initialize figure
fig_cats, axs_cats = plt.subplots(nrows = 1,
                                  ncols = 3,
                                  figsize = get_figure_size(width=246.502, height = 165.5/2),
                                  layout = 'constrained'
                                  )

# set font sizes
set_font_sizes()

# flatten numpy array of axis for easier handling
axs_cats = axs_cats.flatten()

for idx, parameter in enumerate(params_toplot):
    
    ax = axs_cats[idx]
    
    violin = sbn.violinplot(data = plt_df,
                            x = 'Region',
                            y = parameter,
                            bw = 0.3,
                            inner = 'quart',
                            linewidth = 1,
                            ax = ax,
                            order = ['BAOT/MeA', 'MeA', 'BAOT'])

    for l in violin.lines:
        l.set_color(colors_dict['primecolor'])

    for violin in violin.collections:
        violin.set_edgecolor(colors_dict['primecolor'])
        violin.set_facecolor('None')

    swarm = sbn.swarmplot(data = plt_df,
                          x = 'Region',
                          y = parameter, 
                          ax = ax,
                          hue = 'Region', 
                          palette = region_colors,
                          size = 5,
                          color = colors_dict['primecolor'],
                          order = ['BAOT/MeA', 'MeA', 'BAOT'])

    ax.legend().set_visible(False)


    if parameter == 'max_freq':
        ax.set_ylabel('Max frequency (number of spikes) [Hz]')
        ax.set_ylim([0-2, 180+2])
        ax.spines['left'].set_bounds([0, 180])
        ax.set_yticks(np.arange(0, 180+1, 50))
        ax.set_yticks(np.arange(0, 180+1, 10), minor = True)

    elif parameter == 'max_inst_freq':
        ax.set_ylabel('Max instantaneous\nspiking frequency [Hz]')
        ax.set_ylim([0-2, 180+2])
        ax.spines['left'].set_bounds([0, 180])
        ax.set_yticks(np.arange(0, 180+1, 50))
        ax.set_yticks(np.arange(0, 180+1, 10), minor = True)

    elif parameter == 'max_inst_initial_freq':
        ax.set_ylabel('Max initial instantaneous\nspiking frequency [Hz]')
        ax.set_ylim([0-2, 180+2])
        ax.spines['left'].set_bounds([0, 180])
        ax.set_yticks(np.arange(0, 180+1, 50))
        ax.set_yticks(np.arange(0, 180+1, 10), minor = True)



[ax.grid(False) for ax in axs_cats]

# despine
[ax.spines[spine].set_visible(False) for ax in axs_cats for spine in ['top', 'right']]
[ax.spines['bottom'].set_bounds([0, 2]) for ax in axs_cats]

# x tick labels
[ax.set_xlabel('') for ax in axs_cats]
[ax.set_xticklabels(['BAOT/\nMeA', 'MeA', 'BAOT'], rotation = 45) for ax in axs_cats]
# [ax.set_xticklabels(['', '', '']) for ax in axs_cats[:4]]

save_figures(fig_cats, 'ccIF-maxfreq_properties', figure_dir, darkmode_bool)