#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:14:20 2024

@author: moritznesseler
"""

import pandas as pd
from os.path import join
import seaborn as sbn
import matplotlib.pyplot as plt
import numpy as np

from functions.functions_plotting import set_font_sizes, get_colors, get_figure_size, save_figures

from parameters.directories_win import quant_data_dir, cell_descrip_dir, table_file, figure_dir


# load data
sagdelta_df = pd.read_excel(join(cell_descrip_dir, 'cc_sag-sagdelta.xlsx'), index_col = 'cell_ID')

sagdelta_df = sagdelta_df.query('sag_delta.notnull()')

cell_IDs = sagdelta_df.index.to_list()

# %% import meta data sheet

MetaData = pd.read_excel(table_file,
                       sheet_name="MetaData",
                       index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

plt_df = pd.concat([sagdelta_df, MetaData], axis = 1)

# %%


darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)

# initialize figure
fig_cats, axs_cats = plt.subplots(nrows = 2,
                                  ncols = 4,
                                  figsize = get_figure_size(),
                                  layout = 'constrained')

# set font sizes
set_font_sizes()

# flatten numpy array of axis for easier handling
axs_cats = axs_cats.flatten()
fig_cats.delaxes(axs_cats[0])
fig_cats.delaxes(axs_cats[2])

parameters_to_plot = ['sag_delta', 'n_reboundspikes', 'reboundspike_t_peak', 'reboundspike_v_threshold', 'reboundspike_v_amplitude', 'reboundspike_FWHM']


for idx, parameter in enumerate(parameters_to_plot):
    
    ax = axs_cats[idx + 1]
    
    if idx > 0:
        ax = axs_cats[idx + 2]
    
    violin = sbn.violinplot(data = plt_df,
                            x = 'Region',
                            y = parameter,
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


#     if parameter == 'rheobase_abs':
#         ax.set_ylabel('Absolute rheobase [pA]')
#         ax.set_ylim([-52, 202])
#         ax.spines['left'].set_bounds([-50, 200])
#         ax.set_yticks(np.arange(-50, 200+1, 50))
#         ax.set_yticks(np.arange(-50, 200+1, 10), minor = True)
        
        
#     elif parameter == 'rheobase_rel':
#         ax.set_ylabel('Relative rheobase [pA]')
#         ax.set_ylim([-52, 202])
#         ax.spines['left'].set_bounds([-50, 200])
#         ax.set_yticks(np.arange(-50, 200+1, 50))
#         ax.set_yticks(np.arange(-50, 200+1, 10), minor = True)
     
#     elif parameter == 'v_thres_rheobase_spike':
#         ax.set_ylabel('Voltage at threshold\nof rheobase spike [mV]')
#         ax.set_ylim([-87, -28])
#         ax.spines['left'].set_bounds([-85, -30])
#         ax.set_yticks(np.arange(-80, -30+1, 10))
#         ax.set_yticks(np.arange(-85, -30+1, 5), minor = True)
        
#     elif parameter == 'r_input':
#         ax.set_ylabel(r'Input resistance [M$\Omega$]')
#         ax.set_ylim([-20, 1820])
#         ax.spines['left'].set_bounds([0, 1800])
#         ax.set_yticks(np.arange(0, 1750+1, 500))
#         ax.set_yticks(np.arange(0, 1750+1, 100), minor = True)  
        
#     elif parameter == 'tau_mem':
#         ax.set_ylabel('Membrane time constant\n[ms]')
#         ax.set_ylim([-2, 42])
#         ax.spines['left'].set_bounds([0, 40])
#         ax.set_yticks(np.arange(0, 40+1, 10))
#         ax.set_yticks(np.arange(0, 40+1, 5), minor = True)
        
#     elif parameter == 'c_mem':
#         ax.set_ylabel('Membrane capacitance\n[pF]')
#         ax.set_ylim([-2, 102])
#         ax.spines['left'].set_bounds([0, 100])
#         ax.set_yticks(np.arange(0, 100+1, 20))
#         ax.set_yticks(np.arange(0, 100+1, 5), minor = True)

    
[ax.grid(False) for ax in axs_cats]

# despine
[ax.spines[spine].set_visible(False) for ax in axs_cats for spine in ['top', 'right']]
[ax.spines['bottom'].set_bounds([0, 2]) for ax in axs_cats]

# x tick labels
[ax.set_xlabel('') for ax in axs_cats]
[ax.set_xticklabels(['BAOT/\nMeA', 'MeA', 'BAOT'], rotation = 45) for ax in axs_cats[4:]]
axs_cats[1].set_xticklabels(['', '', ''])
axs_cats[3].set_xticklabels(['', '', ''])

# save_figures(fig_cats, 'ccIF-active_passive_properties', figure_dir, darkmode_bool)



# # %%



# sbn.jointplot(data = plt_df[plt_df['Region'] == 'BAOT'], 
#               x = 'r_input', 
#               y = 'v_thres_rheobase_spike',
#               hue = 'Region', 
#               palette = region_colors)

# # plt.xlim([0, 60])
# plt.xlim([0, 1750])