# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 10:00:51 2024

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

# get parameters of first AP
fstAP_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_col = 'cell_ID')
todrop = ['v_peaks', 't_peaks', 't_threshold', 'idx_threshold', 'v_AHP', 't_AHP', 'idx_AHP', 'v_AHP_amplitude', 't_to_AHP', 'v_HM', 't1_HM', 't2_HM', 'SR_ms']
fstAP_df.drop(columns = todrop, inplace = True)

# concatenate active and passive properties
# IF_cat_df = pd.concat([passiv_properties_df, active_properties_df], axis = 1)
IF_cat_df = passiv_properties_df
# IF_cat_df.drop(columns = 'rheobase_step_idx', inplace = True)

# get cell IDs
cell_IDs = IF_cat_df.index.to_list()

# %% import meta data sheet

MetaData = pd.read_excel(table_file,
                      sheet_name="MetaData",
                      index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

plt_df = pd.concat([IF_cat_df, MetaData], axis = 1)

# %%


darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)

# initialize figure
fig_cats, axs_cats = plt.subplots(nrows = 1,
                                  ncols = 3,
                                  figsize = get_figure_size(width = 247.752, height = 73.25),
                                  layout = 'constrained')


# flatten numpy array of axis for easier handling
axs_cats = axs_cats.flatten()


for idx, parameter in enumerate(IF_cat_df.columns):
    
    ax = axs_cats[idx]
    
    violin = sbn.violinplot(data = plt_df,
                            x = 'Region',
                            y = parameter,
                            bw = 0.3,
                            inner = 'quart',
                            linewidth = 1,
                            ax = ax,
                            size = 0.9,
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
                          size = 3,
                          color = colors_dict['primecolor'],
                          order = ['BAOT/MeA', 'MeA', 'BAOT'])

    ax.legend().set_visible(False)


    if parameter == 'rheobase_abs':
        ax.set_ylabel('Absolute rheobase [pA]')
        ax.set_ylim([-52, 202])
        ax.spines['left'].set_bounds([-50, 200])
        ax.set_yticks(np.arange(-50, 200+1, 50))
        ax.set_yticks(np.arange(-50, 200+1, 10), minor = True)
        
        
    elif parameter == 'rheobase_rel':
        ax.set_ylabel('Relative rheobase [pA]')
        ax.set_ylim([-52, 202])
        ax.spines['left'].set_bounds([-50, 200])
        ax.set_yticks(np.arange(-50, 200+1, 50))
        ax.set_yticks(np.arange(-50, 200+1, 10), minor = True)
     
    elif parameter == 'v_thres_rheobase_spike':
        ax.set_ylabel('Voltage at threshold\nof rheobase spike [mV]')
        ax.set_ylim([-87, -28])
        ax.spines['left'].set_bounds([-85, -30])
        ax.set_yticks(np.arange(-80, -30+1, 10))
        ax.set_yticks(np.arange(-85, -30+1, 5), minor = True)
        
    elif parameter == 'r_input':
        ax.set_ylabel(r'Input resistance [M$\Omega$]')
        ax.set_ylim([-20, 1820])
        ax.spines['left'].set_bounds([0, 1800])
        ax.set_yticks(np.arange(0, 1750+1, 500))
        ax.set_yticks(np.arange(0, 1750+1, 100), minor = True)  
        
    elif parameter == 'tau_mem':
        ax.set_ylabel('Membrane time constant\n[ms]')
        ax.set_ylim([-2, 42])
        ax.spines['left'].set_bounds([0, 40])
        ax.set_yticks(np.arange(0, 40+1, 10))
        ax.set_yticks(np.arange(0, 40+1, 5), minor = True)
        
    elif parameter == 'c_mem':
        ax.set_ylabel('Membrane capacitance\n[pF]')
        ax.set_ylim([-2, 102])
        ax.spines['left'].set_bounds([0, 100])
        ax.set_yticks(np.arange(0, 100+1, 20))
        ax.set_yticks(np.arange(0, 100+1, 5), minor = True)

    
[ax.grid(False) for ax in axs_cats]

# despine
[ax.spines[spine].set_visible(False) for ax in axs_cats for spine in ['top', 'right']]
[ax.spines['bottom'].set_bounds([0, 2]) for ax in axs_cats]

# x tick labels
[ax.set_xlabel('') for ax in axs_cats]
[ax.set_xticklabels(['BAOT/\nMeA', 'MeA', 'BAOT'], rotation = 45) for ax in axs_cats]
# [ax.set_xticklabels(['', '', '']) for ax in axs_cats[:3]]

# set font sizes
set_font_sizes(12)

temp_fig_dir = 'C:/Users/nesseler/Desktop/TAC-presentation_data/ePhys'

save_figures(fig_cats, 'ccIF-active_passive_properties', temp_fig_dir, darkmode_bool, figure_format='both')



# %%



# sbn.jointplot(data = plt_df[plt_df['Region'] == 'BAOT'], 
#               x = 'r_input', 
#               y = 'v_thres_rheobase_spike',
#               hue = 'Region', 
#               palette = region_colors)

# # plt.xlim([0, 60])
# plt.xlim([0, 1750])





