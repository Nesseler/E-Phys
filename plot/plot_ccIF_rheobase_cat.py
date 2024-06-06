# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 17:29:07 2024

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

IF_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_col = 'i_input')
IF_inst_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_col = 'i_input')

active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')
passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col = 'cell_ID')

# get parameters of first AP
fstAP_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_col = 'cell_ID')

# cc_rest
activity_df = pd.read_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col = 'cell_ID')

# get cell_ID
cell_IDs = passiv_properties_df.index.to_list()

# import meta data sheet
MetaData = pd.read_excel(table_file, sheet_name="MetaData", index_col='cell_ID')
MetaData = MetaData.loc[cell_IDs, :]

plt_df = pd.concat([active_properties_df, passiv_properties_df, fstAP_df, activity_df, MetaData], axis = 1)

# parameters
params_toplot = ['rheobase_rel', 'n_rheobasespikes', 'delta_vrest_to_vthres', 
                 't_peaks', 'v_threshold', 'v_AHP_amplitude', 'FWHM']

# recalculation
plt_df['t_peaks'] = plt_df['t_peaks'] - cc_IF_parameters['t_pre']
plt_df['t_threshold'] = plt_df['t_threshold'] - cc_IF_parameters['t_pre']
plt_df['delta_vrest_to_vthres'] = plt_df['v_threshold'] - plt_df['v_rest']

# %% figure

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

for idx, parameter in enumerate(params_toplot):
    
    ax = axs_cats[idx + 1]
    
    if parameter == 'v_AHP_amplitude':
        ax.plot([-0.4, 2.4],[0, 0], c = colors_dict['primecolor'], linestyle = '--', lw = 1)
    
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


    if parameter == 'rheobase_rel':
        ax.set_ylabel('Relative rheobase [pA]')
        ax.set_ylim([0-2, 210+2])
        ax.spines['left'].set_bounds([0, 210])
        ax.set_yticks(np.arange(0, 210+1, 50))
        ax.set_yticks(np.arange(0, 210+1, 10), minor = True)

    elif parameter == 'n_rheobasespikes':
        ax.set_ylabel('Number of rheobase spikes [#]')
        ax.set_ylim([0-0.1, 10+0.1])
        ax.spines['left'].set_bounds([0, 10])
        ax.set_yticks(np.arange(0,10+1, 5))
        ax.set_yticks(np.arange(0,10+1, 2), minor = True)

    elif parameter == 'delta_vrest_to_vthres':
        ax.set_ylabel('Delta resting membrane potential\nto rheobase spike threshold [mV]')
        ax.set_ylim([-5-1, 45+1])
        ax.spines['left'].set_bounds([-5, 45])
        ax.set_yticks(np.arange(0, 45+1, 10))
        ax.set_yticks(np.arange(-5, 45+1, 5), minor = True)

    elif parameter == 't_peaks':
        ax.set_ylabel('Time to rheobase spike [ms]')
        ax.set_ylim([0-10, 700+10])
        ax.spines['left'].set_bounds([0, 700])
        ax.set_yticks(np.arange(0,700+1, 100))
        ax.set_yticks(np.arange(0,700+1, 50), minor = True)

    elif parameter == 'v_threshold':
        ax.set_ylabel('Voltage at threshold\nof rheobase spike [mV]')
        ax.set_ylim([-80-2, -35+2])
        ax.spines['left'].set_bounds([-80, -35])
        ax.set_yticks(np.arange(-80, -35+1, 10))
        ax.set_yticks(np.arange(-80, -35+1, 5), minor = True)
        
    elif parameter == 'v_AHP_amplitude':
        ax.set_ylabel('Voltage delta to AHP [mV]')
        ax.set_ylim([-20-1, 20+1])
        ax.spines['left'].set_bounds([-20, 20])
        ax.set_yticks(np.arange(-20, 20+1, 10))
        ax.set_yticks(np.arange(-20, 20+1, 2.5), minor = True)  

    elif parameter == 'FWHM':
        ax.set_ylabel('Rheobase spike FWHM [ms]')
        ax.set_ylim([0.5-0.1, 2+0.1])
        ax.spines['left'].set_bounds([0.5, 2])
        ax.set_yticks(np.arange(0.5, 2+.1, 0.5))
        ax.set_yticks(np.arange(0.5, 2+.1, 0.25), minor = True)




[ax.grid(False) for ax in axs_cats]

# despine
[ax.spines[spine].set_visible(False) for ax in axs_cats for spine in ['top', 'right']]
[ax.spines['bottom'].set_bounds([0, 2]) for ax in axs_cats]

# x tick labels
[ax.set_xlabel('') for ax in axs_cats]
[ax.set_xticklabels(['BAOT/\nMeA', 'MeA', 'BAOT'], rotation = 45) for ax in axs_cats]
[ax.set_xticklabels(['', '', '']) for ax in axs_cats[:4]]

save_figures(fig_cats, 'ccIF-rheobase_properties', figure_dir, darkmode_bool)























