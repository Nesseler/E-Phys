# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 18:12:32 2024

@author: nesseler
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir, figure_dir
from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size

activity_df = pd.read_excel(os.path.join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col = 'cell_ID')


activity_df = activity_df.sort_values('n_spikes', ascending = True)

cell_IDs = list(activity_df.index)
n_cells = len(cell_IDs)

# %% EVENTPLOT + V_REST FIGURE + N_spike

darkmode_bool = True

# get colors for plotting
colors_dict = get_colors(darkmode_bool)

fig_v_rest, axs_v_rest = plt.subplots(1,2,
                                      gridspec_kw={'width_ratios': [3,1]},
                                      figsize = get_figure_size(),
                                      layout = 'tight')

set_font_sizes()

tick_size = 0.9

for cell_idx, cell_ID in enumerate(cell_IDs):
    
    # read time points of spikes as string representation of a list
    t_spikes = activity_df.at[cell_ID, 't_spikes'].strip('][').split(', ')
    
    # convert individual string elements to floats
    # check for empty list not to be converted to float
    if len(t_spikes) > 1:
        t_spikes = [float(t) for t in t_spikes]
    else:
        t_spikes = []
    
    axs_v_rest[0].eventplot(t_spikes,
                            orientation = 'horizontal', 
                            lineoffsets=cell_idx, 
                            linewidth = 1.5,
                            linelengths=0.9, 
                            color = colors_dict['color2'])


axs_v_rest[0].set_ylim([0-(tick_size/2), n_cells-1+(tick_size/2)])
axs_v_rest[0].set_yticks(ticks = np.arange(5 - 1, n_cells+1, 5), 
                         labels = np.arange(5, n_cells + 1, 5))
axs_v_rest[0].set_yticks(ticks = np.arange(0, n_cells, 1), 
                         minor = True)
axs_v_rest[0].set_ylabel('Cells [#]')

axs_v_rest[0].set_xlim([0, 30])
axs_v_rest[0].set_xlabel('Time [s]')
axs_v_rest[0].set_xticks(np.arange(0, 30+1, 10))
axs_v_rest[0].set_xticks(np.arange(0, 30+1, 1), minor = True)

# creating specific color pallet for seaborn plotting functions
color_pal = {'silent' : colors_dict['color1'], 'spiking' : colors_dict['color2']}

violins = sbn.violinplot(data = activity_df, 
                         x = "activity", 
                         y = "v_rest", 
                         inner = 'quart', 
                         ax = axs_v_rest[1], 
                         palette = color_pal,
                         linewidth = 1.5)

swarms = sbn.swarmplot(data = activity_df, 
                       x = "activity", 
                       y = "v_rest", 
                       size = 7, 
                       ax = axs_v_rest[1], 
                       color = colors_dict['primecolor'])

for l in violins.lines:
    l.set_color(colors_dict['primecolor'])

[violin.set_edgecolor(colors_dict['primecolor']) for violin in violins.collections]

axs_v_rest[1].set_ylabel('Resting membrane potential [mV]')
axs_v_rest[1].set_ylim([-100, -40])
axs_v_rest[1].set_yticks(np.arange(-100, -40+1, 10))
axs_v_rest[1].set_yticks(np.arange(-100, -40+1, 5), minor = True)

axs_v_rest[1].set_xlabel('')

[ax.grid(False) for ax in axs_v_rest]

save_figures(fig_v_rest, 'Resting_n_eventplot_nspikes', figure_dir, darkmode_bool)