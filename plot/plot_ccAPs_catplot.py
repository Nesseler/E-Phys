# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 16:59:43 2024

@author: nesseler
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mtl
import seaborn as sbn

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir, figure_dir

# custom functions
from functions.functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes


# %% import 

# FWHM, tpeaks, vamplitude
measure_str = 'mean_' + 'FWHM'

mean_df = pd.read_excel(os.path.join(cell_descrip_dir, f'ccAPs-{measure_str}.xlsx'), index_col='frequencies')

cell_IDs = list(mean_df.columns)

frequencies = list(mean_df.index)
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]


# %% transpose and melt Dataframe to be used with seaborn functions

mean_df = mean_df.transpose()
mean_df_melt = pd.melt(mean_df, var_name='frequency', value_name=measure_str, ignore_index= False)     


# %% color code

norm_min = 1
norm_max = 75
cmap_str = 'plasma'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
cmap.set_array(freqs_int)


# creating specific color pallet for seaborn plotting functions
color_pal = {'1Hz'  : cmap.to_rgba(freqs_int[0]), 
             '5Hz'  : cmap.to_rgba(freqs_int[1]),
             '10Hz' : cmap.to_rgba(freqs_int[2]),
             '30Hz' : cmap.to_rgba(freqs_int[3]),
             '50Hz' : cmap.to_rgba(freqs_int[4]),
             '75Hz' : cmap.to_rgba(freqs_int[5])}


# %% plotting

darkmode_bool = True

colors_dict, _ = get_colors(darkmode_bool)

fig_APcats, axs_APcats = plt.subplots(nrows = 1, 
                                      ncols = 1,
                                      layout = 'constrained',
                                      dpi = 600,
                                      figsize = get_figure_size(),
                                      sharex = 'col')

set_font_sizes()


fig_APcats.colorbar(cmap, 
                  ax = axs_APcats,
                  ticks = freqs_int,
                  label = 'Stimulation frequency [Hz]')


violins1 = sbn.violinplot(data = mean_df_melt, 
                          x = 'frequency', 
                          y = measure_str, 
                          inner = 'quart', 
                          ax = axs_APcats, 
                          linewidth = 1,
                          palette = color_pal, 
                          width = 0.9)

swarms = sbn.swarmplot(data = mean_df_melt, 
                        x = 'frequency', 
                        y = measure_str, 
                        ax = axs_APcats,
                        size = 6,
                        color = colors_dict['primecolor'])


# get positions of all points in swarmplots to plot the connecting lines
positions_df = pd.DataFrame()

for idx, frequency in enumerate(frequencies):
    positions = np.array(axs_APcats.collections[6 + idx].get_offsets())
    
    cur_df = pd.DataFrame({'x' : positions[:, 0],
                           'y' : positions[:, 1],
                           'frequency' : [freqs_int[idx]] * len(positions),
                           'freq_str' : [frequency] * len(positions),
                           'cell_ID' : cell_IDs
                            })

    positions_df = pd.concat([positions_df, cur_df])
    
     

for idx, cell_ID in enumerate(cell_IDs):
    lines = sbn.lineplot(data = positions_df[positions_df.cell_ID == cell_ID], 
                            x = 'x', 
                            y = 'y',
                            ax = axs_APcats,
                            estimator=None,
                            color = colors_dict['primecolor'],
                            alpha = 0.3)

for l in violins1.lines:
    l.set_color(colors_dict['primecolor'])

[violin.set_edgecolor(colors_dict['primecolor']) for violin in violins1.collections]


axs_APcats.set_ylabel(f'{measure_str}')
        
# axs_APcats.set_ylim([0,5])

axs_APcats.grid(False)

plt.show()



save_figures(fig_APcats, f'{measure_str}_cell_all_freq_cats', figure_dir, darkmode_bool)























