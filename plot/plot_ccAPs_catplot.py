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
from parameters.directories_win import cell_descrip_dir, figure_dir, table_file

# custom functions
from functions.functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes


# %% import 

# FWHM, tpeaks, vamplitude
measure_str = 'mean_' + 'vamplitude'

mean_df = pd.read_excel(os.path.join(cell_descrip_dir, f'ccAPs-{measure_str}.xlsx'), index_col='frequencies')

cell_IDs = list(mean_df.columns)

frequencies = list(mean_df.index)
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]


# %% load Metadata

MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]


# %% transpose and melt Dataframe to be used with seaborn functions

mean_df = mean_df.transpose()
mean_df_melt = pd.melt(mean_df, var_name='frequency', value_name=measure_str, ignore_index= False)     


# %% color code

# norm_min = 1
# norm_max = 75
# cmap_str = 'plasma'

# norm = mtl.colors.Normalize(norm_min, norm_max)
# cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
# cmap.set_array(freqs_int)


# # creating specific color pallet for seaborn plotting functions
# color_pal = {'1Hz'  : cmap.to_rgba(freqs_int[0]), 
#              '5Hz'  : cmap.to_rgba(freqs_int[1]),
#              '10Hz' : cmap.to_rgba(freqs_int[2]),
#              '30Hz' : cmap.to_rgba(freqs_int[3]),
#              '50Hz' : cmap.to_rgba(freqs_int[4]),
#              '75Hz' : cmap.to_rgba(freqs_int[5])}

# %% create dataframe for seaborn to handle


region_df = MetaData['Region']

plt_df = mean_df

plt_df = plt_df.transpose()



plt_df_melted = pd.melt(plt_df, var_name='cell_ID', value_name=measure_str, ignore_index= False) 
plt_df_melted.index.name = 'frequency'

plt_df_melted = plt_df_melted.reset_index().set_index('cell_ID')


plt_df_melted['Region'] = region_df



# %% plotting

darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)

regions = ['BAOT/MeA', 'MeA', 'BAOT']

fig_APcats, axs_APcats = plt.subplots(nrows = 1, 
                                      ncols = 1,
                                      layout = 'constrained',
                                      dpi = 600,
                                      figsize = get_figure_size(),
                                      sharex = 'col')

set_font_sizes()


# fig_APcats.colorbar(cmap, 
#                     ax = axs_APcats,
#                     ticks = freqs_int,
#                     label = 'Stimulation frequency [Hz]')


violins = sbn.violinplot(data = plt_df_melted, 
                          x = 'frequency', 
                          y = measure_str, 
                          inner = 'quart',
                          hue = 'Region',
                          hue_order = regions,
                          ax = axs_APcats, 
                          linewidth = 1.5,
                          palette = ['k', 'k', 'k'])

# 6 stimulation frequencies, 3 regions and 3 lines (for quartiles) in violins
region_per_violin = []
[region_per_violin.append(region) for freq in frequencies for region in regions for i in range(3)]

for idx_l, l in enumerate(violins.lines):
    l.set_color(region_colors[region_per_violin[idx_l]])

# 6 stimulation frequencies, 3 regions and 3 lines (for quartiles) in violins
region_per_violin = []
[region_per_violin.append(region) for freq in frequencies for region in regions]

for idx_violin, violin in enumerate(violins.collections):
    violin.set_edgecolor(region_colors[region_per_violin[idx_violin]])



swarms = sbn.swarmplot(plt_df_melted,
                       x = 'frequency',
                       y = measure_str,
                       hue = 'Region',
                       palette = region_colors,
                       dodge = True,
                       size = 5,
                       hue_order = regions)


for region_idx, region in enumerate(regions):
    
    print(region_idx, region)
    
    print(MetaData[MetaData['Region'] == region].index.to_list())
    
    region_cell_IDs = MetaData[MetaData['Region'] == region].index.to_list()

    # get positions of all points in swarmplots to plot the connecting lines
    positions_df = pd.DataFrame()
    
    for idx, frequency in enumerate(frequencies):
        positions = np.array(axs_APcats.collections[(idx*3)+region_idx+18].get_offsets())
        
        cur_df = pd.DataFrame({'x' : positions[:, 0],
                                'y' : positions[:, 1],
                                'frequency' : [freqs_int[idx]] * len(positions),
                                'freq_str' : [frequency] * len(positions),
                                'cell_ID' : region_cell_IDs
                                })
    
        positions_df = pd.concat([positions_df, cur_df])
        


    for idx, cell_ID in enumerate(region_cell_IDs):
        lines = sbn.lineplot(data = positions_df[positions_df.cell_ID == cell_ID], 
                              x = 'x', 
                              y = 'y',
                              ax = axs_APcats,
                              estimator=None,
                              color = region_colors[region],
                              alpha = 0.3)



axs_APcats.set_xlabel('Stimulation frequency [Hz]')
axs_APcats.set_ylabel(f'{measure_str}')

axs_APcats.grid(False)

axs_APcats.get_legend().remove()

plt.show()



save_figures(fig_APcats, f'{measure_str}_cell_all_freq_cats', figure_dir, darkmode_bool)























