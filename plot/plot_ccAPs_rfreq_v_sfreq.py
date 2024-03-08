# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 18:07:01 2024

@author: nesseler
"""


import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtl
import seaborn as sbn


# custom directories & parameters
from parameters.directories_win import cell_descrip_dir, figure_dir, table_file
from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size


# %% get data

mean_ISIs_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccAPs-mean_ISIs.xlsx'), index_col = 'frequencies')
resul_freq_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccAPs-resul_freq.xlsx'), index_col = 'frequencies')

# get maximal resulting frequency
max_resul_freq_df = resul_freq_df.max(axis = 0)


frequencies = list(mean_ISIs_df.index)
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

cell_IDs = list(mean_ISIs_df.columns)


# %% load Metadata

MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

# %% plotting

darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)

set_font_sizes()

fig_freq, axs_freq = plt.subplots(1,2,
                                  layout = 'constrained',
                                  figsize = get_figure_size(),
                                  gridspec_kw = {'width_ratios': [6,2]})


# initialise color code
norm_min = 0
norm_max = 75
cmap_str = 'viridis'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# colorbar
fig_freq.colorbar(cmap, ax = axs_freq[0])

# add unity line
axs_freq[0].axline(xy1 = (0,0), slope = 1, 
                   c = 'gray', 
                   linestyle = '--')

# plot resulted frequency vs stimulated frequency with colorcode
for cell_ID in cell_IDs:
    axs_freq[0].plot(freqs_int, resul_freq_df[cell_ID],
                     marker = 'o',
                     c = cmap.to_rgba(max_resul_freq_df.at[cell_ID]))


axs_freq[0].set_xticks(freqs_int)
axs_freq[0].set_xlim([0,80])
axs_freq[0].set_xlabel('Stimulated frequency [Hz]')

axs_freq[0].set_yticks(np.arange(0, 80, 25))
axs_freq[0].set_yticks(np.arange(0, 80, 5), minor = True)
axs_freq[0].set_ylim([0,80])
axs_freq[0].set_ylabel('Resulting frequency [Hz]')


# max frequency
violin = sbn.violinplot(data = max_resul_freq_df,
                        width = .8,
                        inner = 'quart',
                        linewidth = 1.5,
                        color = colors_dict['seccolor'],
                        ax = axs_freq[1])

[c.set_edgecolor(colors_dict['primecolor']) for c in violin.collections[:2]]

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

sbn.swarmplot(data = max_resul_freq_df,
              c = max_resul_freq_df, 
              cmap = cmap_str, 
              norm = norm,
              ax = axs_freq[1], 
              size = 8)

axs_freq[1].set_xticks([])
axs_freq[1].set_ylim([0,80])
axs_freq[1].set_ylabel('Maximal resulting frequency [Hz]')


[ax.grid(False) for ax in axs_freq]

plt.show()

save_figures(fig_freq, 'rfreq_v_sfreq', figure_dir, darkmode_bool)


# %% plotting + region

darkmode_bool = True

set_font_sizes()

fig_freq, axs_freq = plt.subplots(1,2,
                                  layout = 'constrained',
                                  figsize = get_figure_size(),
                                  gridspec_kw = {'width_ratios': [6,2]})


# initialise color code
norm_min = 0
norm_max = 75
cmap_str = 'viridis'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# colorbar
fig_freq.colorbar(cmap, ax = axs_freq[0])

# add unity line
axs_freq[0].axline(xy1 = (0,0), slope = 1, 
                   c = 'gray', 
                   linestyle = '--')

# plot resulted frequency vs stimulated frequency with colorcode
for cell_ID in cell_IDs:
    axs_freq[0].plot(freqs_int, resul_freq_df[cell_ID],
                     marker = 'o',
                     c = region_colors[MetaData.at[cell_ID, 'Region']])


axs_freq[0].set_xticks(freqs_int)
axs_freq[0].set_xlim([0,80])
axs_freq[0].set_xlabel('Stimulated frequency [Hz]')

axs_freq[0].set_yticks(np.arange(0, 80, 25))
axs_freq[0].set_yticks(np.arange(0, 80, 5), minor = True)
axs_freq[0].set_ylim([0,80])
axs_freq[0].set_ylabel('Resulting frequency [Hz]')


# max frequency
violin = sbn.violinplot(data = max_resul_freq_df,
                        width = .8,
                        inner = 'quart',
                        linewidth = 1.5,
                        color = colors_dict['seccolor'],
                        ax = axs_freq[1])

[c.set_edgecolor(colors_dict['primecolor']) for c in violin.collections[:2]]

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

sbn.swarmplot(data = max_resul_freq_df,
              c = max_resul_freq_df, 
              cmap = cmap_str, 
              norm = norm,
              ax = axs_freq[1], 
              size = 8)

axs_freq[1].set_xticks([])
axs_freq[1].set_ylim([0,80])
axs_freq[1].set_ylabel('Maximal resulting frequency [Hz]')


[ax.grid(False) for ax in axs_freq]

plt.show()







