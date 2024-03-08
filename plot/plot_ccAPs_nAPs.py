# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 16:58:10 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import pandas as pd
import os
import seaborn as sbn

# custom parameters
from parameters.directories_win import cell_descrip_dir, figure_dir

# custom functions
from functions.functions_plotting import save_figures, get_colors, get_figure_size, set_font_sizes



# %% import data frame

nAPs_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccAPs-nAPs.xlsx'), index_col = 'frequencies')

cell_IDs = list(nAPs_df.columns)

frequencies = list(nAPs_df.index)
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

# %%

# dataframe for mean follow-index
fIndex = pd.DataFrame()
fIndex.index.name = 'cell_ID'


for cell_ID in cell_IDs:
    fIndex.loc[cell_ID ,'mean'] = nAPs_df[cell_ID].mean()
    fIndex.loc[cell_ID ,'std'] = nAPs_df[cell_ID].std()



fIndex_melted = fIndex.melt(var_name = 'measurement')
    


# %% split mean and std


dm_bool = True

colors_dict, region_colors = get_colors(dm_bool)

set_font_sizes()


fig_APs, axs_APs = plt.subplots(1, 3,
                                    layout = 'constrained',
                                    dpi = 600,
                                    figsize = get_figure_size(),
                                    gridspec_kw={'width_ratios': [6,2,2]})

# initialise color code
norm_min = 40
norm_max = 100
cmap_str = 'plasma'
c_str = 'mean'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
cmap.set_array(list(fIndex[c_str]))



# colorbar
fig_APs.colorbar(cmap, ax = axs_APs[0])

# plot lines for all cells
for i, cell_ID in enumerate(cell_IDs):
    axs_APs[0].plot(freqs_int, nAPs_df[cell_ID], 
                      marker = '.', 
                      c=cmap.to_rgba(fIndex.at[cell_ID, c_str]),
                      markersize = 8)

# axis configuration
axs_APs[0].set_ylim(0,105)
axs_APs[0].set_ylabel('Number of APs [#]')

axs_APs[0].set_xlim(0,76)
axs_APs[0].set_xticks(freqs_int)
axs_APs[0].set_xlabel('Stimulation frequency [Hz]')


# average n APs

violin = sbn.violinplot(data = fIndex,
                        y = 'mean',
                        width = .9,
                        inner = 'quart',
                        linewidth = 1.5,
                        color = colors_dict['seccolor'],
                        ax = axs_APs[1])

[c.set_edgecolor(colors_dict['primecolor']) for c in violin.collections[:2]]

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

sbn.swarmplot(data = fIndex,
              y = 'mean',
              c = fIndex[c_str], 
              cmap = cmap_str, 
              norm = norm,
              ax = axs_APs[1], 
              size = 8)


# std n APs

violin = sbn.violinplot(data = fIndex,
                        y = 'std',
                        width = .9,
                        inner = 'quart',
                        linewidth = 1.5,
                        color = colors_dict['seccolor'],
                        ax = axs_APs[2])


[c.set_edgecolor(colors_dict['primecolor']) for c in violin.collections[:2]]

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

sbn.swarmplot(data = fIndex,
              y = 'std',
              c = fIndex[c_str], 
              cmap = cmap_str, 
              norm = norm,
              ax = axs_APs[2], 
              size = 8)



axs_APs[1].set_ylim(0, 105)
axs_APs[1].set_ylabel('Mean Number of APs')

axs_APs[1].spines['right'].set_visible(False)
axs_APs[1].spines['top'].set_visible(False)

axs_APs[2].set_ylim(0, 105)
axs_APs[2].set_ylabel('Std Number of APs')

axs_APs[2].spines['right'].set_visible(False)
axs_APs[2].spines['top'].set_visible(False)

axs_APs[1].set_xticks([])


save_figures(fig_APs, f'nAPs-{c_str}_ccoded', figure_dir, dm_bool)