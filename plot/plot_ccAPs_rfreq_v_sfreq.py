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
max_resul_freq_df.name = 'max_resul_freq'


frequencies = list(mean_ISIs_df.index)
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

cell_IDs = list(mean_ISIs_df.columns)


# %% load Metadata

MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

max_resul_freq_n_metadata_df = pd.concat([max_resul_freq_df, MetaData['Region']], axis = 1)

# %% plotting

darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)

set_font_sizes(small_font_size=12)

fig_freq, axs_freq = plt.subplots(1,2,
                                  layout = 'constrained',
                                  figsize = get_figure_size(width = 277.25, height = 139.607),
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
                        bw  = 0.3,
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
              size = 3)

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
violins = sbn.violinplot(data = max_resul_freq_n_metadata_df,
                        y = 'max_resul_freq',
                        width = .8,
                        inner = 'quart',
                        linewidth = 1.5,
                        color = colors_dict['seccolor'],
                        ax = axs_freq[1])

for l in violins.lines:
    l.set_color(colors_dict['primecolor'])

[violin.set_edgecolor(colors_dict['primecolor']) for violin in violins.collections]

sbn.swarmplot(data = max_resul_freq_n_metadata_df,
              y = 'max_resul_freq', 
              ax = axs_freq[1], 
              size = 8,
              hue = 'Region',
              palette = region_colors)

# axs_freq[1].set_xticks([])
axs_freq[1].set_ylim([0,80])
axs_freq[1].set_ylabel('Maximal resulting frequency [Hz]')


[ax.grid(False) for ax in axs_freq]

plt.show()
save_figures(fig_freq, 'rfreq_v_sfreq+regions', figure_dir, darkmode_bool)


# %% plot + separate panels for regions


# initialize figure
ax_keys = ['A', 'B', 'C', 'D']
regions = ['BAOT/MeA', 'MeA', 'BAOT']

fig_regions, axs_regions = plt.subplot_mosaic('AD;BD;CD', 
                                              layout = 'constrained',
                                              figsize = get_figure_size(width = 277.25, height = 137.607),
                                              width_ratios = [2.5, 1],
                                              dpi = 600)

# set_font_sizes(small_font_size=12)

# sfreq vs rfreq

for idx_region, region in enumerate(regions):
    
    # get cell_IDs of dataframe per region
    cell_IDs_region = MetaData[MetaData['Region'] == region].index.to_list()
    
    # get number of cells in region
    n_cells_region = len(cell_IDs_region)
    
    # set axis 
    ax = axs_regions[ax_keys[idx_region]]
    
    # set region name as title
    ax.set_title(region, fontsize = 18)
    
    # add unity line
    ax.plot([1, 75], [1, 75],
            c = 'gray', 
            linestyle = '--',
            lw = 1.5)
    
    for idx_cell, cell_ID in enumerate(cell_IDs_region):
        # plot resulted frequency vs stimulated frequency with colorcode
        # for cell_ID in cell_IDs_region:
        ax.plot(freqs_int, resul_freq_df[cell_ID],
                marker = 'o',
                c = region_colors[MetaData.at[cell_ID, 'Region']],
                markersize = 3,
                lw = .5)

    ax.set_xlim([-1, 80])
    ax.set_xticks(ticks = freqs_int,
                  labels = freqs_int)
    
    ax.set_ylim([-1, 80])
    ax.set_yticks(ticks = np.arange(0, 75+1, 25), 
                  labels = np.arange(0, 75+1, 25))
    ax.set_yticks(np.arange(0, 75+1, 5), minor = True)
    
    # limit spines
    ax.spines['left'].set_bounds([1, 75])
    ax.spines['bottom'].set_bounds([1, 75])

    # despine
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]



# remove ticks between first two activity plots
for i in ['A', 'B']:
    axs_regions[i].set_xticks(ticks = freqs_int, labels = [])
    
axs_regions['C'].set_xlabel('Stimulated frequency [Hz]')  

fig_regions.supylabel('Resulting frequency [Hz]', fontsize = 12)
    
# max frequency
violins = sbn.violinplot(data = max_resul_freq_n_metadata_df,
                        y = 'max_resul_freq',
                        x  = 'Region',
                        width = .9,
                        bw = 0.5,
                        inner = 'quart',
                        linewidth = 1.5,
                        color = colors_dict['seccolor'],
                        ax = axs_regions['D'],
                        bw_adjust=.5,
                        density_norm = "area")

for l in violins.lines:
    l.set_color(colors_dict['primecolor'])

[violin.set_edgecolor(colors_dict['primecolor']) for violin in violins.collections]

swarms = sbn.swarmplot(data = max_resul_freq_n_metadata_df,
                       y = 'max_resul_freq',
                       x  = 'Region', 
                       ax = axs_regions['D'], 
                       size = 6,
                       hue = 'Region',
                       palette = region_colors,
                       order = ['BAOT/MeA', 'MeA', 'BAOT'])

# axs_freq[1].set_xticks([])
axs_regions['D'].set_ylim([0,75])
axs_regions['D'].set_ylabel('Maximal resulting frequency [Hz]')
axs_regions['D'].set_yticks(np.arange(0, 75+1, 25))
axs_regions['D'].set_yticks(np.arange(0, 75+1, 5), minor = True)
axs_regions['D'].set_xticklabels(['BAOT/\nMeA', 'MeA', 'BAOT'])

axs_regions['D'].set_xlabel('')


# limit spines
axs_regions['D'].spines['left'].set_bounds([0, 75])
axs_regions['D'].spines['bottom'].set_bounds([0, 2])

# despine
[axs_regions['D'].spines[spine].set_visible(False) for spine in ['top', 'right']]

[axs_regions[ax_keys].grid(False) for ax_keys in axs_regions]

axs_regions['D'].get_legend().remove()

# set font sizes
small_font_size = 12

plt.rc('font', size = small_font_size)
plt.rc('axes', titlesize = small_font_size, 
               labelsize = small_font_size,
               linewidth = 1)
plt.rc('xtick', labelsize = small_font_size)
plt.rc('ytick', labelsize = small_font_size)
plt.rc('lines', linewidth = 1)

fig_regions.align_labels() 


plt.show()

save_figures(fig_regions, 'rfreq_v_sfreq+sep_regions', figure_dir, darkmode_bool,
             figure_format= 'both',
             dataframe_to_save = max_resul_freq_n_metadata_df, index_label = 'cell_ID', add_measures = True, axis_for_calcs = 0,
             groups_bool= True, groups= ['BAOT/MeA', 'MeA', 'BAOT'], groups_name= 'Region')


resul_freq_df.to_excel(os.path.join(figure_dir, 'rfreq_v_sfreq.xlsx'), index_label = 'cell_ID')  









