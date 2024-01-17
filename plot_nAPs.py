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

from PGFs import cc_APs_parameters

from functions_plotting import save_figures, get_colors, get_figure_size, set_font_sizes
from directories_win import cell_descrip_dir, quant_data_dir, figure_dir
from functions_export import set_df_to_cell_descrips


# %% 

frequencies = list(cc_APs_parameters.keys())

n_APs_df = pd.DataFrame()

quant_dir = os.path.join(quant_data_dir, 'APs')

cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]

def import_all_freq_onecell(cell_ID, frequencies):
    
    # create file path
    cell_data_path = os.path.join(os.path.join(quant_data_dir, 'APs'), cell_ID)
    
    # get number of APs for each stim
    n_APs = {}
    
    for frequency in frequencies:
        
        # create file path
        file_path = os.path.join(cell_data_path, f'{cell_ID}_{frequency}.xlsx')
        
        # read excel file
        AP_all_params = pd.read_excel(file_path)
        
        # run query for where spikes were detected
        only_AP_all_params = AP_all_params.query('v_peaks.notnull()')
        
        # get number of spikes
        n_APs_this_frequency = len(only_AP_all_params.index)
        
        # write to dictionary
        n_APs[frequency] = n_APs_this_frequency
    
    return n_APs



for cell_ID in cell_IDs:
    n_APs = import_all_freq_onecell(cell_ID, frequencies)
    
    n_APs_df[cell_ID] = n_APs

nAPs_path = os.path.join(cell_descrip_dir, 'nAPs.xlsx')
n_APs_df.to_excel(nAPs_path, index_label='cell_ID')


# %%

# dataframe for mean follow-index
fIndex = pd.DataFrame()
fIndex.index.name = 'cell_ID'


for cell_ID in n_APs_df.columns:
    fIndex.loc[cell_ID ,'mean'] = n_APs_df[cell_ID].mean()
    fIndex.loc[cell_ID ,'std'] = n_APs_df[cell_ID].std()



fIndex_melted = fIndex.melt(var_name = 'measurement')
    

# %% add mean and std to n_APs_df

# write excel file
set_df_to_cell_descrips(n_APs_df, add_header_ext = 'nAPs')

set_df_to_cell_descrips(fIndex, add_header_ext='nAPs')


# %%

dm_bool = True

colors_dict = get_colors(dm_bool)

set_font_sizes()


fig_n_APs, axs_n_APs = plt.subplots(1, 2,
                                    layout = 'constrained',
                                    dpi = 600,
                                    figsize = get_figure_size(),
                                    gridspec_kw={'width_ratios': [6,4]})

# fig_n_APs.set_constrained_layout_pads(wspace=0.05, w_pad=0.0,
#                                       hspace=0.05, h_pad=0.0) 

# create list of integers for stimulation frequencies
freqs = [int(freq.replace('Hz', '')) for freq in n_APs_df.index]

# create normalization vector with min and max
# norm = mtl.colors.CenteredNorm(mean_fIndex.min(), mean_fIndex.max())
# norm = mtl.colors.Normalize(0, 100)

# lc = get_colorcode(freqs, n_APs_df, mean_fIndex['fIndex'], norm = norm, cmap='seismic')

# line = axs_n_APs[0].add_collection(lc)

norm_min = 0
norm_max = 50
cmap_str = 'plasma'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

cmap.set_array(list(fIndex['mean']))


fig_n_APs.colorbar(cmap, 
                   ax = axs_n_APs[0])

for i, cell_ID in enumerate(n_APs_df.columns):
    axs_n_APs[0].plot(freqs, n_APs_df[cell_ID], 
                      marker = '.', 
                      c=cmap.to_rgba(fIndex.at[cell_ID, 'mean']),
                      markersize = 8)


axs_n_APs[0].set_ylim(0,105)
axs_n_APs[0].set_ylabel('Number of APs [#]')

axs_n_APs[0].set_xlim(0,76)
axs_n_APs[0].set_xticks(freqs)
axs_n_APs[0].set_xlabel('Stimulation frequency [Hz]')



# axs_n_APs[1].scatter([1]*len(n_APs_df.columns), mean_fIndex['fIndex'], c=mean_fIndex['fIndex'], cmap=cmap_str, norm = norm)


# violin = sbn.violinplot(data = fIndex, 
#                         y = 'mean_fIndex',
#                         width = 1.3,
#                         inner = 'quart',
#                         linewidth = 1.5,
#                         color = colors_dict['seccolor'],
#                         ax = axs_n_APs[1])

# violin.collections[0].set_edgecolor(colors_dict['primecolor'])

# for l in violin.lines:
#     l.set_color(colors_dict['primecolor'])



# sbn.swarmplot(x = [1.]*len(mean_fIndex.index), 
#               y = mean_fIndex['fIndex'],
#               c = mean_fIndex['fIndex'], 
#               cmap = cmap_str, 
#               norm = norm,
#               ax = axs_n_APs[1], 
#               size = 8)


violin = sbn.violinplot(data = fIndex_melted,
                        x = 'measurement',
                        y = 'value',
                        width = .9,
                        inner = 'quart',
                        linewidth = 1.5,
                        color = colors_dict['seccolor'],
                        ax = axs_n_APs[1])

[c.set_edgecolor(colors_dict['primecolor']) for c in violin.collections[:2]]

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])



sbn.swarmplot(data = fIndex_melted,
              x = 'measurement',
              y = 'value',
              c = fIndex['mean'], 
              cmap = cmap_str, 
              norm = norm,
              ax = axs_n_APs[1], 
              size = 8)


axs_n_APs[1].set_ylim(0, 105)
axs_n_APs[1].set_ylabel('Mean Number of APs')

axs_n_APs[1].set_xlim(-0.5,1.5)
axs_n_APs[1].set_xticks([])
# axs_n_APs[1].set_xlabel('Stimulation frequency [Hz]')


plt.show()    

save_figures(fig_n_APs, 'nAPs', figure_dir, dm_bool)


# %% split mean and std


dm_bool = True

colors_dict = get_colors(dm_bool)

set_font_sizes()


fig_APs, axs_APs = plt.subplots(1, 3,
                                    layout = 'constrained',
                                    dpi = 600,
                                    figsize = get_figure_size(),
                                    gridspec_kw={'width_ratios': [6,2,2]})

# initialise color code
norm_min = 0
norm_max = 50
cmap_str = 'plasma'
c_str = 'std'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
cmap.set_array(list(fIndex[c_str]))



# colorbar
fig_APs.colorbar(cmap, ax = axs_APs[0])

# plot lines for all cells
for i, cell_ID in enumerate(n_APs_df.columns):
    axs_APs[0].plot(freqs, n_APs_df[cell_ID], 
                      marker = '.', 
                      c=cmap.to_rgba(fIndex.at[cell_ID, c_str]),
                      markersize = 8)

# axis configuration
axs_APs[0].set_ylim(0,105)
axs_APs[0].set_ylabel('Number of APs [#]')

axs_APs[0].set_xlim(0,76)
axs_APs[0].set_xticks(freqs)
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

# axs_APs[1].set_xlim(-0.5,1)
axs_APs[1].set_xticks([])
# axs_n_APs[1].set_xlabel('Stimulation frequency [Hz]')


save_figures(fig_APs, f'nAPs-{c_str}_ccoded', figure_dir, dm_bool)