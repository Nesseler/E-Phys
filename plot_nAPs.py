# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 16:58:10 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import pandas as pd
import numpy as np
import os
import parameters
from PGFs import cc_APs_parameters
import directories_win as directories
from plotting_functions import save_figures, get_colors, get_colorcode, get_figure_size, set_font_sizes
import seaborn as sbn


from directories_win import cell_descrip_dir



# %% 

frequencies = list(cc_APs_parameters.keys())

n_APs_df = pd.DataFrame()

quant_dir = os.path.join(directories.quant_data_dir, 'APs')

cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]

def import_all_freq_onecell(cell_ID, frequencies):
    
    # create file path
    cell_data_path = os.path.join(os.path.join(directories.quant_data_dir, 'APs'), cell_ID)
    
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



# write excel file
n_APs_df_path = os.path.join(cell_descrip_dir, 'n_APs.xlsx')

n_APs_df.to_excel(n_APs_df_path, index_label='frequency')


# %%

# dataframe for mean follow-index
mean_fIndex = pd.DataFrame()

for cell_ID in n_APs_df.columns:
    mean_fIndex.loc[cell_ID ,'fIndex'] = n_APs_df[cell_ID].mean()


    

# %%

dm_bool = True

colors_dict = get_colors(dm_bool)

set_font_sizes()


fig_n_APs, axs_n_APs = plt.subplots(1, 2,
                                    layout = 'constrained',
                                    dpi = 600,
                                    figsize = get_figure_size(),
                                    gridspec_kw={'width_ratios': [8,2]})

# fig_n_APs.set_constrained_layout_pads(wspace=0.05, w_pad=0.0,
#                                       hspace=0.05, h_pad=0.0) 

# create list of integers for stimulation frequencies
freqs = [int(freq.replace('Hz', '')) for freq in n_APs_df.index]

# create normalization vector with min and max
# norm = mtl.colors.CenteredNorm(mean_fIndex.min(), mean_fIndex.max())
# norm = mtl.colors.Normalize(0, 100)

# lc = get_colorcode(freqs, n_APs_df, mean_fIndex['fIndex'], norm = norm, cmap='seismic')

# line = axs_n_APs[0].add_collection(lc)

norm_min = 40
norm_max = 100
cmap_str = 'plasma'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
cmap.set_array(list(mean_fIndex['fIndex']))


fig_n_APs.colorbar(cmap, 
                   ax = axs_n_APs[0])

for i, cell_ID in enumerate(n_APs_df.columns):
    axs_n_APs[0].plot(freqs, n_APs_df[cell_ID], 
                      marker = '.', 
                      c=cmap.to_rgba(mean_fIndex.at[cell_ID, 'fIndex']),
                      markersize = 8)


axs_n_APs[0].set_ylim(0,105)
axs_n_APs[0].set_ylabel('Number of APs [#]')

axs_n_APs[0].set_xlim(0,76)
axs_n_APs[0].set_xticks(freqs)
axs_n_APs[0].set_xlabel('Stimulation frequency [Hz]')



# axs_n_APs[1].scatter([1]*len(n_APs_df.columns), mean_fIndex['fIndex'], c=mean_fIndex['fIndex'], cmap=cmap_str, norm = norm)


violin = sbn.violinplot(x = [1.]*len(mean_fIndex.index), 
                        y = mean_fIndex['fIndex'],
                        width = 1.3,
                        inner = 'quart',
                        linewidth = 1.5,
                        color = colors_dict['seccolor'],
                        ax = axs_n_APs[1])

violin.collections[0].set_edgecolor(colors_dict['primecolor'])

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])



sbn.swarmplot(x = [1.]*len(mean_fIndex.index), 
              y = mean_fIndex['fIndex'],
              c = mean_fIndex['fIndex'], 
              cmap = cmap_str, 
              norm = norm,
              ax = axs_n_APs[1], 
              size = 8)



axs_n_APs[1].set_ylim(0, 105)
axs_n_APs[1].set_ylabel('Mean Number of APs')

axs_n_APs[1].set_xlim(-1,1)
axs_n_APs[1].set_xticks([])
# axs_n_APs[1].set_xlabel('Stimulation frequency [Hz]')


plt.show()    

save_figures(fig_n_APs, 'nAPs', directories.figure_dir, dm_bool)


# %%





