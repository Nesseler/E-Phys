# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 18:07:01 2024

@author: nesseler
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:27:54 2024

@author: nesseler
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtl
import seaborn as sbn


# custom directories & parameters
from PGFs import cc_APs_parameters, cc_APs_n_stims

from directories_win import cell_descrip_dir, quant_data_dir, figure_dir, vplot_dir

from functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size
from functions_export import set_df_to_cell_descrips


# %%

frequencies = list(cc_APs_parameters.keys())
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

# initialise dataframe to populate in loop
tAPs_df = pd.DataFrame()

# set directory to APs subfolder to get parameters of action potentials
quant_dir = os.path.join(quant_data_dir, 'APs')

# create a list of cell-IDs that have subfolder in the APs folder i.e. that have been analysed
cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]



def import_AP_measurement_all_freqs(frequencies, parameter, cell_ID):   

    # create file path
    cell_data_path = os.path.join(quant_dir, cell_ID)
    
    # get number of APs for each stim
    measurement_df = pd.DataFrame()
    
    for frequency in frequencies:
    # frequency = '1Hz'
            
        # create file path
        file_path = os.path.join(cell_data_path, f'{cell_ID}_{frequency}.xlsx')
        
        # read excel file
        AP_all_params = pd.read_excel(file_path, index_col = 'idx_step')
        
        # write to dictionary
        measurement_df[frequency] = AP_all_params[parameter]

    return measurement_df



# %% get data

# start with one cell
cell_ID = 'E-092'

mean_ISIs_df = pd.DataFrame(index = frequencies)

for cell_ID in cell_IDs:
    t_APs = import_AP_measurement_all_freqs(frequencies, 't_peaks', cell_ID)
      
    # create list with interstimulus interval per stimulation frequency
    ISIs = [sum(cc_APs_parameters[f].values()) for f in frequencies]
    total_dur = [ISI * cc_APs_n_stims for ISI in ISIs]
    
    # create dataframe with all stimuation time of every frequency
    t_stims_df = pd.DataFrame()
    
    for idx_freq, freq in enumerate(frequencies):
        
        t_pre = cc_APs_parameters[freq]['t_pre']
        
        t_stims = [t_pre + ISIs[idx_freq] * i for i in np.arange(cc_APs_n_stims)]
    
        t_stims_df[freq] = t_stims
    
    # add both dataframes and create spike times in continues timeseries
    t_APs_cont = t_APs.add(t_stims_df)
    
    
    # calc mean ISI from continues spike times
    mean_ISIs_df[cell_ID] = [np.mean(t_APs_cont[freq].dropna().diff()) for freq in frequencies]

  
# %% convert mean ISI to frequency

# ms to s
resul_freq_df = mean_ISIs_df.div(1000)

# ISI to freq
resul_freq_df = resul_freq_df.rdiv(1)

# get maximal resulting frequency
max_resul_freq_df = resul_freq_df.max(axis = 0)

# %% export measurements 

# save measurements to excel file
mean_ISIs_df.to_excel(os.path.join(cell_descrip_dir, 'mean_ISIs.xlsx'), index_label = 'frequencies')
resul_freq_df.to_excel(os.path.join(cell_descrip_dir, 'resul_freq.xlsx'), index_label = 'frequencies')

# %% plotting

darkmode_bool = True

colors_dict = get_colors(darkmode_bool)

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










