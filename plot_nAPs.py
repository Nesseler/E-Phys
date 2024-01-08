# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 16:58:10 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import parameters
from PGFs import cc_APs_parameters
import directories_win as directories
from plotting_functions import save_figures, get_colors



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
n_APs_df_path = os.path.join(directories.quant_data_dir, 'APs', 'n_APs.xlsx')

n_APs_df.to_excel(n_APs_df_path)

    

# %%

dm_bool = False

colors_dict = get_colors(dm_bool)


fig_n_APs, axs_n_APs = plt.subplots(1, 1,
                                    layout = 'constrained',
                                    sharey='all', sharex='all',
                                    dpi = 600)

# fig_n_APs.set_constrained_layout_pads(wspace=0.05, w_pad=0.0,
#                                       hspace=0.05, h_pad=0.0) 



freqs = [int(freq.replace('Hz', '')) for freq in n_APs_df.index]

axs_n_APs.plot(freqs, n_APs_df, marker = '.')


axs_n_APs.set_ylim(0,105)
axs_n_APs.set_ylabel('Number of APs [#]')

axs_n_APs.set_xlim(0,76)
axs_n_APs.set_xticks(freqs)
axs_n_APs.set_xlabel('Stimulation frequency [Hz]')

#axs_n_APs.set_xscale('log')

plt.show()    

save_figures(fig_n_APs, 'nAPs', directories.figure_dir, dm_bool)