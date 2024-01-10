# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 18:14:12 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import directories_win as directories
from PGFs import cc_APs_parameters
import pandas as pd
from useful_functions import calc_time_series, butter_filter
import os
from cc_IF_functions import get_IF_data
from plotting_functions import get_colors, save_figures, get_figure_size, set_font_sizes


# %%

table = pd.read_excel(directories.table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')


frequencies = ['1Hz', '30Hz', '75Hz']

# loop to create string to include all frequencies in query
query_str = ''

for idx, frequency in enumerate(frequencies):
    PGF = 'cc_APs_' + frequency
    
    if idx > 0:
        query_str = query_str + ' and '
        
    query_str = query_str + f'{PGF}.notnull()'
    

# limit lookup table
lookup_table = table.query(query_str)

# %% 

# test cell E-092
cell_ID = 'E-092'


vf_dict = dict.fromkeys(frequencies)
t_dict = dict.fromkeys(frequencies)

for frequency in frequencies:
    
    # PGF to load
    PGF = 'cc_APs_' + frequency
    PGF_parameters = cc_APs_parameters[frequency]
    
    # lookup_table = table.query(f'{PGF}.notnull()')
    
    
    # get indices of current cell with the dataframe containing all indices    
    group_idx = int(lookup_table.at[cell_ID, 'group'])-1
    series_idx = int(lookup_table.at[cell_ID, f'{PGF}'])-1
    
    # construct traceIndex with indices
    traceIndex = [group_idx, series_idx, 0, 0]
    
    # call on data file with indices from dataframe above
    current_file = lookup_table.at[cell_ID, 'file']
    
    data_file_path = os.path.join(directories.raw_data_dir, current_file + '.dat')
    
    data_file_path_str = fr"{data_file_path}"
    
    # get IF data form file
    i, v, t, SR, n_steps = get_IF_data(data_file_path_str, traceIndex, 'ms')
    
    # sampling rate in ms
    SR_ms = SR / 1e3
    
    # concatenate individual steps
    n_points = int(np.shape(i)[0] * np.shape(i)[1])
    
    i_concat = i.flatten() #.reshape([n_points], order='F')
    
    v_concat = v.flatten() #.reshape([n_points], order='F')
    
    t_concat = calc_time_series(v_concat, SR)
    
    # plt.plot(v_concat[0:15000])
    # plt.show()
    
    
    # filter voltage (to vf)
    vf = butter_filter(v_concat, 
                       order = 3,
                       cutoff = 1e3,
                       sampling_rate = SR)
    
    vf_dict[frequency] = vf
    t_dict[frequency] = t_concat


# %%

dm_bool = True

colors_dict = get_colors(dm_bool)

set_font_sizes(small_font_size = 14, large_font_size = 16)


fig_spiketrains, axs_spiketrains = plt.subplots(3, 3,
                                                layout = 'constrained',
                                                sharey = 'col',
                                                dpi = 600,
                                                figsize = get_figure_size(),
                                                gridspec_kw={'width_ratios': [2.5,5,3]})

v_range = [-100, 50]

## full time frame

axs_spiketrains[0][1].plot(t_dict['1Hz'], vf_dict['1Hz'])
axs_spiketrains[1][1].plot(t_dict['30Hz'], vf_dict['30Hz'])
axs_spiketrains[2][1].plot(t_dict['75Hz'], vf_dict['75Hz'])


# y axis formatting
for col in np.arange(1, 3, 1):
    axs_spiketrains[0][col].set_ylim(v_range)
    axs_spiketrains[0][col].set_yticks(np.arange(v_range[0], v_range[1] + 1, 50))
    axs_spiketrains[0][col].set_yticks(np.arange(v_range[0], v_range[1], 25), minor = True)
    
    for i, frequency in enumerate(frequencies):
        axs_spiketrains[i][col].set_ylabel('Voltage [mV]')
    
# x axis formatting
# 1 Hz
axs_spiketrains[0][1].set_xlim(0, t_dict['1Hz'][-1]+1)
axs_spiketrains[0][1].set_xticks(ticks = np.arange(0, t_dict['1Hz'][-1]+1, 20000),
                                 labels = np.arange(0, 100 + 1, 20, dtype = int))
axs_spiketrains[0][1].set_xticks(ticks = np.arange(0, t_dict['1Hz'][-1]+1, 5000),
                                 minor = True)

# 30 Hz
axs_spiketrains[1][1].set_xlim(0, t_dict['30Hz'][-1]+1)
axs_spiketrains[1][1].set_xticks(ticks = np.arange(0, t_dict['30Hz'][-1]+1, 1000),
                                 labels = np.arange(0, (t_dict['30Hz'][-1] / 1e3), 1, dtype = int))
axs_spiketrains[1][1].set_xticks(ticks = np.arange(0, t_dict['30Hz'][-1]+1, 500),
                                 minor = True)

# 75 Hz
axs_spiketrains[2][1].set_xlim(0, t_dict['75Hz'][-1]+1)
axs_spiketrains[2][1].set_xticks(ticks = np.arange(0, t_dict['75Hz'][-1]+1, 1000),
                                 labels = np.arange(0, (t_dict['75Hz'][-1] / 1e3), 1, dtype = int))
axs_spiketrains[2][1].set_xticks(ticks = np.arange(0, t_dict['75Hz'][-1]+1, 50),
                                 minor = True)

axs_spiketrains[2][1].set_xlabel('Time [s]')
axs_spiketrains[2][1].set_xlabel('Time [ms]')


## insets markers

for i, frequency in enumerate(frequencies):
    
    inset_start = 0.002 * 1e3   # ms
    inset_length = 0.098 * 1e3  # ms
    
    if i == 0:
        inset_start = 500 - inset_length / 2
    
    inset_marker = mtl.patches.Rectangle(xy = (inset_start, v_range[0]*0.99), 
                                         width=inset_length, 
                                         height=(np.abs(v_range[0]-v_range[1]))*0.99,
                                         facecolor = 'none',
                                         edgecolor='r',
                                         linewidth = 3)
    
    axs_spiketrains[i][1].add_patch(inset_marker)

## inset time frame

    start_idx = int(inset_start * SR_ms)
    stop_idx = int((inset_start + inset_length) * SR_ms)
    
    v_inset = vf_dict[frequency][start_idx:stop_idx]
    t_inset = t_dict[frequency][start_idx:stop_idx]

    axs_spiketrains[i][2].plot(t_inset, v_inset)
    

    axs_spiketrains[i][2].set_xlim(inset_start, inset_start + inset_length)
    
    if i == 0:
        axs_spiketrains[i][2].set_xticks(np.arange(450, 550 + 1, 50))
    elif i > 0:
        axs_spiketrains[i][2].set_xticks(np.arange(0, inset_length + 3, 50))


# remove ticks in first column

for i, frequency in enumerate(frequencies):
    axs_spiketrains[i][0].tick_params(axis = 'both',labelsize = 0,  size = 0)
    axs_spiketrains[i][0].grid(False)


save_figures(fig_spiketrains, 'stim_frequencies_examples', directories.figure_dir, dm_bool)

plt.show()























