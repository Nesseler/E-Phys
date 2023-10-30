#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:02:40 2023

@author: moritznesseler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from useful_functions import calc_time_series, calc_dvdt, butter_filter, save_figures, get_sampling_rate, get_data
import os.path


from plotting_functions import get_figure_size



## n/3 by 3 subplots with each a 30 sec resting activity plot

## eventplot with all cells


lookup_table = pd.read_csv('/Users/moritznesseler/Desktop/local E-Phys/cc_rest.csv', 
                           delimiter=';',
                           index_col='cell_ID')


data_folder = '/Users/moritznesseler/Desktop/local E-Phys'



n_cells = len(lookup_table)
n_cols = 3
n_rows = int(np.ceil(n_cells / n_cols))


plt_idc = [] 

# construct indices array for all plots in subplots
for row in np.arange(n_rows):
    for col in np.arange(n_cols):
        plt_idc.append((row, col))


# initialise figure
fig_rest, axs_rest = plt.subplots(n_rows, n_cols, 
                                  sharex=True, sharey=True)
                                  #figsize = get_figure_size())


# convert indices of dataframe to list to loop through
all_cell_IDs = lookup_table.index.to_list()

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    
    # get indices of current cell with the dataframe containing all indices    
    group_idx = int(lookup_table.at[cell_ID, 'group_idx'])
    series_idx = int(lookup_table.at[cell_ID, 'series_idx'])

    # construct traceIndex with indices
    traceIndex = [group_idx, series_idx, 0, 0]
    
    # call on data file with indices from dataframe above
    current_file = lookup_table.at[cell_ID, 'file']
    
    data_file_path = os.path.join(data_folder, current_file + '.dat')

    data_file_path_str = fr"{data_file_path}"
    
    
    # get data with file path & trace index
    i, v, t, SR = get_data(data_file_path_str, traceIndex, scale='s')
    
    v_filtered = butter_filter(v, order=3, cutoff=1e3, sampling_rate=SR)
    
    axs_rest[plt_idc[cell_idx][0]][plt_idc[cell_idx][1]].plot(t, v_filtered)
    
    

axs_rest[-1][-1].set_xlim([t[0], t[-1]])
axs_rest[-1][-1].set_ylim([-100, 40])

fig_rest.supylabel('Membrane potential [mV]')
fig_rest.supxlabel('Time [s]')

