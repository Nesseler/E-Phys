# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 18:14:12 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import numpy as np
import directories_win as directories
from PGFs import cc_APs_parameters
import pandas as pd
from useful_functions import calc_time_series, butter_filter
import os
from cc_IF_functions import get_IF_data
from plotting_functions import get_colors, save_figures


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

dm_bool = False

colors_dict = get_colors(dm_bool)


fig_spiketrains, axs_spiketrains = plt.subplots(3, 3,
                                                layout = 'constrained',
                                                sharey = 'col',
                                                dpi = 600)


axs_spiketrains[0][1].plot(t_dict['1Hz'], vf_dict['1Hz'])




plt.show()























