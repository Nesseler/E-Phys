# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 14:26:47 2023

@author: nesseler
"""

# E-092: use for testing

import directories_win as directories
import pandas as pd
import os
from cc_IF_functions import get_IF_data
import matplotlib.pyplot as plt
import numpy as np
from useful_functions import calc_time_series, butter_filter
import parameters
from PGFs import cc_APs_parameters
import scipy as sc

table = pd.read_excel(directories.table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')

# PGF to load
frequency = '30Hz'
PGF = 'cc_APs_' + frequency
PGF_parameters = cc_APs_parameters[frequency]

lookup_table = table.query(f'{PGF}.notnull()')


# test cell
cell_ID = 'E-092'


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

# plt.plot(t_concat, v_concat)
# plt.show()


# limit the array to stimulus time frame + set time post stimulus to accomodate late APs
t_post_stim = parameters.cc_APs_t_post_stim

idc_stim = np.arange(start = PGF_parameters['t_pre']*SR_ms,
                     stop = PGF_parameters['t_pre']*SR_ms + PGF_parameters['t_stim']*SR_ms + t_post_stim*SR_ms,
                     dtype=int)

# limit voltage array
v_stim = v[:, idc_stim]


# initialize lists
t_spikes = []

for step_idx in np.arange(n_steps):
    
    # filter voltage (to vf)
    vf = butter_filter(v[step_idx], 
                      order = 3,
                      cutoff = 1e3,
                      sampling_rate = SR)
    
    # limit data array to just the stimulus
    vl = vf[idc_stim]
    
    # find peaks as spikes
    idx_peaks, dict_peak = sc.signal.find_peaks(vl, 
                                                prominence = parameters.min_peak_prominence, 
                                                distance = parameters.min_peak_distance * (SR_ms))
    
    # get times of spikes
    t_peaks = np.divide(idx_peaks, (SR_ms))
    t_spikes.append(t_peaks)

    ### opt: PLOT ###
    plt.plot(vl)
    plt.eventplot(idx_peaks, color = 'r', lineoffsets=30, linelengths=5)
    plt.ylim([-100, 50])
    plt.show()
    

