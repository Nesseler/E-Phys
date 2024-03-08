#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:02:40 2023

@author: moritznesseler
"""

import pandas as pd
import numpy as np
import os.path
import scipy as sc
import warnings 

# custom directories & parameters
from parameters.directories_win import raw_data_dir, figure_dir, cell_descrip_dir, table_file
from parameters.parameters import min_peak_prominence, min_peak_distance

# custom functions
from functions.functions_useful import butter_filter, get_data




# %%

## n/3 by 3 subplots with each a 30 sec resting activity plot

table = pd.read_excel(table_file,
                      sheet_name="PGFs",
                      index_col='cell_ID')

lookup_table = table.query('cc_rest.notnull()')

darkmode_bool = True

# %%

# initialize dataframes to populate in the loop
i_df = pd.DataFrame()
v_df = pd.DataFrame()
v_filtered_df = pd.DataFrame()
t_df = pd.DataFrame()
SR_ls = []



# convert indices of dataframe to list to loop through
all_cell_IDs = lookup_table.index.to_list()

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    
    # get indices of current cell with the dataframe containing all indices    
    group_idx = int(lookup_table.at[cell_ID, 'group'])-1
    series_idx = int(lookup_table.at[cell_ID, 'cc_rest'])-1

    # construct traceIndex with indices
    traceIndex = [group_idx, series_idx, 0, 0]
    
    # call on data file with indices from dataframe above
    current_file = lookup_table.at[cell_ID, 'file']
    
    data_file_path = os.path.join(raw_data_dir, current_file + '.dat')

    data_file_path_str = fr"{data_file_path}"
    
    
    # get data with file path & trace index
    i, v, t, SR = get_data(data_file_path_str, traceIndex, scale='s')
    
    # edge case when file exceeds the 30 sec recording (E-069)
    if len(v) > (30 * SR):
        warnings.warn(str(cell_ID) + ' exceeds 30 sec and will be cut down.')
        i = i[0:(30*SR)]
        v = v[0:(30*SR)]
        t = t[0:(30*SR)]
        
    # filter all data with 1kHz cutoff
    v_filtered = butter_filter(v, order=3, cutoff=1e3, sampling_rate=SR)

    # populate the dataframes & lists
    i_df[cell_idx] = i
    v_df[cell_idx] = v
    v_filtered_df[cell_idx] = v_filtered
    t_df[cell_idx] = t
    SR_ls.append(SR)
    

n_cells = len(lookup_table)



# %% find peaks

spiketimes_df = pd.DataFrame(columns=['t_spikes'])
spiketimes_ls = []
n_peaks = []

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    idx_peaks, dict_peak = sc.signal.find_peaks(v_filtered_df[cell_idx], 
                                                prominence = min_peak_prominence, 
                                                distance = min_peak_distance * (SR/1e3))

    t_peaks = idx_peaks / (SR_ls[cell_idx])
    
    spiketimes_ls.append(t_peaks)
    spiketimes_df.at[cell_ID, 't_spikes'] = list(t_peaks)
    
    n_peaks.append(len(idx_peaks))
 
    
# create dataframe for browsable data
n_spikes_df = pd.DataFrame({'n_spikes' : n_peaks},
                            index = all_cell_IDs)

# %%



# sort spike times list for length and therefore number of APs
spiketimes_ls.sort(key=len)


# %% V_REST

# TODO 
 # try baseline function & filtered
 # edge cases, when spike is near start or stop for nan replacement


v_rest = []


for cell_idx, cell_ID in enumerate(all_cell_IDs):
        
    v_cell = v_filtered_df[cell_idx].to_list()
    
    idc_peaks, dict_peak = sc.signal.find_peaks(v_cell, 
                                                prominence = min_peak_prominence, 
                                                distance = min_peak_distance * (SR/1e3))

    # introduce NaN value at spike times for more accurate v_rest calculation
    
    # time before and after peak to exclude (in ms)
    t_pre = 5
    t_post = 40
    t_total = t_pre + t_post
    
    for i, idx_peak in enumerate(idc_peaks):
        idx_pre = idx_peak - int((t_pre * (SR_ls[cell_idx]/1e3)))
        idx_post = idx_peak + int((t_post * (SR_ls[cell_idx]/1e3)))
        
        v_cell[idx_pre:idx_post] = [np.nan] * int((t_total * (SR_ls[cell_idx] / 1e3)))
       
    v_rest.append(np.nanmean(v_cell))

# create dataframe for browsable data
v_rest_df = pd.DataFrame({'v_rest' : v_rest},
                         index = all_cell_IDs)

# save dataframe as csv file
v_rest_path = os.path.join(figure_dir, 'v_rest.csv')

v_rest_df.to_csv(v_rest_path, header = ['v_rest'])



# %% V_REST mean & std

v_rest_mean = np.mean(v_rest_df['v_rest'])
v_rest_std = np.std(v_rest_df['v_rest'])



# %% create dict with all, active and non-active cells

# concatenate v_rest and n_spikes dataframes
activity_df = pd.concat([v_rest_df, n_spikes_df, spiketimes_df], axis = 1)

# add activity column for categorical plots´with silent as default value
activity_df['activity'] = 'silent'

# change activity value of spiking cells with n_spike > 0 to 'spiking'
activity_df.loc[activity_df['n_spikes'] > 0, 'activity'] = 'spiking'



# save activity dataframe to quant data folder
activity_df.to_excel(os.path.join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_label='cell_ID')









