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
from tqdm import tqdm

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir
from parameters.parameters import min_peak_prominence, min_peak_distance

# custom functions
from functions.functions_useful import butter_filter, calc_time_series, calc_dvdt_padded
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol


# define protocol
PGF = 'cc_rest'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = 'PGFs_Syn')

# get number of cells
n_cells = len(cell_IDs)

# init plotting
from functions.initialize_plotting import * # analysis:ignore

# define output
activity_df = pd.DataFrame(index = cell_IDs, columns = ['v_rest', 'n_spikes', 't_spikes'])

# verification plots
vplots = True


# %% data loading

# initialize dataframes to populate in the loop
v_df = pd.DataFrame(columns = cell_IDs)
vf_df = pd.DataFrame(columns = cell_IDs)
SR_df = pd.DataFrame(index = cell_IDs, columns = ['SR'])


print('loading ...')

for cell_idx, cell_ID in enumerate(tqdm(cell_IDs)):
            
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name = 'PGFs_Syn')
    
    # get data with file path & trace index
    i, v, t, SR, n_step = get_cc_data(file_path, traceIndex, scale='s')
    
    # get first and only step of protocol
    v = v[0]
    
    # edge case when file exceeds the 30 sec recording (E-069)
    if len(v) > (30 * SR):
        warnings.warn(str(cell_ID) + ' exceeds 30 sec and will be cut down.')
        v = v[0:(30*SR)]
        
    # filter all data with 1kHz cutoff
    vf = butter_filter(v, order=3, cutoff=1e3, sampling_rate=SR)

    # populate the dataframes & lists
    v_df[cell_ID] = v
    vf_df[cell_ID] = vf
    SR_df.at[cell_ID, 'SR'] = SR
    

# check if all protocols have the same sampling rate
if len(SR_df['SR'].unique()) != 1:
    warnings.warn('Not all protocols have the same sampling rate!')

else:
    # calculate time
    t = calc_time_series(v, sampling_rate=SR, scale = 's')
    t_ms = calc_time_series(v, sampling_rate=SR, scale = 'ms')


# %% find spikes

cell_ID = 'E-247'

# define dataframe
spiketimes_df = pd.DataFrame(index = cell_IDs, columns=['t_spikes'])

# get filtered voltage and sampling rate
vf = vf_df[cell_ID] #[25:]
SR = SR_df.at[cell_ID, 'SR']

# find peaks
idc_spikes, dict_peak = sc.signal.find_peaks(vf, 
                                             prominence = min_peak_prominence, 
                                             distance = min_peak_distance * (SR/1e3))
    
# calculate spike times in seconds
t_spikes = np.divide(idc_spikes, SR)
n_spikes = len(t_spikes)

# write to dataframe
activity_df.at[cell_ID, 't_spikes'] = list(t_spikes)
activity_df.at[cell_ID, 'n_spikes'] = n_spikes

plt.plot(vf[idc_spikes[0]-2000:idc_spikes[0]+2000])


# %% v_rest

from functions.functions_extractspike import extract_spike
from functions.functions_spiketrains import calc_vmem_at_spiketrain

# cut out spikes for resting membrane potential calc
if n_spikes > 0:
    
    # calc first derivative
    dvdt = calc_dvdt_padded(vf, t_ms)
    
    extract_spike(t = t_ms, 
                  v = vf, 
                  dvdt = dvdt, 
                  idx_peak = idc_spikes[0])


    
# %%







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
    
    print(f'Started: {cell_ID}')
        
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
# v_rest_path = os.path.join(figure_dir, 'v_rest.csv')

# v_rest_df.to_csv(v_rest_path, header = ['v_rest'])



# %% V_REST mean & std

# v_rest_mean = np.mean(v_rest_df['v_rest'])
# v_rest_std = np.std(v_rest_df['v_rest'])



# %% create dict with all, active and non-active cells

# concatenate v_rest and n_spikes dataframes
# activity_df = pd.concat([v_rest_df, n_spikes_df, spiketimes_df], axis = 1)

# add activity column for categorical plots´with silent as default value
activity_df['activity'] = 'silent'

# change activity value of spiking cells with n_spike > 0 to 'spiking'
activity_df.loc[activity_df['n_spikes'] > 0, 'activity'] = 'spiking'



# save activity dataframe to quant data folder
activity_df.to_excel(os.path.join(cell_descrip_dir, 'cc_rest-syn-activity.xlsx'), index_label='cell_ID')



print('Finished!')





