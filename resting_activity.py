#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:02:40 2023

@author: moritznesseler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os.path
import scipy as sc
import seaborn as sbn
import warnings 

# custom directories & parameters
from directories_win import raw_data_dir, figure_dir
from parameters import min_peak_prominence, min_peak_distance

# custom functions
from functions_useful import butter_filter, get_data
from functions_plotting import get_figure_size, get_colors, save_figures, set_font_sizes
from functions_export import set_df_to_cell_descrips



# %%

## n/3 by 3 subplots with each a 30 sec resting activity plot

table = pd.read_excel('//Fileserver/AG Spehr/File transfer/Moritz_transfer/InVitro_Database.xlsx',
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

spiketimes_ls = []
n_peaks = []

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    idx_peaks, dict_peak = sc.signal.find_peaks(v_filtered_df[cell_idx], 
                                                prominence = min_peak_prominence, 
                                                distance = min_peak_distance * (SR/1e3))

    t_peaks = idx_peaks / (SR_ls[cell_idx])
    
    spiketimes_ls.append(t_peaks)
    
    n_peaks.append(len(idx_peaks))
 
    
# create dataframe for browsable data
n_spikes_df = pd.DataFrame({'n_spikes' : n_peaks},
                            index = all_cell_IDs)

# %%

# get colors for plotting
colors_dict = get_colors(darkmode_bool)

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
activity_df = pd.concat([v_rest_df, n_spikes_df], axis = 1)

# add activity column for categorical plots´with silent as default value
activity_df['activity'] = 'silent'

# change activity value of spiking cells with n_spike > 0 to 'spiking'
activity_df.loc[activity_df['n_spikes'] > 0, 'activity'] = 'spiking'



# save activity dataframe to quant data folder
# activity_df.to_excel(cell_descrip_file, index_label='cell_ID')

set_df_to_cell_descrips(activity_df)


# %% EVENTPLOT + V_REST FIGURE + N_spike

fig_v_rest, axs_v_rest = plt.subplots(1,2,
                                      gridspec_kw={'width_ratios': [3,1]},
                                      figsize = get_figure_size(),
                                      layout = 'tight')

set_font_sizes()

tick_size = 0.9

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    
    axs_v_rest[0].eventplot(spiketimes_ls[cell_idx],
                            orientation = 'horizontal', 
                            lineoffsets=cell_idx, 
                            linewidth = 1.5,
                            linelengths=0.9, 
                            color = colors_dict['color2'])


axs_v_rest[0].set_ylim([0-(tick_size/2), n_cells-1+(tick_size/2)])
axs_v_rest[0].set_yticks(ticks = np.arange(5 - 1, n_cells+1, 5), 
                         labels = np.arange(5, n_cells + 1, 5))
axs_v_rest[0].set_yticks(ticks = np.arange(0, n_cells, 1), 
                         minor = True)
axs_v_rest[0].set_ylabel('Cells [#]')

axs_v_rest[0].set_xlim([0, 30])
axs_v_rest[0].set_xlabel('Time [s]')
axs_v_rest[0].set_xticks(np.arange(0, 30+1, 10))
axs_v_rest[0].set_xticks(np.arange(0, 30+1, 1), minor = True)

# creating specific color pallet for seaborn plotting functions
color_pal = {'silent' : colors_dict['color1'], 'spiking' : colors_dict['color2']}

violins = sbn.violinplot(data = activity_df, 
                         x = "activity", 
                         y = "v_rest", 
                         inner = 'quart', 
                         ax = axs_v_rest[1], 
                         palette = color_pal,
                         linewidth = 1.5)

swarms = sbn.swarmplot(data = activity_df, 
                       x = "activity", 
                       y = "v_rest", 
                       size = 7, 
                       ax = axs_v_rest[1], 
                       color = colors_dict['primecolor'])

for l in violins.lines:
    l.set_color(colors_dict['primecolor'])

[violin.set_edgecolor(colors_dict['primecolor']) for violin in violins.collections]

axs_v_rest[1].set_ylabel('Resting membrane potential [mV]')
axs_v_rest[1].set_ylim([-100, -40])
axs_v_rest[1].set_yticks(np.arange(-100, -40+1, 10))
axs_v_rest[1].set_yticks(np.arange(-100, -40+1, 5), minor = True)

axs_v_rest[1].set_xlabel('')

[ax.grid(False) for ax in axs_v_rest]

save_figures(fig_v_rest, 'Resting_n_eventplot_nspikes', figure_dir, darkmode_bool)





