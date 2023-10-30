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
import scipy as sc
import seaborn as sbn


from plotting_functions import get_figure_size, get_colors


# %%

## n/3 by 3 subplots with each a 30 sec resting activity plot




lookup_table = pd.read_csv('/Users/moritznesseler/local E-Phys/cc_rest.csv', 
                           delimiter=';',
                           index_col='cell_ID')


data_folder = '/Users/moritznesseler/local E-Phys'

figure_dir = '/Users/moritznesseler/local E-Phys'

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
    
    i_df[cell_idx] = i
    v_df[cell_idx] = v
    v_filtered_df[cell_idx] = v_filtered
    t_df[cell_idx] = t
    SR_ls.append(SR)
    
   
# %%

n_cells = len(lookup_table)
n_cols = 3
n_rows = int(np.ceil(n_cells / n_cols))


plt_idc = [] 

# construct indices array for all plots in subplots
for row in np.arange(n_rows):
    for col in np.arange(n_cols):
        plt_idc.append((row, col))
        
        
# get colors for plotting
colors_dict = get_colors(darkmode_bool)

# initialise figure
fig_rest, axs_rest = plt.subplots(n_rows, n_cols, 
                                  sharex=True, sharey=True)
                                  #figsize = get_figure_size())


for cell_idx, cell_ID in enumerate(all_cell_IDs):
    axs_rest[plt_idc[cell_idx][0]][plt_idc[cell_idx][1]].plot(t_df[cell_idx], v_filtered_df[cell_idx], color=colors_dict['prime_color'])
    
   

axs_rest[-1][-1].set_xlim([t[0], t[-1]])
axs_rest[-1][-1].set_ylim([-100, 40])

fig_rest.supylabel('Membrane potential [mV]')
fig_rest.supxlabel('Time [s]')

# ToDo
    # tick marks
    # darkmode
    # save figure
    
save_figures(fig_rest, 'resting_subplots', '/Users/moritznesseler/Desktop/local E-Phys', darkmode_bool)


# %% EVENTPLOT

## eventplot with all cells

min_peak_prominence = 40 #(mV)
min_peak_distance = 1 #ms

spiketimes_ls = []

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    idx_peaks, dict_peak = sc.signal.find_peaks(v_filtered_df[cell_idx], prominence = min_peak_prominence, distance = min_peak_distance * (SR/1e3))

    t_peaks = idx_peaks / (SR_ls[cell_idx])
    
    spiketimes_ls.append(t_peaks)


fig_rest_event, ax_rest_event = plt.subplots(1,1)
                                            #figsize = get_figure_size())

tick_size = 0.9

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    
    ax_rest_event.eventplot(spiketimes_ls[cell_idx],
                            orientation = 'horizontal', 
                            lineoffsets=cell_idx, 
                            linelengths=0.9, 
                            color = colors_dict['color2'])


ax_rest_event.set_xlim([0, 30])
fig_rest_event.supxlabel('Time [s]')


ax_rest_event.set_ylim([0-(tick_size/2), n_cells+(tick_size/2)])
ax_rest_event.set_yticks(np.arange(0,n_cells+1,3))
ax_rest_event.set_yticks(np.arange(0,n_cells+1), minor=True)
fig_rest_event.supylabel('Cells [#]')

save_figures(fig_rest_event, 'Resting_eventplot', figure_dir, darkmode_bool)



# %% V_REST

# TODO 
 # try baseline function & filtered
 # edge cases, when spike is near start or stop for nan replacement


v_rest = []

# initialise figure
# fig_rest, axs_rest = plt.subplots(n_rows, n_cols, 
#                                   sharex=True, sharey=True)
#                                   #figsize = get_figure_size())


for cell_idx, cell_ID in enumerate(all_cell_IDs):
        
    v_cell = v_filtered_df[cell_idx].to_list()
    
    idc_peaks, dict_peak = sc.signal.find_peaks(v_cell, 
                                                prominence = min_peak_prominence, 
                                                distance = min_peak_distance * (SR/1e3))

    # introduce NaN value at spike times for more accurate v_rest calculation
    
    # time before and after peak to exclude (in ms)
    t_pre = 2
    t_post = 5
    t_total = t_pre + t_post
    
    for i, idx_peak in enumerate(idc_peaks):
        idx_pre = idx_peak - int((t_pre * (SR_ls[cell_idx]/1e3)))
        idx_post = idx_peak + int((t_post * (SR_ls[cell_idx]/1e3)))
        
        v_cell[idx_pre:idx_post] = [np.nan] * int((t_total * (SR_ls[cell_idx] / 1e3)))
    
    row_idx = plt_idc[cell_idx][0]
    col_idx = plt_idc[cell_idx][1]
    
    # axs_rest[row_idx][col_idx].plot(t_df[cell_idx], 
    #                                 v_cell, 
    #                                 color=colors_dict['prime_color'])
    
    # axs_rest[row_idx][col_idx].hlines(np.mean(v_filtered_df[cell_idx]), 
    #                                   0, 
    #                                   30, 
    #                                   linestyle='--', color='r')
    
    v_rest.append(np.nanmean(v_cell))

# axs_rest[-1][-1].set_xlim([22.63, 22.64])
# axs_rest[-1][-1].set_ylim([-100, 40])

# fig_rest.supylabel('Membrane potential [mV]')
# fig_rest.supxlabel('Time [s]')

# %% V_REST FIGURE

fig_v_rest, axs_v_rest = plt.subplots(1,1)


# axs_v_rest.scatter([1]*len(v_rest), v_rest)


sbn.swarmplot(v_rest, ax=axs_v_rest)

axs_v_rest.set_ylim([-100, -40])



