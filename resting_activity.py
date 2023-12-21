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
import warnings 


from plotting_functions import get_figure_size, get_colors


# %%

## n/3 by 3 subplots with each a 30 sec resting activity plot




# lookup_table = pd.read_csv('/Users/moritznesseler/local E-Phys/cc_rest.csv', 
#                            delimiter=';',
#                            index_col='cell_ID')
# data_folder = '/Users/moritznesseler/local E-Phys'
# figure_dir = '/Users/moritznesseler/local E-Phys/figures'


## excel file loading routine for group & series indices


table = pd.read_excel('//Fileserver/AG Spehr/File transfer/Moritz_transfer/InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')

lookup_table = table.query('cc_rest.notnull()')

data_folder = 'C:/Users/nesseler/Desktop/local E-Phys'

figure_dir = 'C:/Users/nesseler/Desktop/local E-Phys/figures'


darkmode_bool = False

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
    
    data_file_path = os.path.join(data_folder, current_file + '.dat')

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
    
   
# %% SUBPLOTS RESTING ACTIVITY

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
                                  sharex=True, sharey=True,
                                  figsize = get_figure_size(),
                                  layout = 'tight')

small_font_size = 14
large_font_size = 16

plt.rc('font', size = small_font_size)
plt.rc('axes', titlesize = large_font_size, 
               labelsize = large_font_size,
               linewidth = 0.5)
plt.rc('xtick', labelsize = large_font_size)
plt.rc('ytick', labelsize = large_font_size)
plt.rc('lines', linewidth = 1)

plots_alone = n_cells % 3


for cell_idx, cell_ID in enumerate(all_cell_IDs):
    axs_rest[plt_idc[cell_idx][0]][plt_idc[cell_idx][1]].plot(t_df[cell_idx], 
                                                              v_filtered_df[cell_idx], 
                                                              color=colors_dict['primecolor'])
    
   
axs_rest[-1][-1].set_xlim([t[0], t[-1]])
axs_rest[-1][-1].set_xticks(np.arange(0, 30+1, 10))
axs_rest[-1][-1].set_xticks(np.arange(0, 30+1, 5), minor = True)

axs_rest[-1][-1].set_ylim([-100, 50])
axs_rest[-1][-1].set_yticks([])
axs_rest[-1][-1].set_yticks(np.arange(-100, 50 + 1, 25), minor = True)


if plots_alone != 0:
    for ax in axs_rest[-1][plots_alone:]:
        ax.remove()


fig_rest.supylabel('Membrane potential [mV]')
fig_rest.supxlabel('Time [s]')

[ax.grid(False) for rows in axs_rest for ax in rows]
    
save_figures(fig_rest, 'resting_subplots', figure_dir, darkmode_bool)


# %% EVENTPLOT



## eventplot with all cells

min_peak_prominence = 40 #(mV)
min_peak_distance = 1 #ms

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

fig_rest_event, ax_rest_event = plt.subplots(1,1)
                                            #figsize = get_figure_size())

tick_size = 0.9

spiketimes_ls.sort(key=len)

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    
    ax_rest_event.eventplot(spiketimes_ls[cell_idx],
                            orientation = 'horizontal', 
                            lineoffsets=cell_idx, 
                            linelengths=tick_size, 
                            linewidths = 0.75,
                            color = colors_dict['color2'])


ax_rest_event.set_xlim([0, 30])
fig_rest_event.supxlabel('Time [s]')

ax_rest_event.set_ylim([0-(tick_size/2), n_cells+(tick_size/2)])
ax_rest_event.set_yticks(np.arange(0,n_cells+1,10))
ax_rest_event.set_yticks(np.arange(0,n_cells+1), minor=True)
fig_rest_event.supylabel('Cells [#]')

plt.grid(False)

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
    t_pre = 5
    t_post = 40
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

# create dataframe for browsable data
v_rest_df = pd.DataFrame({'v_rest' : v_rest},
                         index = all_cell_IDs)

# save dataframe as csv file
v_rest_path = os.path.join(figure_dir, 'v_rest.csv')

v_rest_df.to_csv(v_rest_path, header = ['v_rest'])


# %% V_REST mean & std

v_rest_mean = np.mean(v_rest_df['v_rest'])
v_rest_std = np.std(v_rest_df['v_rest'])

# %% EVENTPLOT + V_REST FIGURE



fig_v_rest, axs_v_rest = plt.subplots(1,2,
                                      gridspec_kw={'width_ratios': [4,1]},
                                      figsize = get_figure_size(),
                                      layout = 'tight')

small_font_size = 14
large_font_size = 16

plt.rc('font', size = small_font_size)
plt.rc('axes', titlesize = large_font_size, 
               labelsize = large_font_size,
               linewidth = 0.5)
plt.rc('xtick', labelsize = large_font_size)
plt.rc('ytick', labelsize = large_font_size)
plt.rc('lines', linewidth = 2)


tick_size = 0.9

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    
    axs_v_rest[0].eventplot(spiketimes_ls[cell_idx],
                            orientation = 'horizontal', 
                            lineoffsets=cell_idx, 
                            linelengths=0.9, 
                            color = colors_dict['color2'])


# ax_rest_event.set_xlim([0, 30])
# fig_rest_event.supxlabel('Time [s]')


axs_v_rest[0].set_ylim([0-(tick_size/2), n_cells+(tick_size/2)])
axs_v_rest[0].set_yticks(np.arange(0,n_cells+1,3))
axs_v_rest[0].set_yticks(np.arange(0,n_cells+1,1), minor=True)
axs_v_rest[0].set_ylabel('Cells [#]')

axs_v_rest[0].set_xlim([0, 30])
axs_v_rest[0].set_xlabel('Time [s]')
axs_v_rest[0].set_xticks(np.arange(0, 30+1, 10))
axs_v_rest[0].set_xticks(np.arange(0, 30+1, 1), minor = True)


# axs_v_rest.scatter([1]*len(v_rest), v_rest)
# axs_v_rest[1].violinplot(v_rest_df['v_rest'],
#                       positions = [0])#,
#                       # facecolor = 'None',
#                       # edgecolor = 'w')

# sbn.violinplot(x=['v_rest']*19,
#                y=v_rest_df['v_rest'], ax=axs_v_rest[1],
#                fill = False,
#                color='w')

sbn.swarmplot(x=[0]*len(v_rest_df),
              y=v_rest_df['v_rest'], 
              ax=axs_v_rest[1],
              color = colors_dict['primecolor'], 
              size = 7)


axs_v_rest[1].errorbar(x = 0, 
                       y = v_rest_mean,
                       yerr = v_rest_std,
                       marker = '_',
                       markersize = 20,
                       color = 'r',
                       capsize = 5,
                       capthick = 3,
                       ecolor = 'r',
                       linewidth = 3)


axs_v_rest[1].set_ylabel('Resting membrane potential [mV]')
axs_v_rest[1].set_ylim([-100, -40])
axs_v_rest[1].set_yticks(np.arange(-100, -40+1, 10))
axs_v_rest[1].set_yticks(np.arange(-100, -40+1, 5), minor = True)

axs_v_rest[1].set_xlim([-1, 1])
axs_v_rest[1].set_xlabel('')
axs_v_rest[1].set_xticks([])

[ax.grid(False) for ax in axs_v_rest]

save_figures(fig_v_rest, 'Resting_n_eventplot', figure_dir, darkmode_bool)




# %% EVENTPLOT + V_REST FIGURE + N_spike

fig_v_rest, axs_v_rest = plt.subplots(1,2,
                                      gridspec_kw={'width_ratios': [3,1]},
                                      figsize = get_figure_size(),
                                      layout = 'tight')

small_font_size = 14
large_font_size = 16

plt.rc('font', size = small_font_size)
plt.rc('axes', titlesize = large_font_size, 
               labelsize = large_font_size,
               linewidth = 0.5)
plt.rc('xtick', labelsize = large_font_size)
plt.rc('ytick', labelsize = large_font_size)
plt.rc('lines', linewidth = 1.5)


tick_size = 0.9

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    
    axs_v_rest[0].eventplot(spiketimes_ls[cell_idx],
                            orientation = 'horizontal', 
                            lineoffsets=cell_idx, 
                            linewidth = 1,
                            linelengths=0.9, 
                            color = colors_dict['color2'])


axs_v_rest[0].set_ylim([0-(tick_size/2), n_cells-1+(tick_size/2)])
axs_v_rest[0].set_yticks(np.arange(0,n_cells,5))
axs_v_rest[0].set_yticks(np.arange(0,n_cells,1), minor=True)
axs_v_rest[0].set_ylabel('Cells [#]')

axs_v_rest[0].set_xlim([0, 30])
axs_v_rest[0].set_xlabel('Time [s]')
axs_v_rest[0].set_xticks(np.arange(0, 30+1, 10))
axs_v_rest[0].set_xticks(np.arange(0, 30+1, 1), minor = True)


sbn.swarmplot(x=[0]*len(v_rest_df),
              y=v_rest_df['v_rest'], 
              ax=axs_v_rest[1],
              color = colors_dict['primecolor'], 
              size = 7)


spiking_indices = n_spikes_df.query('n_spikes > 0').index
silent_indices = n_spikes_df.query('n_spikes == 0').index

sbn.swarmplot(x=[1]*len(spiking_indices),
              y=v_rest_df['v_rest'].loc[spiking_indices], 
              ax=axs_v_rest[1],
              color = colors_dict['color2'], 
              size = 7)

sbn.swarmplot(x=[2]*len(silent_indices),
              y=v_rest_df['v_rest'].loc[silent_indices], 
              ax=axs_v_rest[1],
              color = colors_dict['color1'], 
              size = 7)


axs_v_rest[1].set_ylabel('Resting membrane potential [mV]')
axs_v_rest[1].set_ylim([-100, -40])
axs_v_rest[1].set_yticks(np.arange(-100, -40+1, 10))
axs_v_rest[1].set_yticks(np.arange(-100, -40+1, 5), minor = True)

axs_v_rest[1].set_xlim([-1, 3])
axs_v_rest[1].set_xlabel('')
axs_v_rest[1].set_xticks([0, 1, 2], ['all', 'active', 'non\nactive'],
                         rotation = 90)

[ax.grid(False) for ax in axs_v_rest]

save_figures(fig_v_rest, 'Resting_n_eventplot_nspikes', figure_dir, darkmode_bool)






