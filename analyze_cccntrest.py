# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:41:30 2024

@author: nesseler
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

# custom directories & parameters
from parameters.directories_win import table_dir
from parameters.parameters import min_peak_distance, min_peak_prominence, min_spike_in_burst

# custom functions
from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_useful import calc_time_series, butter_filter
from functions.functions_plotting import get_colors, get_figure_size


table = pd.read_excel(table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')


# loop to create string to include all frequencies in query
PGF = 'cc_cnt_rest'  

# limit lookup table
lookup_table = table.query(f'{PGF}.notnull()')

# cell IDs 
cell_ID = 'E-113'


# %%


    
print(f'Started: {cell_ID}')

# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)

# get IF data form file
i, v, t, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')

# sampling rate in ms
SR_ms = SR / 1e3

# concatenate individual steps
n_points = int(np.shape(i)[0] * np.shape(i)[1])

i_concat = i.flatten() 

v_concat = v.flatten() 

t_ms = calc_time_series(v_concat, SR)
t_s = calc_time_series(v_concat, SR, scale = 's')

# plt.plot(v_concat[0:15000])
# plt.show()


# filter voltage (to vf)
vf = butter_filter(v_concat, 
                   order = 3,
                   cutoff = 1e3,
                   sampling_rate = SR)



# find peaks
idx_peaks, dict_peak = sc.signal.find_peaks(vf, 
                                            prominence = min_peak_prominence, 
                                            distance = min_peak_distance * (SR/1e3))

t_peaks_s = np.divide(idx_peaks, SR)

spikes_df = pd.DataFrame({'t_spikes' : t_peaks_s})

# %% calc ISIs

# calculate ISIs in seconds
ISIs_s = np.diff(t_peaks_s)

# pad array with nan for first spike
ISIs_s = np.pad(ISIs_s,
                pad_width = (1, 0),
                mode = 'constant',
                constant_values = (np.nan,))

# write to spikes_df
spikes_df['ISIs'] = ISIs_s



mean_ISI = spikes_df['ISIs'].mean()
median_ISI = spikes_df['ISIs'].median()

mean_FR = 1 / mean_ISI


pos_ISI = 600 / len(spikes_df.index)


# %% poisson distribution

n_spikes = len(t_peaks_s)

uni_t_spikes = np.random.uniform(0, 600, n_spikes)
uni_t_spikes = np.sort(uni_t_spikes)
uni_ISIs = np.diff(uni_t_spikes)
uni_mean_ISI = np.mean(uni_ISIs)
uni_std_ISI = np.std(uni_ISIs)

x = np.arange(0, 20, 0.1)
x1 = np.arange(0, 20, 0.1)

poi_ISI_pmf = sc.stats.norm.sf(x = x, loc = uni_mean_ISI, scale = uni_std_ISI)

poi_ISI_cdf = sc.stats.norm.cdf(x = x, loc = uni_mean_ISI, scale = uni_std_ISI)



# plt.hist(poi_ISI_dist, density=True, edgecolor='black')

plt.plot(x, poi_ISI_pmf, 'r')
plt.plot(x, poi_ISI_cdf, 'gray')

bins = np.arange(0, 20, 0.1)

plt.hist(ISIs_s, bins, density=True)

plt.xlim([-1, 5])
         
plt.show()

# %% burst or not 

# set non burst condition
spikes_df['burst'] = [0] * len(spikes_df.index)
spikes_df['burst_id'] = [np.nan] * len(spikes_df.index)


burst_id = 0
next_burst = False



# moving window approach to burst classification
for idx, t_spike in enumerate(t_peaks_s):
    
    if idx >= min_spike_in_burst:
        ISI_mov_window = ISIs_s[idx-3:idx+1]
        
        if all([ISI <= uni_mean_ISI for ISI in ISI_mov_window]):
            spikes_df.loc[idx-min_spike_in_burst:idx, 'burst'] = 1
            spikes_df.loc[idx-min_spike_in_burst:idx, 'burst_id'] = burst_id
            
            next_burst = True          
        else:
            if next_burst:
                burst_id += 1
                next_burst = False



# %%

colors_dict = get_colors(False)


fig_trace, ax_trace = plt.subplots(nrows = 3,
                                   ncols = 1,
                                   layout = 'constrained',
                                   sharex = 'col',
                                   sharey = 'row',
                                   gridspec_kw={'height_ratios': [1,1,1]},
                                   figsize = get_figure_size())




ax_trace[0].plot(t_s, vf)

# ax_trace[0].set_xlim([430,470])

ax_trace[0].eventplot(spikes_df.loc[spikes_df['burst'] == 1,'t_spikes'],
                      orientation = 'horizontal', 
                      lineoffsets=60, 
                      linewidth = .5,
                      linelengths=10, 
                      color = [colors_dict['color2']])

ax_trace[0].eventplot(spikes_df.loc[spikes_df['burst'] == 0,'t_spikes'],
                      orientation = 'horizontal', 
                      lineoffsets=60, 
                      linewidth = .5,
                      linelengths=10, 
                      color = [colors_dict['color1']])


# plot ISIs

ax_trace[1].plot(t_peaks_s, ISIs_s, '-o')



# %% histogram


fig_hist, ax_hist = plt.subplots(nrows = 2,
                                 ncols = 2,
                                 layout = 'constrained',
                                 sharex = 'col',
                                 sharey = 'row',
                                 gridspec_kw={'height_ratios': [1, 3], 'width_ratios': [2, 1]},
                                 figsize = get_figure_size())

fig_hist.suptitle(cell_ID)

bin_width = .5 # mV
bins_min = -100
bins_max = 40

ylim_border = 0.1

bins = np.arange(bins_min, bins_max + bin_width, bin_width)

print(bins)

N, bins, patches = ax_hist[1][0].hist(vf, bins = bins, density = True)


ax_hist[1][0].set_ylim([0, ylim_border])
ax_hist[1][0].set_xticks(np.arange(bins_min, bins_max + bin_width, 20))
ax_hist[1][0].set_xticks(bins, minor = True)



N, bins, patches = ax_hist[0][0].hist(vf, bins = bins, density = True)


ax_hist[0][0].set_ylim([ylim_border, 1])


N, bins, patches = ax_hist[0][1].hist(vf, bins = bins, density = True)
N, bins, patches = ax_hist[1][1].hist(vf, bins = bins, density = True)

ax_hist[1][1].set_xticks(np.arange(bins_min, bins_max + bin_width, 20))
ax_hist[1][1].set_xticks(bins, minor = True)
ax_hist[1][1].set_xlim([-90, -60])

[ax.grid(False) for axs in ax_hist for ax in axs]








