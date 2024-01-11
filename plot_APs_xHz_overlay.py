# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 17:37:32 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import directories_win as directories
from PGFs import cc_APs_parameters, cc_th1Ap_parameters
import pandas as pd
from useful_functions import calc_time_series, butter_filter
import os
from cc_IF_functions import get_IF_data
from plotting_functions import get_colors, save_figures, get_figure_size, set_font_sizes, return_segments
import scipy as sc
import parameters

from matplotlib.collections import LineCollection

# %%

lookup_table = pd.read_excel(directories.table_dir + 'InVitro_Database.xlsx',
                             sheet_name="PGFs",
                             index_col='cell_ID')

# test cell E-092
cell_ID = 'E-092'

frequency_int = 50

frequency = str(frequency_int) + 'Hz'
    
# PGF to load
PGF = f'cc_APs_{frequency}'
PGF_parameters = cc_APs_parameters[frequency]

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

t_ms = calc_time_series(v_concat, SR)
t_s = calc_time_series(v_concat, SR, scale = 's')

# plt.plot(v_concat[0:15000])
# plt.show()


# filter voltage (to vf)
vf = butter_filter(v_concat, 
                   order = 3,
                   cutoff = 1e3,
                   sampling_rate = SR)


# %% find peaks

idx_peaks, dict_peak = sc.signal.find_peaks(vf, 
                                            prominence = parameters.min_peak_prominence, 
                                            distance = parameters.min_peak_distance * (SR_ms))



# %% create current array

i_hold = 0 # pA

i_stim = 100 # pA

i_start = cc_th1Ap_parameters['i_start']
i_delta = cc_th1Ap_parameters['i_delta']

i_steps= [i_stim] * n_steps

i = [None] * n_steps

for idx, i_stim in enumerate(i_steps):
    i_pre = np.full(int((SR_ms * PGF_parameters['t_pre'])), i_hold)
    i_stim = np.full(int((SR_ms * PGF_parameters['t_stim'])), i_stim)
    i_post = np.full(int((SR_ms * PGF_parameters['t_post'])-1), i_hold)
    
    i_step = np.concatenate((i_pre, i_stim, i_post))
    i[idx] = i_step



# %% split concatenate arrays back to steps wise 
# needs to occurs after filtering because of the filtering artifact

v = [None] * n_steps
t = [None] * n_steps
peaks = [None] * n_steps

step_dur = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']
step_points = step_dur * SR_ms

for idx in np.arange(0, n_steps, 1):
    start_idx = int(step_points * idx)
    stop_idx = int(start_idx + step_points - 1)

    v[idx] = vf[start_idx:stop_idx]
    t[idx] = t_ms[start_idx:stop_idx]
    peaks[idx] = [(idx_peak / SR_ms) - (step_dur * idx) for idx_peak in idx_peaks if idx_peak > start_idx and idx_peak < stop_idx]

# reduce dimensionality again to suit the eventplot function
peaks = [peak for peak in peaks]

# %%

darkmode_bool = True

color_dict = get_colors(darkmode_bool)


# Initialise figure
fig_ov, ax_ov = plt.subplots(1, 1, 
                             figsize = get_figure_size(width=78.418, height=52.833),
                             sharex = 'col')

# set font sizes
set_font_sizes()


# plot parameters
## initialise a line plot
segs_v = return_segments(t[0], v[:])
linesegments_v = LineCollection(segs_v, colors = 'grey', lw = 1.5)
linecollection = ax_ov.add_collection(linesegments_v)

line_v, = ax_ov.plot(t[0], v[-1], lw = 1.5, color = color_dict['color1'])


## initialise eventmarker
# events, = ax_ov.eventplot(peaks, color = 'grey', lineoffsets=40, linelengths=5, linewidth = 3)

## initialise text for threshold label
# text = ax_ov[0].text(x = 375,
#                       y = 100,
#                       s = '', 
#                       fontsize = 14,
#                       verticalalignment='center')

# axis settings
v_range = [-100, 50]

# ax_ov.set_xlabel('Time [ms]')
ax_ov.set_xlim([0, step_dur])
ax_ov.set_xticks(ticks = np.linspace(0, step_dur, 3))
ax_ov.set_xticks(ticks = np.linspace(0, step_dur, 5), minor = True)

# ax_ov.set_ylabel('Voltage [mV]')
ax_ov.set_ylim(v_range)

ax_ov.set_yticks(ticks = np.arange(v_range[0], v_range[1] + 1, 50),
                 labels = [None]*4)
ax_ov.set_yticks(np.arange(v_range[0], v_range[1] + 1, 25), minor = True)

ax_ov.grid(False)


fig_ov.tight_layout()


plt.show()


save_figures(fig_ov, f'{frequency}_APs_ov', directories.figure_dir, darkmode_bool)


print('Done!')











