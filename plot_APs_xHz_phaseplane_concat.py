# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 19:01:11 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import pandas as pd
import os
import scipy as sc

# custom directories & parameters
import parameters.directories_win as directories
from parameters.parameters import min_peak_prominence_ccIF, min_peak_distance_ccIF
from parameters.PGFs import cc_APs_parameters, cc_th1Ap_parameters

from functions.functions_ccIF import get_IF_data
from functions.functions_useful import calc_time_series, butter_filter, calc_dvdt
from functions.functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes
from functions.functions_spiketrains import get_colorcode

# %%

lookup_table = pd.read_excel(directories.table_file,
                             sheet_name="PGFs",
                             index_col='cell_ID')

# test cell E-092
cell_ID = 'E-092'

frequency_int = 75

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
                                            prominence = min_peak_prominence_ccIF, 
                                            distance = min_peak_distance_ccIF * (SR_ms))



# %% create current array

i_hold = 0 # pA

i_stim = 100 # pA

i_start = cc_th1Ap_parameters['i_start']
i_delta = cc_th1Ap_parameters['i_delta']

i_steps= [i_stim] * n_steps

i = []

for idx, i_stim in enumerate(i_steps):
    i_pre = np.full(int((SR_ms * PGF_parameters['t_pre'])), i_hold)
    i_stim = np.full(int((SR_ms * PGF_parameters['t_stim'])), i_stim)
    i_post = np.full(int((SR_ms * PGF_parameters['t_post'])), i_hold)
    
    i_step = np.concatenate((i_pre, i_stim, i_post))
    i = np.concatenate((i, i_step))


#dvdt
dvdt = calc_dvdt(vf, t_ms)

# %% split concatenate arrays back to steps wise 
# needs to occurs after filtering because of the filtering artifact

# v = [None] * n_steps
# dvdt_s = [None] * n_steps
# t = [None] * n_steps
# peaks = [None] * n_steps

step_dur = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']
total_dur = step_dur * n_steps
# step_points = step_dur * SR_ms

# for idx in np.arange(0, n_steps, 1):
#     start_idx = int(step_points * idx)
#     stop_idx = int(start_idx + step_points - 1)

#     v[idx] = vf[start_idx:stop_idx]
#     t[idx] = t_ms[start_idx:stop_idx]
#     dvdt_s[idx] = dvdt[start_idx:stop_idx]
#     peaks[idx] = [(idx_peak / SR_ms) - (step_dur * idx) for idx_peak in idx_peaks if idx_peak > start_idx and idx_peak < stop_idx]

# # reduce dimensionality again to suit the eventplot function
# peaks = [peak for peak in peaks]



# %% limit dataset to stim



# stim_start_idx = int((PGF_parameters['t_pre'] - 10) * SR_ms)
# stim_stop_idx = int((PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + 10) * SR_ms) + 1

# if stim_start_idx < 0:
#     stim_start_idx = 0

# v = v[0][stim_start_idx:stim_stop_idx]
# i = i[0][stim_start_idx:stim_stop_idx]
# t = t[0][stim_start_idx:stim_stop_idx]
# dvdt = dvdt_s[0][stim_start_idx:stim_stop_idx]
data_fc = t_ms[200:]
v = vf[200:]
t = t_ms[200:]
i = i[200:]
dvdt = dvdt[200:]



# %%

darkmode_bool = True

color_dict, _ = get_colors(darkmode_bool)

cmap = 'viridis_r'


v_range = [-100, 50]

cc_phaseplane, cc_ppaxs = plt.subplots(2,2, 
                                       layout = 'constrained',
                                       figsize = get_figure_size(),
                                       gridspec_kw={'height_ratios': [1,4]},
                                       sharex='col')

cc_phaseplane.delaxes(cc_ppaxs[0][1])

set_font_sizes()



# function to colorcode time of timeseries









# norm = mtl.colors.Normalize(0, data_fc.max())
norm = mtl.colors.Normalize(0, total_dur)

# i vs t
cc_ppaxs[0][0].plot(t, i, color_dict['color2'])

cc_ppaxs[0][0].set_ylabel('Current [pA]')
cc_ppaxs[0][0].set_ylim([-50, 200])
cc_ppaxs[0][0].set_yticks(ticks = np.arange(-50, 200+1, 50), 
                     labels = [None, 0, None, 100, None, 200])
cc_ppaxs[0][0].set_yticks(np.arange(-50, 200+1, 50), minor = True)

# v vs t
# cc_ppaxs[1][0].plot(t, v, color_dict['color1'])

lc = get_colorcode(t, v, data_fc, norm = norm, cmap=cmap)
line = cc_ppaxs[1][0].add_collection(lc)

cc_ppaxs[1][0].set_xlim(0, t.max()+1)
cc_ppaxs[1][0].set_xlabel('Time [ms]')
cc_ppaxs[1][0].set_ylim(v_range)
cc_ppaxs[1][0].set_ylabel('Membrane potential [mV]')


# phase plane
# cc_ppaxs[1][1].plot(v[1:], dvdt, color_dict['color1'])

lc = get_colorcode(v[1:], dvdt, data_fc, norm = norm, cmap=cmap)
line = cc_ppaxs[1][1].add_collection(lc)

cc_phaseplane.colorbar(line, ax=cc_ppaxs[1][0], orientation = 'horizontal',
                       ticks = np.linspace(0, total_dur, 5))

cc_ppaxs[1][1].set_ylim([-100, 250])
cc_ppaxs[1][1].set_xlim(v_range)
cc_ppaxs[1][1].set_ylabel('Rate of membrane potential change\n[mV/ms]')
cc_ppaxs[1][1].set_xlabel('Membrane potential [mV]')

[ax.grid(False) for row_axs in cc_ppaxs for ax in row_axs]

plt.show()

save_figures(cc_phaseplane, f'{frequency}_APs_phaseplane_concat', directories.figure_dir, darkmode_bool)






















