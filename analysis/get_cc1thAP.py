# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 18:41:12 2024

@author: nesseler
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

# custom directories & parameters
from parameters.PGFs import cc_th1Ap_parameters
from parameters.parameters import min_peak_prominence, min_peak_distance

from functions.functions_ccIF import get_IF_data
from functions.functions_import import get_traceIndex_n_file
from functions.functions_useful import calc_time_series, butter_filter

# %%

# PGF to load
PGF = 'cc_th1AP'
cell_ID = 'E-069'

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


# %% find peaks

idx_peaks, dict_peak = sc.signal.find_peaks(vf, 
                                            prominence = min_peak_prominence, 
                                            distance = min_peak_distance * (SR_ms))



# %% create current array

i_hold = 0 # pA

i_start = cc_th1Ap_parameters['i_start']
i_delta = cc_th1Ap_parameters['i_delta']

i_steps= np.arange(i_start, i_start + (i_delta * n_steps), i_delta)

i = [None] * n_steps

for idx, i_stim in enumerate(i_steps):
    i_pre = np.full(int((SR_ms * cc_th1Ap_parameters['t_pre'])), i_hold)
    i_stim = np.full(int((SR_ms * cc_th1Ap_parameters['t_stim'])), i_stim)
    i_post = np.full(int((SR_ms * cc_th1Ap_parameters['t_post'])-1), i_hold)
    
    i_step = np.concatenate((i_pre, i_stim, i_post))
    i[idx] = i_step



# %% split concatenate arrays back to steps wise 
# needs to occurs after filtering because of the filtering artifact

v = [None] * n_steps
t = [None] * n_steps
peaks = [None] * n_steps

step_dur = cc_th1Ap_parameters['t_pre'] + cc_th1Ap_parameters['t_stim'] + cc_th1Ap_parameters['t_post']
step_points = step_dur * SR_ms

for idx in np.arange(0, n_steps, 1):
    start_idx = int(step_points * idx)
    stop_idx = int(start_idx + step_points - 1)
    
    v[idx] = vf[start_idx:stop_idx]
    t[idx] = t_ms[start_idx:stop_idx]
    peaks[idx] = [(idx_peak / SR_ms) - (step_dur * idx) for idx_peak in idx_peaks if idx_peak > start_idx and idx_peak < stop_idx]





plt.plot(v[7])
    
    
    
    
    
    
    
    