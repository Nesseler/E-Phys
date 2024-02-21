# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:16:29 2024

@author: nesseler
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

# custom directories & parameters
from parameters.directories_win import table_file
from parameters.parameters import min_peak_prominence, min_peak_distance, dvdt_threshold, AP_parameters
from parameters.PGFs import cc_IF_parameters
from getter.get_cell_IDs import get_cell_IDs_one_protocol, get_cell_IDs_all_ccAPfreqs

from functions.functions_constructors import construct_current_array
from functions.functions_ccIF import get_IF_data
from functions.functions_import import get_traceIndex_n_file
from functions.functions_useful import calc_time_series, butter_filter, calc_dvdt, calc_dvdt_padded
from functions.functions_plotting import get_colors, get_figure_size, save_figures, plot_t_vs_v
from functions.functions_extractspike import get_AP_parameters


# %% settings

vplot_bool = True

darkmode_bool = False
colors_dict, _ = get_colors(darkmode_bool)


# %% load data

# protocol 
PGF = 'cc_IF'

# get cell IDs
# cell_IDs = get_cell_IDs_one_protocol(PGF)
cell_IDs = get_cell_IDs_all_ccAPfreqs()

# get hold current as table
I_hold_table = pd.read_excel(table_file, sheet_name="V_or_I_hold", index_col='cell_ID').loc[cell_IDs, :]


# test cekk 
cell_IDs = ['E-109']

for cell_ID in cell_IDs:
    
    print(f'Started: {cell_ID}')
    
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)
    
    # get IF data form file
    i, v, t, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')
    
    # sampling rate in ms
    SR_ms = SR / 1e3
    
    # concatenate individual steps
    n_points = int(np.shape(i)[0] * np.shape(i)[1])
    
    v_concat = v.flatten() 
    
    t_ms = calc_time_series(v_concat, SR)
    t_s = calc_time_series(v_concat, SR, scale = 's')
    
    # filter voltage (to vf)
    vf = butter_filter(v_concat, 
                       order = 3,
                       cutoff = 1e3,
                       sampling_rate = SR)

    ### construct voltage dataframe
    i_hold = I_hold_table.at[cell_ID, PGF]
    i = construct_current_array(i_hold = i_hold,
                                n_steps = n_steps,
                                parameters_dict = cc_IF_parameters,
                                SR_ms = SR_ms)
    
    # split concatenate arrays back to steps wise 
    # needs to occurs after filtering because of the filtering artifact
    
    v = [None] * n_steps
    t = [None] * n_steps
    peaks = [None] * n_steps
    idx_peaks_s = [None] * n_steps
    
    step_dur = cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim'] + cc_IF_parameters['t_post']
    step_points = step_dur * SR_ms
    
    pre_points = int(cc_IF_parameters['t_pre'] * SR_ms)
    pre_n_stim_points = int((cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim']) * SR_ms)
    
    
    # loop through steps to limit voltage trace
    for idx in np.arange(0, n_steps, 1):
        
        # calc start and stop indices for step
        start_idx = int(step_points * idx)
        stop_idx = int(start_idx + step_points)
        
        # set voltage trace of step
        v_step = vf[start_idx:stop_idx]
        
        # assign voltage and time traces to arrays
        v[idx] = v_step
        t[idx] = t_ms[start_idx:stop_idx]
    
        # find peaks
        idc_peaks, dict_peak = sc.signal.find_peaks(v_step, 
                                                    prominence = min_peak_prominence, 
                                                    distance = min_peak_distance * (SR_ms))

            
        # limit spike indices to stimulation time period
        idc_peaks = [idx_peak for idx_peak in idc_peaks if idx_peak > pre_points and idx_peak <= pre_n_stim_points]
        
        # verification plot for detection of spikes in each step
        if vplot_bool:
            plt.plot(v_step, linewidth = 1, c = colors_dict['primecolor'])
            plt.eventplot(idc_peaks, lineoffsets=60, colors = 'r', linelengths=5)
            plt.ylim([-100, 75])
            plt.grid(False)
            plt.show()        
    

    
















