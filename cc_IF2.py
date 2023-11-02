#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:28:20 2023

@author: moritznesseler
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
from useful_functions import calc_time_series, calc_dvdt, butter_filter, save_figures, get_sampling_rate, get_data
import os.path
import scipy as sc
import seaborn as sbn
import warnings 


from plotting_functions import get_figure_size, get_colors

from cc_IF_functions import get_IF_data


# %%

# TODO
    # function for input current
        # get holding current
        # form relative input current array starting at -50 pA
    
    # R_input & tau_mem
        # function for exponential fit
        # option: figure output of every fit
    
    # tau_mem
    
    # rheobase
    
    # IF curve
    
    
    
    
lookup_table = pd.read_csv('/Users/moritznesseler/local E-Phys/cc_IF.csv', 
                           delimiter=';',
                           index_col='cell_ID')


data_folder = '/Users/moritznesseler/local E-Phys'

figure_dir = '/Users/moritznesseler/local E-Phys'

darkmode_bool = True


### parameters of step_IF

pre_post_dur = 0.250 #s
pulse_dur = 1 #s
i_stepdelta = 5 #pA
i_start = -50 #pA


### SPIKES

# set parameters to find peaks
min_peak_prominence = 20 #(mV)
min_peak_distance = 1 #ms


### INPUT RESISTANCE AND MEMBRANE TIME CONSTANT ###
# time window for exponential fit: 200 ms pre & post stim begin
time_interval = 0.150 #s

# define function of exponential fit
def exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c

# set a guess for exponential fit
popt_guess = [10, 0.005, -85]

# convert indices of dataframe to list to loop through
all_cell_IDs = lookup_table.index.to_list()


### LOOP THROUGH ALL CELLS ###

### SINGLE VALUES PER IF ###
# R_input, tau_mem, rheobase, f_max
IF_values = pd.DataFrame()

### ARRAY LIKE PER IF ###
    # i_inputs, n_spikes, f_spikes
n_spikes_df = pd.DataFrame()
f_spikes_df = pd.DataFrame()
t_spikes_df = pd.DataFrame()

i_max = 400 # pA
possible_i_inputs = np.arange(start = i_start,
                              stop = i_max + 1,
                              step = i_stepdelta)


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
    
    # get IF data form file
    i, v, t, SR, n_steps = get_IF_data(data_file_path_str, traceIndex, 'ms')
    
    # calculate pre, pulse, and post indices
    pre_idx = np.arange(0, int(pre_post_dur * SR))
    pulse_idx = np.arange(pre_idx[-1]+1, int((pre_post_dur + pulse_dur) * SR))
    post_idx = np.arange(pulse_idx[-1]+1, int((pre_post_dur + pulse_dur + pre_post_dur) * SR))
    
    # calc current steps from holding potential
    i_hold = np.mean(i[:,:int(0.25*SR)])
    i_hold_r = 5 * round(i_hold/5)
    
    #subtract holding potential from start current
    i_start_wo_hold = i_start - i_hold_r
    i_stop = i_start + (i_stepdelta * n_steps)
    
    i_inputs = np.arange(start = i_start_wo_hold, 
                         stop = i_stop + 1,
                         step = i_stepdelta)
    
    ### SPIKES ###
    
    # initialize lists
    r_inputs = []
    tau_mems = []
    t_spikes = []
    n_spikes = []
    f_spikes = []
    
    # indices for limiting data for exponential fit
    start_idx = int((pre_post_dur - time_interval) * SR)
    step_start_idx = int(pre_post_dur * SR)
    stop_idx = int((pre_post_dur + time_interval) * SR)
    
    # stim indices 
    idx_pulse_start = int(pre_post_dur * SR)
    idx_pulse_stop = int(((pre_post_dur + pulse_dur) * SR ) + 1)
                        

    # number of first steps to use for calculation of R_input and tau_mem
    n_steps_for_R = 3
    
    # loop through all steps
        # get spike times
        # get number of spikes
        # get spike frequency
        # get R_input
        # get tau_mem
    
    for step_idx in np.arange(n_steps):
        
        # filter voltage (to vf)
        vf = butter_filter(v[step_idx], 
                          order = 3,
                          cutoff = 1e3,
                          sampling_rate = SR)
        
        # limit data array to just the stimulus
        vl = vf[idx_pulse_start:idx_pulse_stop]
        
        # find peaks as spikes
        idx_peaks, dict_peak = sc.signal.find_peaks(vl, 
                                                    prominence = min_peak_prominence, 
                                                    distance = min_peak_distance * (SR/1e3))
        
        # get times of spikes
        t_peaks = np.divide(idx_peaks, (SR/1e3))
        t_spikes.append(t_peaks)
        
        # get number of spikes
        n_peaks = len(idx_peaks)
        n_spikes.append(n_peaks)
        
        # get frequency over the 
        if n_peaks > 0:
            f_peaks = n_peaks / pulse_dur
        else:
            f_peaks = 0   
        f_spikes.append(f_peaks)
            
        ### opt: PLOT ###
        # plt.plot(vf)
        # plt.eventplot(idx_peaks, color = 'r')
        # plt.pause(0.4)
        # plt.show()
        
        
        # loop through steps
        if step_idx < n_steps_for_R:
        
            # define data pre and post stim begin
            v_pre = vf[start_idx:step_start_idx]
            v_post = vf[step_start_idx:stop_idx]
        
            # calc x for v_post
            x_post = np.arange(len(v_post))
        
            # exponential fit
            try:
                popt, pcov = sc.optimize.curve_fit(exp_func, x_post, v_post, p0=popt_guess, maxfev = 1000)
           
                ### R_INPUT ###
                #R_input = ∆U/∆I
                #∆U = a = popt[0], ∆I = 20 (for first step)
                v_pre_mean = np.mean(v_pre)
                delta_v = (popt[2] - np.mean(v_pre))
            
                #delta I
                delta_i = i_inputs[step_idx]
            
                # calc r_input for current step
                r_input = ( delta_v / delta_i ) * 1e3 #in MOhm
            
                # append to list
                r_inputs.append(r_input)
            
                ### MEMBRANE TIME CONSTANT
                # tau_mem
                # time it takes the potential to reach 1 - (1/e) (~63%) of the max voltage
            
                # calc 1 - 1/e
                tau_perc_value = 1-(1/np.exp(1))
                
                # calc max voltage delta
                delta_v_63 = delta_v * tau_perc_value
                
                # calc 63 % of max voltage
                v_tau = v_pre_mean+delta_v_63
                
                # calc time (indices first) it takes to reach v_tau
                idx_63 = - (np.log((v_tau - popt[2]) / (popt[0]))) / (popt[1])
                tau_mem = idx_63 / (SR / 1e3) 
                
                # append to list
                tau_mems.append(tau_mem)
            
            except RuntimeError:
                print(f'{cell_ID} step number {step_idx+1} has been omitted')
        
            ### PLOT ###
            # fig1_R_input, ax1_R_input = plt.subplots(1,1)
            # ax1_R_input.plot(v_post)
            # plt.plot(x_post, exp_func(x_post, *popt), 'r--')
            # ax1_R_input.hlines(v_tau, 0, idx_63, colors='k', linestyle='--')
            # ax1_R_input.vlines(idx_63, v_pre_mean, v_tau, colors='k', linestyle='--')
            # ax1_R_input.set_ylim([-110, -80])
            # plt.pause(0.5)
            # plt.show()
     
    
    # get maximal spiking frequency
    f_max = np.max(f_spikes)
    
    # calc r_input and tau_mem as mean of first 3 steps
    R_input = np.mean(r_inputs)
    Tau_mem = np.mean(tau_mems)
    
    # get first index in number of spikes where there more than 0 spikes
    rheobase_idx = next(idx for idx, n_spike in enumerate(n_spikes) if n_spike > 0)
    
    # get rheobase relative to holding
    rheobase_rel = i_inputs[rheobase_idx]
    
    # add holding current
    rheobase = rheobase_rel + i_hold_r
    

    ### SINGLE VALUES PER IF ###
    # R_input, tau_mem, rheobase, f_max
    IF_values.insert(cell_idx, cell_ID, pd.Series([R_input, Tau_mem, rheobase, f_max], 
                                                  index=['R_input', 'Tau_mem', 'rheobase', 'f_max']))
    
    ### ARRAY LIKE PER IF ###
        # i_inputs, n_spikes, f_spikes
    
    # f_spikes_df = pd.DataFrame({'i_inputs' : possible_i_inputs},
    #                            index = possible_i_inputs)
    
    n_pad_before = abs(int(i_hold_r / i_stepdelta))
    n_pad_after = abs(int(len(possible_i_inputs) - (n_steps + n_pad_before)))
    
    # convert i_inputs to float to pad with NaNs
    i_inputs = [float(i) for i in i_inputs]
    
    # numbers of spikes
    n_spikes = [float(n) for n in n_spikes]
    n_spikes = np.pad(n_spikes, 
                      pad_width = (n_pad_before, n_pad_after), 
                      mode = 'constant', 
                      constant_values = (np.nan,))
    
    n_spikes_df.insert(cell_idx, cell_ID, pd.Series(n_spikes, index = possible_i_inputs))
    
    # spiking frequencies
    f_spikes = np.pad(f_spikes, 
                      pad_width = (n_pad_before, n_pad_after), 
                      mode = 'constant', 
                      constant_values = (np.nan,))
    
    f_spikes_df.insert(cell_idx, cell_ID, pd.Series(f_spikes, index = possible_i_inputs))


    print(f'{cell_idx+1} of {len(all_cell_IDs)}')


IF_values = IF_values.transpose()


# %%

f_spikes_df.plot()
plt.xlim([i_start, i_max])
plt.show()

# %% 

sbn.swarmplot(x=[1.]*len(IF_values), y = IF_values['rheobase'])

plt.show()
