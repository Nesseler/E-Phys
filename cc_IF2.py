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
import statistics as stats


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

figure_dir = '/Users/moritznesseler/local E-Phys/figures'

darkmode_bool = True

# %%


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

i_max = 350 # pA
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
    
    # convert i_inputs to float to pad with NaNs
    i_inputs = [float(i) for i in i_inputs]
    i_input_max = i_inputs[-1]


    ### SINGLE VALUES PER IF ###
    # R_input, tau_mem, rheobase, f_max
    IF_values.insert(cell_idx, cell_ID, pd.Series([R_input, Tau_mem, rheobase, f_max, i_input_max, i_hold_r], 
                                                  index=['R_input', 'Tau_mem', 'rheobase', 'f_max', 'i_input_max', 'i_hold']))
    
    ### ARRAY LIKE PER IF ###
        # i_inputs, n_spikes, f_spikes
    
    # f_spikes_df = pd.DataFrame({'i_inputs' : possible_i_inputs},
    #                            index = possible_i_inputs)
    n_pad_before = abs(int(i_hold_r / i_stepdelta))
    n_pad_after = abs(int(len(possible_i_inputs) - (n_steps + n_pad_before)))

    
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

sbn.swarmplot(x=[1.]*len(IF_values), y = IF_values['f_max'])

plt.show()


# %%

column_names = list(IF_values.columns)

f_cutoff = 50

fast_spiking = IF_values.query(f'f_max > {f_cutoff}')
slow_spiking = IF_values.query(f'f_max <= {f_cutoff}')


# f_spikes_df[fast_spiking.index].plot()
# plt.xlim([i_start, i_max])
# plt.show()

# f_spikes_df[slow_spiking.index].plot()
# plt.xlim([i_start, i_max])
# plt.show()

test = n_spikes_df.last_valid_index

n_steps_pre_max = 5

indices_depol_block = []
indices_no_depol_block = []


for col_head in n_spikes_df[slow_spiking.index]:
    # print(n_spikes_df[col_head])
    
    # print(np.max(n_spikes_df[col_head]))
    
    print(col_head)
    
    i_max_index = IF_values.at[col_head ,'i_input_max']
    
    i_pre_index = i_max_index - (n_steps_pre_max * i_stepdelta)
    
    n_spikes_mean_i_max = round(np.mean(n_spikes_df[col_head].loc[i_pre_index:i_max_index]))
    
    n_max = np.max(n_spikes_df[col_head])
    
    print(n_max, n_spikes_mean_i_max, n_spikes_mean_i_max - n_max)
    
   
    if (n_max - n_spikes_mean_i_max) > 10:
        indices_depol_block.append(col_head)
    else:
        indices_no_depol_block.append(col_head)



# f_spikes_df[indices_depol_block].plot()
# plt.xlim([i_start, i_max])
# plt.show()



# %%

darkmode_bool = True

colors_dict = get_colors(darkmode_bool)

add_parameter_ls = ['f_max', 'R_input', 'Tau_mem', 'rheobase']

for add_parameter in add_parameter_ls:
    
    fig_IFs, axs_IFs = plt.subplots(2, 3,
                                    figsize = get_figure_size(),
                                    layout = 'tight')
                                    # sharey='row')
    
    small_font_size = 14
    large_font_size = 16
    
    plt.rc('font', size = small_font_size)
    plt.rc('axes', titlesize = large_font_size, 
                   labelsize = large_font_size,
                   linewidth = 0.5)
    plt.rc('xtick', labelsize = large_font_size)
    plt.rc('ytick', labelsize = large_font_size)
    plt.rc('lines', linewidth = 2)
    
    
    
    gs = axs_IFs[0, 1].get_gridspec()
    # remove the underlying axes
    for ax in axs_IFs[0][:2]:
        ax.remove()
        
    ax_all_IFs = fig_IFs.add_subplot(gs[0,:2])
    
    
    
    
    ax_all_IFs.plot(f_spikes_df[fast_spiking.index], color = colors_dict['color1'])
    ax_all_IFs.plot(f_spikes_df[indices_depol_block], color = colors_dict['color2'])
    ax_all_IFs.plot(f_spikes_df[indices_no_depol_block], color = colors_dict['color3'])
    
    
    
    # axs_IFs[0][2].scatter(x = [1] * len(IF_values), 
    #                       y = IF_values[add_parameter])
    
    colors_fs = {i : colors_dict['color1'] for i in list(fast_spiking.index)}
    colors_ssdb = {i : colors_dict['color2'] for i in list(indices_depol_block)}
    colors_ssndb = {i : colors_dict['color3'] for i in list(indices_no_depol_block)}
    
    colors_all = {}
    
    for i, cell in enumerate(all_cell_IDs):
        if cell in fast_spiking.index:
            colors_all[cell] = colors_dict['color1']
        elif cell in indices_depol_block:
            colors_all[cell] = colors_dict['color2']
        elif cell in indices_no_depol_block:
            colors_all[cell] = colors_dict['color3']
    
    
    x_values = [1 + stats.random.random() for i in all_cell_IDs]
    
    # axs_IFs[0][2].scatter(x = x_values, 
    #                       y = IF_values[add_parameter],
    #                       c = colors_all.values())
    
    
    # sbn.swarmplot(x = [1.]*len(IF_values), 
    #               y = IF_values[add_parameter],
    #               ax = axs_IFs[0][2], 
    #               palette = colors_all.values())
    
    
    sbn.swarmplot(x=[0.]*len(IF_values), 
                  y = IF_values[add_parameter], 
                  ax = axs_IFs[0][2], 
                  color = colors_dict['primecolor'])
    
    sbn.swarmplot(x=[1.]*len(fast_spiking.index), 
                  y = IF_values[add_parameter].loc[fast_spiking.index], 
                  ax = axs_IFs[0][2], 
                  color = colors_dict['color1'])
    
    sbn.swarmplot(x=[2.]*len(indices_depol_block), 
                  y = IF_values[add_parameter].loc[indices_depol_block], 
                  ax = axs_IFs[0][2], 
                  color = colors_dict['color2'])
    
    sbn.swarmplot(x=[3.]*len(indices_no_depol_block), 
                  y = IF_values[add_parameter].loc[indices_no_depol_block], 
                  ax = axs_IFs[0][2], 
                  color = colors_dict['color3'])
    
    
    axs_IFs[1][0].plot(f_spikes_df[fast_spiking.index], color = colors_dict['color1'])
    axs_IFs[1][1].plot(f_spikes_df[indices_depol_block], color = colors_dict['color2'])
    axs_IFs[1][2].plot(f_spikes_df[indices_no_depol_block], color = colors_dict['color3'])
    
    # axs_IFs[0][2].plot(f_spikes_df[indices_no_depol_block])
    
    
    # axs_IFs[0][0].get_shared_x_axes().join(axs_IFs[1][0],
    #                                        axs_IFs[1][0],
    #                                        axs_IFs[1][1],
    #                                        axs_IFs[1][2])
    
    
    IF_axes = [ax_all_IFs, axs_IFs[1][0], axs_IFs[1][1], axs_IFs[1][2]]
    IF_titles = ['Fast spiking', 'Slow spiking\nDepolarisation block', 'Slow spiking\n No depolarisation block']
    
    for i, axis in enumerate(IF_axes):
        axis.set_ylim([0, 80])
        axis.set_yticks(np.arange(0, 80+1, 20))
        axis.set_yticks(np.arange(0, 80+1, 5), minor = True)
        axis.set_ylabel('AP freq. [Hz]')
        
        axis.set_xlim([-50, i_max])
        axis.set_xticks(np.arange(0, i_max, 100))
        axis.set_xticks(np.arange(0, i_max+1, 25), minor = True)
        axis.set_xlabel('Injected current [pA]')
        
        if i > 0:
            axis.set_title(IF_titles[i-1])
    
    # fig_IFs.supxlabel('Injected current [pA]')
        
    
    ax_all_IFs.set_xticks(np.arange(-50, i_max + 1, 50), minor = True)
    
    fig_IFs.align_labels(axs_IFs[1])
    
    # parameter plot
    
    if add_parameter == 'f_max':
        # f_max
        axs_IFs[0][2].set_ylim([0, 80])
        axs_IFs[0][2].set_ylabel('Maximal AP freq. [Hz]')
        
    elif add_parameter == 'R_input':
        axs_IFs[0][2].set_ylim([0, 1300])
        axs_IFs[0][2].set_ylabel('Input resistance [MΩ]')
    
    elif add_parameter == 'Tau_mem':
        axs_IFs[0][2].set_ylim([0, 50])
        axs_IFs[0][2].set_ylabel('Membrane time\nconstant [ms]')
    
    elif add_parameter == 'rheobase':
        axs_IFs[0][2].set_ylim([i_start, 150])
        axs_IFs[0][2].set_yticks(np.arange(i_start, 150+1, 50))
        axs_IFs[0][2].set_ylabel('Rheobase [pA]')
    
    
    axs_IFs[0][2].set_xticks([0, 1, 2, 3],['all', 'fs', 'ss\ndb', 'ss\nndb'])
     
    [ax.grid(False) for rows in axs_IFs for ax in rows]
    ax_all_IFs.grid(False)
    
    save_figures(fig_IFs, f'IF_categories_{add_parameter}', figure_dir, darkmode_bool)








