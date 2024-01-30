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
from parameters.PGFs import cc_th1Ap_parameters, cc_APs_parameters
from parameters.parameters import min_peak_prominence, min_peak_distance, dvdt_threshold, AP_parameters
from parameters.directories_win import table_dir, vplot_dir, cell_descrip_dir, quant_data_dir

from functions.functions_ccIF import get_IF_data
from functions.functions_import import get_traceIndex_n_file
from functions.functions_useful import calc_time_series, butter_filter, calc_dvdt, calc_dvdt_padded
from functions.functions_spiketrains import get_AP_parameters
from functions.functions_plotting import get_colors, get_figure_size, save_figures


# %%

vplot_bool = False

# %% get cell IDs

table = pd.read_excel(table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')


# loop to create string to include all frequencies in query
query_str = ''

frequencies = list(cc_APs_parameters.keys())

for idx, frequency in enumerate(frequencies):
    PGF = 'cc_APs_' + frequency
    
    if idx > 0:
        query_str = query_str + ' and '
        
    query_str = query_str + f'{PGF}.notnull()'

query_str = query_str + ' and cc_th1AP.notnull()'    

# limit lookup table
lookup_table = table.query(query_str)

# cell IDs 
cell_IDs = list(lookup_table.index)


# %% initialize dataframes to populate

i_th_df = pd.DataFrame(columns = ['i_th_abs', 'i_th_rel'], index = cell_IDs)
firstAP_parameters_df = pd.DataFrame(columns = AP_parameters, index = cell_IDs)

# %%

PGF = 'cc_th1AP'

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
    
    
    # get i_hold
    lookup_table = pd.read_excel(table_dir + 'InVitro_Database.xlsx',
                                 sheet_name="V_or_I_hold",
                                 index_col='cell_ID')
    
    i_hold = lookup_table.at[cell_ID, PGF] # pA
        
    i_start = cc_th1Ap_parameters['i_start']
    i_delta = cc_th1Ap_parameters['i_delta']
    
    i_steps = np.arange(i_start, i_start + (i_delta * n_steps), i_delta)
    
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
    idx_peaks_s = [None] * n_steps
    
    step_dur = cc_th1Ap_parameters['t_pre'] + cc_th1Ap_parameters['t_stim'] + cc_th1Ap_parameters['t_post']
    step_points = step_dur * SR_ms
    pre_points = int(cc_th1Ap_parameters['t_pre'] * SR_ms)
    poststim_points = int((cc_th1Ap_parameters['t_stim'] + 10) * SR_ms)
    
    for idx in np.arange(0, n_steps, 1):
        start_idx = int(step_points * idx) + pre_points
        stop_idx = int(start_idx + poststim_points)
        
        v[idx] = vf[start_idx:stop_idx]
        t[idx] = t_ms[start_idx:stop_idx]
        idx_peaks_s[idx] = [int(idx_peak - (step_points * idx)) for idx_peak in idx_peaks if idx_peak > start_idx and idx_peak < stop_idx]
        peaks[idx] = [(idx_peak / SR_ms) - (step_dur * idx) for idx_peak in idx_peaks if idx_peak > start_idx and idx_peak < stop_idx]
    
    
    # %% where is the threshold
    
    idx_th = next(p for p, idx_peak in enumerate(peaks) if len(idx_peak) > 0)
    
    
    
    i_th_abs = i_steps[idx_th]
    
    i_th_df.at[cell_ID, 'i_th_abs'] = i_th_abs
    i_th_df.at[cell_ID, 'i_th_rel'] = i_th_abs - i_hold
    
    
    
    
    # %% build dataframe with v steps to threshold
    
    idx_start = int(cc_th1Ap_parameters['t_pre']) * int(SR_ms)
    idx_stop = int(cc_th1Ap_parameters['t_pre'] + cc_th1Ap_parameters['t_stim'] + 10.) * int(SR_ms)
    
    v_th = v[idx_th]
    i_th = i[idx_th][idx_start:idx_stop]
    t = calc_time_series(data = v_th, sampling_rate=SR)
    
    dvdt_th = calc_dvdt_padded(v[idx_th], t)
    
    
    
    v_df = pd.DataFrame(index = t)
    dvdt_df = pd.DataFrame(index = t)
    
    
    
    for idx in np.arange(idx_th):
        # limit v to the stim and 10 ms post stim
        v_df[idx] = v[idx]
        
        dvdt_cur = calc_dvdt_padded(v[idx], t)
     
        dvdt_df[idx] = dvdt_cur
     
     
    
    # %% get first AP parameters
    
    th_spike_idx = [int(idx - (cc_th1Ap_parameters['t_pre'] * SR_ms)) for idx in idx_peaks_s[idx_th]]
    
    AP_params = get_AP_parameters(v_th, th_spike_idx, 
                                  SR = SR,
                                  dvdt_threshold = dvdt_threshold,
                                  t_pre = 2,
                                  t_post = 10)
    
    firstAP_parameters_df.loc[cell_ID] = AP_params.iloc[0]
    
    firstAP_parameters_df.at[cell_ID, 'SR_ms'] = SR_ms
    
    
    # %% build first AP dataframe
    
    firstAP_df = pd.DataFrame({'i' : i_th,
                               'v' : v_th,
                               'dvdt' : dvdt_th},
                              index = t)
    
        
    
    # %% test plots
    
    if vplot_bool:
    
        darkmode_bool = True
        
        colors_dict = get_colors(darkmode_bool)
        
        fig_1stAP, ax_1stAP = plt.subplots(1, 2, 
                                           layout = 'constrained',
                                           figsize = get_figure_size())
        fig_1stAP.suptitle(cell_ID)
        
        for idx in np.arange(0, idx_th):
            ax_1stAP[0].plot(t, v_df[idx],
                             c = 'grey')
        
        ax_1stAP[0].plot(t, v_th,
                         c = colors_dict['color1'])
            
        ax_1stAP[0].hlines(y = AP_params['v_HM'][0],
                            xmin = AP_params['t1_HM'][0],
                            xmax = AP_params['t2_HM'][0], 
                            colors = 'r')
        
        ax_1stAP[0].vlines(x = AP_params['t_peaks'][0],
                            ymin = AP_params['v_threshold'][0],
                            ymax = AP_params['v_threshold'][0] + AP_params['v_amplitude'][0], 
                            colors = 'r')
        
        ax_1stAP[0].set_xlim([0, 20])
        ax_1stAP[0].set_xticks(np.linspace(0, 20, 3))
            
        ax_1stAP[0].set_ylim([-100, 60])
        ax_1stAP[0].set_yticks(np.arange(-100, 60 + 1, 20))
        ax_1stAP[0].set_yticks(np.arange(-100, 60 + 1, 5), minor = True)
        
        ax_1stAP[0].set_ylabel('Membrane potential [mV]')
        
        ax_1stAP[0].set_xlabel('Time [ms]')
        
        
        ax_1stAP[1].plot(v_df, dvdt_df,
                         c = 'grey')
        
        ax_1stAP[1].plot(v_th, dvdt_th,
                         c = colors_dict['color1'])
        
        ax_1stAP[1].set_ylabel('Rate of membrane potential change\n[mV/ms]')
        ax_1stAP[1].set_ylim([-150, 250]) 
        
        ax_1stAP[1].set_xlabel('Membrane potential [mV]')
        ax_1stAP[1].set_xlim([-100, 60])   
        
        [ax.grid(False) for ax in ax_1stAP]    
           
        # save figure as verification plot
        fig_path = os.path.join(vplot_dir, 'th1AP')
        if not os.path.exists(fig_path):
            os.mkdir(fig_path)
        
        save_figures(fig_1stAP, f'cc_th1AP-{cell_ID}-1st_AP-phaseplane_plot', fig_path, darkmode_bool)
  
        plt.show()
        

# %% save AP parameters


    fst_AP_dir = os.path.join(quant_data_dir, '1stAP')
    if not os.path.exists(fst_AP_dir):
        os.mkdir(fst_AP_dir)
    
    firstAP_df.to_excel(os.path.join(fst_AP_dir, f'ccth1AP-1stAP-{cell_ID}.xlsx'), index_label = 't')


firstAP_parameters_df.to_excel(os.path.join(cell_descrip_dir, 'ccth1AP-fst_AP_parameters.xlsx'), index_label = 'cell_ID')

i_th_df.to_excel(os.path.join(cell_descrip_dir, 'ccth1AP-fst_AP_i.xlsx'), index_label = 'cell_ID')




        