# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 18:01:22 2024

@author: nesseler
"""

# E-092: use for testing

import directories_win as directories
import pandas as pd
import os
from cc_IF_functions import get_IF_data
import matplotlib.pyplot as plt
import numpy as np
from useful_functions import calc_time_series, butter_filter, calc_dvdt
import parameters
from PGFs import cc_APs_parameters
import scipy as sc
from plotting_functions import get_colors, save_figures, get_figure_size
from spiketrains_functions import plot_vt_n_dvdtv_colorcoded, get_colorcode, plot_voltage_v_time, phase_plane_plot, get_AP_parameters

# %%

table = pd.read_excel(directories.table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')


frequencies = list(cc_APs_parameters.keys())

# loop to create string to include all frequencies in query
query_str = ''

for idx, frequency in enumerate(frequencies):
    PGF = 'cc_APs_' + frequency
    
    if idx > 0:
        query_str = query_str + ' and '
        
    query_str = query_str + f'{PGF}.notnull()'
    

# limit lookup table
lookup_table = table.query(query_str)

# %% 

# test cell E-092
cell_ID = 'E-092'

# def export_all_freqs_and_AP_parameters(cell_ID, lookup_table):
    
# plotting specifications    
darkmode_bool = False
colors_dict = get_colors(darkmode_bool)

# table = pd.read_excel(directories.table_dir + 'InVitro_Database_copy.xlsx',
#                       sheet_name="PGFs",
#                       index_col='cell_ID')

frequencies = list(cc_APs_parameters.keys())

for frequency in frequencies:
    
    # PGF to load
    PGF = 'cc_APs_' + frequency
    PGF_parameters = cc_APs_parameters[frequency]
    
    # lookup_table = table.query(f'{PGF}.notnull()')
    
    
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
    
    t_concat = calc_time_series(v_concat, SR)
    
    # plt.plot(v_concat[0:15000])
    # plt.show()
    
    
    # filter voltage (to vf)
    vf = butter_filter(v_concat, 
                       order = 3,
                       cutoff = 1e3,
                       sampling_rate = SR)
    
    
    # limit the array to stimulus time frame + set time post stimulus to accomodate late APs
    t_post_stim = parameters.cc_APs_t_post_stim
    
    # use concatenated voltage trace
        # some protocols have a ISI shorter than 5 ms
        # the time used to accomodate late APs
        # also to filter concatenated trace
    
    # initialise empty list for indices for each step
    idc_stim_ls = [None] * n_steps
    
    # time of each step in ms
    t_step = sum(PGF_parameters.values())
    
    # number of points in pre & stim
    n_points_step = t_step * SR_ms
    n_points_pre = PGF_parameters['t_pre']*SR_ms
    n_points_stim = PGF_parameters['t_stim']*SR_ms
    n_points_post_stim = t_post_stim * SR_ms
    
    # # number of points in entire PGF
    # n_points_total = n_steps * t_step * SR_ms
    
    # create array with indices for stim and time after
    for idx in np.arange(n_steps):
        start_idx = int(n_points_step * idx + n_points_pre)
        stop_idx = int(start_idx + n_points_stim + n_points_post_stim)
        
        # conditional statement since last stimulus cannot be fully accommodated with post stim t
        fill_bool = False
        
        if stop_idx >= n_points:
            fill_bool = True
            n_points_to_fill = stop_idx - n_points
            stop_idx = n_points
            
    
        idc_stim_ls[idx] = np.arange(start = start_idx, 
                                     stop = stop_idx,
                                     dtype = int)
    
    
    # for step_idx in np.arange(n_steps):
    #     plt.plot(vf[idc_stim_ls[step_idx]])
    #     plt.ylim([-100, 50])
    #     plt.show()
    
    
    
    # %%
    
    # create array of the limited voltage steps
    v_ar = np.empty([n_steps, int(n_points_stim + n_points_post_stim)])
    t_ar = np.empty([n_steps, int(n_points_stim + n_points_post_stim)])
    
    for step_idx in np.arange(n_steps):
        
        if step_idx == (n_steps-1) and fill_bool:
            nan_fill = np.empty(n_points_to_fill)
            nan_fill[:] = np.nan
            v_last = vf[idc_stim_ls[step_idx]]
            v_last = np.append(v_last, nan_fill)
            v_ar[step_idx] = v_last
            t_ar[step_idx] = calc_time_series(v_last, SR)
        else:
            v_ar[step_idx] = vf[idc_stim_ls[step_idx]]
            t_ar[step_idx] = calc_time_series(vf[idc_stim_ls[step_idx]], SR)
        
    
    # %%
    
    
    # initialize lists
    t_spikes = []
    idx_spikes = []
    
    AP_all_params = pd.DataFrame(columns = ['v_peaks',
                                            't_peaks',
                                            'v_threshold',
                                            't_threshold',
                                            'idx_threshold',
                                            'v_amplitude',
                                            't_toPeak',
                                            'v_AHP',
                                            't_AHP',
                                            'idx_AHP',
                                            'v_AHP_amplitude',
                                            't_to_AHP',
                                            't_rise',
                                            'FWHM',
                                            'v_HM',
                                            't1_HM',
                                            't2_HM'])
    
    for step_idx in np.arange(n_steps):
        
        # limit data array to just the stimulus
        vs = v_ar[step_idx]
        # vl = vf
        
        # find peaks as spikes
        idx_peaks, dict_peak = sc.signal.find_peaks(vs, 
                                                    prominence = parameters.min_peak_prominence, 
                                                    distance = parameters.min_peak_distance * (SR_ms))
        
        idx_spikes.append(idx_peaks)
 
        ### opt: PLOT ###
        # plt.plot(vs)
        # plt.eventplot(idx_peaks, color = 'r', lineoffsets=30, linelengths=5)
        # plt.ylim([-100, 50])
        # plt.show()    
 
        # get times of spikes
        t_peaks = np.divide(idx_peaks, (SR_ms))
        t_spikes.append(t_peaks)
        
        AP_params = get_AP_parameters(vs, idx_peaks, 
                                      SR = SR,
                                      dvdt_threshold = 20,
                                      t_pre = 2,
                                      t_post = 10)
        
        AP_params['idx_step'] = step_idx
           
        AP_all_params = pd.concat([AP_all_params, AP_params])
    
    
    AP_all_params = AP_all_params.set_index('idx_step', drop = True, verify_integrity = True)
    
    
    
    # %% verification plot
    
    # opt: verification plot #
    
    n_rows = 10
    n_cols = 10
    
    darkmode_bool = True
    color_dict = get_colors(darkmode_bool)
    
    plt_idc = []
    
    # construct indices array for all plots in subplots
    for row in np.arange(n_rows):
        for col in np.arange(n_cols):
            plt_idc.append((row, col))
    
    
    
    fig_steps, axs_steps = plt.subplots(n_rows, n_cols,
                                        layout = 'constrained',
                                        sharey='all', sharex='all',
                                        dpi = 600,
                                        figsize = get_figure_size())
    
    fig_steps.set_constrained_layout_pads(wspace=0.05, w_pad=0.0,
                                          hspace=0.05, h_pad=0.0) 
    
    
    for step_idx, v_step in enumerate(v_ar):
        
        cur_t_step = t_ar[step_idx]
        cur_ax = axs_steps[plt_idc[step_idx][0]][plt_idc[step_idx][1]]
        
        if len(t_spikes[step_idx]):
            color = 'c'
            spike_bool = True
        else:
            color = 'm'
            spike_bool = False
       
        
        # if spike_bool:
        #     cur_ax.scatter(AP_all_params.at[step_idx, 't_threshold'], 
        #                     AP_all_params.at[step_idx, 'v_threshold'],
        #                     marker = 'x', 
        #                     c = 'r',
        #                     linewidth = 1,
        #                     s = 5)
    
        #     cur_ax.vlines(t_spikes[step_idx], 
        #                   AP_all_params.at[step_idx, 'v_threshold'], 
        #                   AP_all_params.at[step_idx, 'v_peaks'], 
        #                   color = 'r',
        #                   linewidth = 1)
    
            # cur_ax.hlines(AP_all_params.at[step_idx, 'v_HM'], 
            #               AP_all_params.at[step_idx, 't1_HM'], 
            #               AP_all_params.at[step_idx, 't2_HM'], 
            #               color = 'r',
            #               linewidth = 1)
            
            # cur_ax.hlines(40, 
            #               AP_all_params.at[step_idx, 't1_HM'], 
            #               AP_all_params.at[step_idx, 't2_HM'], 
            #               color = 'r',
            #               linewidth = 1)
    
        cur_ax.plot(cur_t_step, v_step, c = color, linewidth = 2)
    
    n_points_stim_n_post = int(n_points_stim + n_points_post_stim)
    t_stim_n_post = n_points_stim_n_post / SR_ms
    
    axs_steps[-1][-1].set_xlim([0,t_stim_n_post])
    axs_steps[-1][-1].set_xticks([])
    
    axs_steps[-1][-1].set_ylim([-100, 50])
    axs_steps[-1][-1].set_yticks([-100, 50])
    
    axs_steps[-1][-1].set_xticklabels([])
    axs_steps[-1][-1].set_yticklabels([])
    
    for row in np.arange(n_rows):
        for col in np.arange(n_cols):
            axs_steps[row][col].tick_params(axis = 'y', size = 0)
            axs_steps[row][col].tick_params(axis = 'x', size = 0)
    
    # path to save verification plot
    vplot_path = os.path.join(directories.vplot_dir, 'APs', cell_ID)
    
    if not os.path.exists(vplot_path):
        os.mkdir(vplot_path) 
        
    save_figures(fig_steps, f'All_AP_plot-{cell_ID}_{frequency}', vplot_path, darkmode_bool)
    
    plt.show()
    
    


# %% export values to excel files
       
    cell_path = os.path.join(directories.quant_data_dir, 'APs', cell_ID)
    
    if not os.path.exists(cell_path):
        os.mkdir(cell_path)
    
    table_path = os.path.join(cell_path, f'{cell_ID}_{frequency}.xlsx')

    AP_all_params.to_excel(table_path)










