# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 16:10:53 2024

@author: nesseler

This script loads all cc_APs_xHz procotol indexed in the InVitro_Database.xlsx
and writes all AP parameters that can be analyzed to the qData folders.

Option: set vplots_bool to True in order to plot the corresponding verification
plots and save them to the vplots directory.

"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

# custom directories & parameters
from parameters.directories_win import table_dir, raw_data_dir, vplot_dir, quant_data_dir, cell_descrip_dir
from parameters.PGFs import cc_APs_parameters, cc_APs_t_stims_df, cc_APs_total_dur
from parameters.parameters import cc_APs_t_post_stim, min_peak_prominence, min_peak_distance

# custom functions
from functions.functions_ccIF import get_IF_data
# from functions.functions_spiketrains import get_AP_parameters
from functions.functions_useful import calc_time_series, butter_filter, calc_dvdt_padded
from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size
from functions.functions_extractspike import get_AP_parameters


# %% verification plots option

vplots_bool = False

# %% load InVitro database and all PGF indices to be loaded

table = pd.read_excel(table_dir + 'InVitro_Database.xlsx',
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

# cell IDs 
cell_IDs = list(lookup_table.index)

# %% plotting specifications if vplots_bool is given
if vplots_bool:
    darkmode_bool = False
    colors_dict = get_colors(darkmode_bool)
    set_font_sizes()


# %% 

nAPs_df = pd.DataFrame(index = frequencies)
mean_ISIs_df = pd.DataFrame(index = frequencies)
resul_freq_df = pd.DataFrame(index = frequencies)

mean_FWHM_df = pd.DataFrame(index = frequencies)
mean_tpeaks_df = pd.DataFrame(index = frequencies)
mean_vamplitude_df = pd.DataFrame(index = frequencies)


# %% get all AP parameters for one cell

# test cell E-092
# cell_IDs = ['E-119', 'E-122' ,'E-123' ,'E-124' ,'E-125' , 'E-126', 'E-129']

for cell_ID in cell_IDs:

    # initialise dataframes of parameters to be summarized
    t_APs = pd.DataFrame()
    n_APs = pd.DataFrame()
    
    # loop through all frequencies to be analyzed
    for frequency in frequencies:
        
        # PGF to load
        PGF = 'cc_APs_' + frequency
        PGF_parameters = cc_APs_parameters[frequency]
    
        # get indices of current cell with the dataframe containing all indices    
        group_idx = int(lookup_table.at[cell_ID, 'group'])-1
        series_idx = int(lookup_table.at[cell_ID, f'{PGF}'])-1
        
        # construct traceIndex with indices
        traceIndex = [group_idx, series_idx, 0, 0]
        
        # call on data file with indices from dataframe above
        current_file = lookup_table.at[cell_ID, 'file']
        data_file_path = os.path.join(raw_data_dir, current_file + '.dat')  
        data_file_path_str = fr"{data_file_path}"
        
        # get IF data form file
        i, v, t, SR, n_steps = get_IF_data(data_file_path_str, traceIndex, 'ms')
        
        # sampling rate in ms
        SR_ms = SR / 1e3
        
        # concatenate individual steps
        n_points = int(np.shape(i)[0] * np.shape(i)[1])      
        v_concat = v.flatten()    
        t = calc_time_series(v_concat, SR)
        
        
        # filter voltage (to vf)
        vf = butter_filter(v_concat, 
                           order = 3,
                           cutoff = 1e3,
                           sampling_rate = SR)
        
        # calc dvdt (for concatenated trace)
        dvdt = calc_dvdt_padded(vf, t)
        
        
        # limit the array to stimulus time frame + set time post stimulus to accomodate late APs
        t_post_stim = cc_APs_t_post_stim
        
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
        
        # create array of the limited voltage steps
        v_ar = np.empty([n_steps, int(n_points_stim + n_points_post_stim)])
        t_ar = np.empty([n_steps, int(n_points_stim + n_points_post_stim)])
        dvdt_ar = np.empty([n_steps, int(n_points_stim + n_points_post_stim)])
        
        for step_idx in np.arange(n_steps):
            
            # for last stimulation (therefore v_last (step)) 
            # pad array with nan
            if step_idx == (n_steps-1) and fill_bool:
                nan_fill = np.empty(n_points_to_fill)
                nan_fill[:] = np.nan
                
                # pad v trace
                v_last = vf[idc_stim_ls[step_idx]]
                v_last = np.append(v_last, nan_fill)
                v_ar[step_idx] = v_last
                
                # pad dvdt trace
                dvdt_last = dvdt[idc_stim_ls[step_idx]]
                dvdt_last = np.append(dvdt_last, nan_fill)
                dvdt_ar[step_idx] = dvdt_last
                
                # calc time series
                t_ar[step_idx] = calc_time_series(v_last, SR)
                
            else:
                v_ar[step_idx] = vf[idc_stim_ls[step_idx]]
                t_ar[step_idx] = calc_time_series(vf[idc_stim_ls[step_idx]], SR)
                dvdt_ar[step_idx] = dvdt[idc_stim_ls[step_idx]]
        
        
        # loop through all steps
        # analyze spikes
        for step_idx in np.arange(n_steps):
            
            # limit data array to just the stimulus
            vs = v_ar[step_idx]
            ts = t_ar[step_idx]
            dvdts = dvdt_ar[step_idx]
          
            # find peaks as spikes
            idx_peaks, dict_peak = sc.signal.find_peaks(vs, 
                                                        prominence = min_peak_prominence, 
                                                        distance = min_peak_distance * (SR_ms))
     
            # exception handling: spike threshold not in analyized window
            # extract_spike functions throws a ValueError with specific error message
            # if this error message is caugth AP_params functions is run again
            # with empty peak index array              
            try:
                AP_params = get_AP_parameters(t_spiketrain = ts,
                                              v_spiketrain = vs,
                                              dvdt_spiketrain = dvdts,
                                              idc_spikes = idx_peaks,
                                              SR = SR)
                
                
            except ValueError as e:
                if str(e) == 'AP threshold not crossed':
                    # print(step_idx, 'caught')
                    
                    AP_params = get_AP_parameters(t_spiketrain = ts,
                                                  v_spiketrain = vs,
                                                  dvdt_spiketrain = dvdts,
                                                  idc_spikes = [],
                                                  SR = SR)
                        
                
            
            AP_params['idx_step'] = step_idx
            
            # for first step the dataframe with all AP parameters represents 'all' AP
            # afterwards the dataframes are concatenated
            if step_idx == 0:
                AP_all_params = AP_params
            else:
                AP_all_params = pd.concat([AP_all_params, AP_params])
                
            
        
        # set the step index as index of pandas array
        AP_all_params = AP_all_params.set_index('idx_step', drop = True, verify_integrity = True)
        
        # get time points of APs
        t_APs[frequency] = AP_all_params['t_peaks']
        
        # get number of APs for given stimulation frequency
        n_APs[frequency] = [len(AP_all_params.query('v_peaks.notnull()').index)]
        
        # calc average FWHM, t_peaks, and v_amplitude
        mean_FWHM_df.at[frequency, cell_ID] = AP_all_params['FWHM'].mean()
        mean_tpeaks_df.at[frequency, cell_ID] = AP_all_params['t_peaks'].mean()
        mean_vamplitude_df.at[frequency, cell_ID] = AP_all_params['v_amplitude'].mean()
        
        
        # opt: verification plot #
        if vplots_bool:  
            n_rows = 10
            n_cols = 10
            
            darkmode_bool = False
            color_dict = get_colors(darkmode_bool)
            
            plt_idc = []
            
            # construct indices array for all plots in subplots
            for row in np.arange(n_rows):
                for col in np.arange(n_cols):
                    plt_idc.append((row, col))
            
            
            
            fig_steps, axs_steps = plt.subplots(n_rows, n_cols,
                                                layout = 'constrained',
                                                sharey='all', sharex='all',
                                                dpi = 600)
            
            fig_steps.set_constrained_layout_pads(wspace=0.05, w_pad=0.0,
                                                  hspace=0.05, h_pad=0.0) 
            
            
            for step_idx, v_step in enumerate(v_ar):
                
                cur_t_step = t_ar[step_idx]
                cur_ax = axs_steps[plt_idc[step_idx][0]][plt_idc[step_idx][1]]
                
                if ~np.isnan(AP_all_params.at[step_idx, 'v_peaks']):
                    color = 'c'
                else:
                    color = 'm'
            
                cur_ax.plot(cur_t_step, v_step, c = color, linewidth = 1)
            
            n_points_stim_n_post = int(n_points_stim + n_points_post_stim)
            t_stim_n_post = n_points_stim_n_post / SR_ms
            
            # figure title
            fig_steps.suptitle(f'{cell_ID} {frequency}')        
            
            axs_steps[-1][-1].set_xlim([0,t_stim_n_post])
            axs_steps[-1][-1].set_xticks([])
            
            axs_steps[-1][-1].set_ylim([-100, 50])
            axs_steps[-1][-1].set_yticks([-100, -50, 0, 50])
            
            axs_steps[-1][-1].set_xticklabels([])
            axs_steps[-1][-1].set_yticklabels([])
            
            for row in np.arange(n_rows):
                for col in np.arange(n_cols):
                    axs_steps[row][col].tick_params(axis = 'y', size = 0)
                    axs_steps[row][col].tick_params(axis = 'x', size = 0)
            
            # path to save verification plot
            vplot_path = os.path.join(vplot_dir, 'APs', cell_ID)
            
            if not os.path.exists(vplot_path):
                os.mkdir(vplot_path) 
                
            save_figures(fig_steps, f'All_AP_plot-{cell_ID}_{frequency}', vplot_path, darkmode_bool)
            
            plt.show()
        
            
            # trace vplot
            nrows = 4
            fig_trace, axs_trace = plt.subplots(nrows = nrows,
                                                ncols = 1,
                                                sharex = True,
                                                sharey = True,
                                                layout = 'constrained',
                                                figsize = get_figure_size())
    
            # create row indices
            idc_row = np.arange(0, n_points + 1, n_points / nrows, dtype = int)
    
            for row in range(nrows):
                v_row = vf[idc_row[row]:idc_row[row+1]]
                t_row = calc_time_series(v_row, SR)
                axs_trace[row].plot(t_row, v_row)
                
            axs_trace[0].set_ylim([-100, 60])
            fig_trace.supylabel('Voltage [mV]')
            
            axs_trace[-1].set_xlabel('Time [ms]')
            axs_trace[-1].set_xlim([t_row[0], t_row[-1]])
            
            save_figures(fig_trace, f'Trace-{cell_ID}_{frequency}', vplot_path, darkmode_bool)
    
            plt.show()
    
        # export values to excel files
           
        cell_path = os.path.join(quant_data_dir, 'APs', cell_ID)
        
        if not os.path.exists(cell_path):
            os.mkdir(cell_path)
        
        table_path = os.path.join(cell_path, f'{cell_ID}_{frequency}.xlsx')
    
        AP_all_params.to_excel(table_path)
        
    
    # %% combine all frequencies for one cell
    
    nAPs_df[cell_ID] = n_APs.transpose()
    
    # add both dataframes and create spike times in continues timeseries
    t_APs_cont = t_APs.add(cc_APs_t_stims_df)
    
    # calc mean ISI from continues spike times
    mean_ISIs_df[cell_ID] = [np.mean(t_APs_cont[freq].dropna().diff()) for freq in frequencies]
    
      
    # %% convert mean ISI to frequency
    
    # ms to s
    resul_freq_df[cell_ID] = mean_ISIs_df[cell_ID].div(1000)
    
    # ISI to freq
    resul_freq_df[cell_ID] = resul_freq_df[cell_ID].rdiv(1)
    
    
    # %% verification plot
    
    if vplots_bool:
        ax_keys = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
        
        fig_event, axs_event = plt.subplot_mosaic('AGHI;BGHI;CGHI;DGHI;EGHI;FGHI', 
                                                  layout = 'constrained',
                                                  figsize = get_figure_size(),
                                                  gridspec_kw = {'width_ratios': [6,2,2,2]})
        
        
        for idx_freq, freq in enumerate(frequencies):
            axs_event[ax_keys[idx_freq]].eventplot(t_APs_cont[freq].dropna(), colors = colors_dict['primecolor'])
            axs_event[ax_keys[idx_freq]].set_xlim([0, cc_APs_total_dur[freq]])
            axs_event[ax_keys[idx_freq]].set_xticks(np.linspace(0, cc_APs_total_dur[freq], 6))
            axs_event[ax_keys[idx_freq]].set_ylabel(freq)
            
            # eventplot for stimulations
            axs_event[ax_keys[idx_freq]].eventplot(cc_APs_t_stims_df[freq] + 5, lineoffsets = 2.25, linelength = 0.75, colors = 'grey')
        
        fig_event.suptitle(cell_ID)
        
        [axs_event[ax_i].set_yticks([]) for ax_i in ax_keys[:6]]
        
        [axs_event[ax_i].spines[spine].set_visible(False) for ax_i in ax_keys[:6] for spine in ['left', 'top', 'right']]
        
        [axs_event[ax_i].grid(False) for ax_i in ax_keys[:]]
        
        fig_event.supxlabel('Time [ms]')
        
        
        # number of APs
        axs_event['G'].plot(nAPs_df[cell_ID], nAPs_df.index)
        
        axs_event['G'].set_xlim([0, 100])
        axs_event['G'].set_xlabel('Number of APs\n[#]')
        
        # mean ISI
        axs_event['H'].plot(mean_ISIs_df[cell_ID], mean_ISIs_df.index)
        
        axs_event['H'].set_xlim([0, 1000])
        axs_event['H'].set_xlabel('Mean ISI\n[ms]')
        
        
        # resulting frequency
        axs_event['I'].plot(resul_freq_df[cell_ID], resul_freq_df.index)
        
        axs_event['I'].set_xlim([0, 75])
        axs_event['I'].set_xlabel('Resulting freq.\n[Hz]')
        
        # save figure as verification plot
        fig_path = os.path.join(vplot_dir, 'APs', cell_ID)
        if not os.path.exists(fig_path):
            os.mkdir(fig_path)
        save_figures(fig_event, f'{cell_ID}_ccAPs', fig_path, darkmode_bool)
        
        plt.show()


    print(f'Finished: {cell_ID}')


# %%

mean_nAPs_df = pd.DataFrame({'mean_APs' : nAPs_df.mean(axis=0)})
std_nAPs_df = pd.DataFrame({'std_APs' : nAPs_df.std(axis=0)})


# %% write dataframes to cell description folder

# save measurements to excel file
nAPs_df.to_excel(os.path.join(cell_descrip_dir, 'ccAPs-nAPs.xlsx'), index_label = 'frequencies')
mean_ISIs_df.to_excel(os.path.join(cell_descrip_dir, 'ccAPs-mean_ISIs.xlsx'), index_label = 'frequencies')
resul_freq_df.to_excel(os.path.join(cell_descrip_dir, 'ccAPs-resul_freq.xlsx'), index_label = 'frequencies')

mean_FWHM_df.to_excel(os.path.join(cell_descrip_dir, 'ccAPs-mean_FWHM.xlsx'), index_label = 'frequencies')
mean_tpeaks_df.to_excel(os.path.join(cell_descrip_dir, 'ccAPs-mean_tpeaks.xlsx'), index_label = 'frequencies')
mean_vamplitude_df.to_excel(os.path.join(cell_descrip_dir, 'ccAPs-mean_vamplitude.xlsx'), index_label = 'frequencies')

mean_nAPs_df.to_excel(os.path.join(cell_descrip_dir, 'ccAPs-mean_nAPs.xlsx'), index_label = 'frequencies')
std_nAPs_df.to_excel(os.path.join(cell_descrip_dir, 'ccAPs-std_nAPs.xlsx'), index_label = 'frequencies')


    
    
    