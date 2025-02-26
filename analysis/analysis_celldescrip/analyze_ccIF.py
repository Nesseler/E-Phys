# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:16:29 2024

@author: nesseler
"""
import warnings
warnings.warn('Script not up-to-date!')
# %%

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sbn
import numpy as np
import scipy as sc

# custom directories & parameters
from parameters.directories_win import table_file, quant_data_dir, cell_descrip_dir, vplot_dir
from parameters.parameters import min_peak_prominence_ccIF, min_peak_distance_ccIF, min_max_peak_width_ccIF, dvdt_threshold, AP_parameters, t_expo_fit, popt_guess, r_squared_thresh, n_APs_initial_inst_freq
from parameters.PGFs import cc_IF_parameters
from functions.get_cell_IDs import get_cell_IDs_one_protocol, get_cell_IDs_all_ccAPfreqs

from functions.functions_constructors import construct_current_array
# from functions.functions_ccIF import get_IF_data
from functions.functions_import import get_traceIndex_n_file, get_cc_data
from functions.functions_useful import calc_time_series, butter_filter, calc_dvdt, calc_dvdt_padded, round_to_base, exp_func, calc_rsquared_from_exp_fit
from functions.functions_plotting import get_colors, get_figure_size, save_figures, set_font_sizes
from functions.functions_extractspike import get_AP_parameters

from analysis.analyze_ccsag_for_rinput_taumem import get_rinput_n_taumem_from_cc_sag

# %% settings

vplot_bool = True

darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)


# %% load data

# protocol 
PGF = 'cc_IF'

# get cell IDs
cell_IDs = get_cell_IDs_one_protocol(PGF)
# cell_IDs = get_cell_IDs_all_ccAPfreqs()

# get hold current as table
I_hold_table = pd.read_excel(table_file, sheet_name="V_or_I_hold", index_col='cell_ID').loc[cell_IDs, :]

# create dataframe for firing frequency
IF_df = pd.DataFrame(columns=cell_IDs, index = np.arange(-100, 400 + 1, 5))
IF_inst_df = pd.DataFrame(columns=cell_IDs, index = np.arange(-100, 400 + 1, 5))
IF_inst_initial_df = pd.DataFrame(columns=cell_IDs, index = np.arange(-100, 400 + 1, 5))

# create dataframe for other parameters
active_properties_df = pd.DataFrame(columns=['rheobase_abs', 'rheobase_rel', 'v_thres_rheobase_spike'], index = cell_IDs)
passiv_properties_df = pd.DataFrame(columns=['r_input', 'tau_mem'], index = cell_IDs)
fstAP_df = pd.DataFrame(columns = AP_parameters, index = cell_IDs)
step_idx_df = pd.DataFrame(columns=['rheobase_step_idx', 'maxfreq_step_idx', 'maxinstfreq_step_idx', 'maxinstinitialfreq_step_idx'], index = cell_IDs)


cells_todrop = []

# test cell 
# cell_IDs = ['E-185']

cell_IDs.reverse()

# cell_IDs = ['E-082', 'E-137', 'E-138', 'E-140']

for cell_ID in cell_IDs:
    
    print(f'Started: {cell_ID}')
    
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)
    
    # get IF data form file
    i, v, t, SR, n_steps = get_cc_data(file_path, traceIndex, 'ms')
    
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

    ### construct current dataframe
    i_hold = I_hold_table.at[cell_ID, PGF]
    
    # calculate current steps relative to I_hold
    ## rounded to nearest 5
    i_hold_rounded = round_to_base(i_hold, 5)
    
    # get current arrays and list of input current relative to i_hold
    i, i_input = construct_current_array(i_hold = i_hold_rounded,
                                         n_steps = n_steps,
                                         parameters_dict = cc_IF_parameters,
                                         SR_ms = SR_ms)
    
    # split concatenate arrays back to steps wise 
    # needs to occurs after filtering because of the filtering artifact
    
    v = [None] * n_steps
    t = [None] * n_steps
    peak_idc = [None] * n_steps
    idx_peaks_s = [None] * n_steps
    
    step_dur = cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim'] + cc_IF_parameters['t_post']
    step_points = step_dur * SR_ms
    
    pre_points = int(cc_IF_parameters['t_pre'] * SR_ms)
    pre_n_stim_points = int((cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim']) * SR_ms)
    
    
    # loop through steps to limit voltage trace
    for step_idx in np.arange(0, n_steps, 1):
        
        # calc start and stop indices for step
        start_idx = int(step_points * step_idx)
        stop_idx = int(start_idx + step_points)
        
        # set voltage trace of step
        v_step = vf[start_idx:stop_idx]
        
        # assign voltage and time traces to arrays
        v[step_idx] = v_step
        t[step_idx] = t_ms[start_idx:stop_idx]
        i_input_step = i_input[step_idx]
    
        # find peaks
        idc_peaks, dict_peak = sc.signal.find_peaks(v_step, 
                                                    prominence = min_peak_prominence_ccIF, 
                                                    distance = min_peak_distance_ccIF * (SR_ms),
                                                    width = np.multiply(min_max_peak_width_ccIF, SR_ms))

            
        # limit spike indices to stimulation time period
        idc_peaks = [idx_peak for idx_peak in idc_peaks if idx_peak > pre_points and idx_peak <= pre_n_stim_points]
        
        # write peak indices to array for all steps
        peak_idc[step_idx] = idc_peaks
        
        # get times of spikes
        t_spikes = np.divide(idc_peaks, SR_ms)
        
        # verification plot for detection of spikes in each step
        if False:
            plt.plot(v_step, linewidth = 1, c = colors_dict['primecolor'])
            plt.eventplot(idc_peaks, lineoffsets=60, colors = 'r', linelengths=5)
            plt.title(f'{cell_ID}: step #: {step_idx}')
            plt.ylim([-140, 75])
            plt.grid(False)
            plt.show()
            
        # calc frequency over entire step 
        n_spikes = len(idc_peaks)
        
        # write to dataframe
        IF_df.at[i_input_step, cell_ID] = n_spikes / (cc_IF_parameters['t_stim'] / 1e3)
    
        # calculate ISI
        if n_spikes >= 2:
            
            # get ISI as difference of spike times
            ISIs = np.diff(t_spikes)
            
            # calculate average ISI
            mean_ISI = np.mean(ISIs)
    
            # calculate the instantaneous firing frequency as inverse of the firing frequency
            inst_freq = (1 / mean_ISI ) * 1e3
            
            # write to dataframe
            IF_inst_df.at[i_input_step, cell_ID] = inst_freq
            
            # calc instant firing rate with first # APs
            if n_spikes >= n_APs_initial_inst_freq:
                # get only first three ISIs
                ISIs_subset = ISIs[:3]
                
                # calculate average ISI
                mean_ISI_subset = np.mean(ISIs_subset)
            
                # calculate the instantaneous firing frequency as inverse of the firing frequency
                initial_inst_freq = (1 / mean_ISI_subset ) * 1e3
            
                # write to dataframe
                IF_inst_initial_df.at[i_input_step, cell_ID] = initial_inst_freq
            
    ### rheobase ###
    # get first index in number of spikes where there more than 0 spikes
    rheobase_idx = next(idx for idx, n_spike in enumerate(IF_df[cell_ID].dropna()) if n_spike > 0)
    
    # get rheobase relative to holding
    rheobase_rel = i_input[rheobase_idx]
    
    # add holding current
    rheobase = rheobase_rel + i_hold_rounded

    # populate dataframe
    active_properties_df.at[cell_ID, 'rheobase_abs'] = rheobase
    active_properties_df.at[cell_ID, 'rheobase_rel'] = rheobase_rel
    step_idx_df.at[cell_ID, 'rheobase_step_idx'] = rheobase_idx
    
    # rheobase as voltage at threshold of first spike
    # get voltage trace with rheobase step
    v_rheobase = v[rheobase_idx]
    t_rheobase = t[0]
    dvdt_rheobase = calc_dvdt_padded(v_rheobase, t_rheobase)
    
    # get index of first spike
    idc_rheobase_spikes = peak_idc[rheobase_idx]
    
    # get number of spikes at rheobase
    n_rheobase_spikes = len(idc_rheobase_spikes)
    
    # get only first indec
    if n_rheobase_spikes > 1:
        idx_rheobase_spike = [idc_rheobase_spikes[0]]
    elif n_rheobase_spikes == 1:
        idx_rheobase_spike = idc_rheobase_spikes
    else:
        raise ValueError('size of list for rheobase spike')


    # rheobase step verification plot
    if vplot_bool:
        set_font_sizes()
        fig_rheobase = plt.figure()
        plt.plot(v_rheobase, linewidth = 1, c = colors_dict['primecolor'])
        plt.eventplot(idx_rheobase_spike, lineoffsets=60, colors = 'r', linelengths=5)
        plt.title(f'{cell_ID}: rheobase step #: {rheobase_idx}')
        plt.ylim([-140, 75])
        plt.grid(False)
        plt.show()
        vplots_path_rheobase = os.path.join(vplot_dir, 'cc_IF', 'rheobase')
        save_figures(fig_rheobase, f'{cell_ID}-{PGF}-rheobase_step', vplots_path_rheobase, darkmode_bool,
                     figure_format = 'png')
        
    
    # get AP parameters of first spike
    rheobase_spike_params, rheobase_spike_v = get_AP_parameters(t_rheobase, v_rheobase, dvdt_rheobase, idx_rheobase_spike)
    
    # write to active properties dataframe
    active_properties_df.at[cell_ID, 'v_thres_rheobase_spike'] = rheobase_spike_params.at[0, 'v_threshold']
    active_properties_df.at[cell_ID, 'n_rheobasespikes'] = n_rheobase_spikes
    
    # concatenate all AP parameters of first spike to dataframe
    fstAP_df.loc[cell_ID] = rheobase_spike_params.iloc[0]
    fstAP_df.at[cell_ID, 'SR_ms'] = SR_ms
    
    # calc time to rheobase spike
    t_rheobasespike = rheobase_spike_params['t_peaks']
    t_torheobasespike = t_rheobasespike - (cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim'])
    
    
    
    # rheobase plot with dvdt
    
    if vplot_bool:
        # calc dvdt for rheobase step
        fig_rheo_2, axs_rheo_2 = plt.subplots(nrows = 1, 
                                            ncols = 2, 
                                            layout = 'constrained',
                                            figsize = get_figure_size(),
                                            width_ratios= [2, 2]
                                            )
        
        set_font_sizes()
        
        fig_rheo_2.suptitle(f'{cell_ID} rheobase')
        
        v_range = [-100, 70]
        
        # voltage v time
        axs_rheo_2[0].plot(t[0], v_rheobase, lw = 1, c = colors_dict['primecolor'])
        # axs_rheo_2[0].eventplot(t_rheobasespike, lineoffsets=60, colors = 'r', linelengths=5)
        
        #x
        axs_rheo_2[0].set_xlabel('Time [ms]')
        axs_rheo_2[0].set_xlim([0, 1500]) 
        axs_rheo_2[0].set_xticks(np.arange(0, 1500 + 1, 250))
        axs_rheo_2[0].set_xticks(np.arange(0, 1500 + 1, 50), minor = True)    
        
        #y
        axs_rheo_2[0].set_ylabel('Voltage [mV]')
        axs_rheo_2[0].set_ylim(v_range) 
        axs_rheo_2[0].set_yticks(np.arange(v_range[0], v_range[1] + 1, 20))
        axs_rheo_2[0].set_yticks(np.arange(v_range[0], v_range[1] + 1, 5), minor = True)
        
        
        # inset marker
        box_tpad_pre = rheobase_spike_params.at[0, 'FWHM'] * 3
        box_tpad_post = rheobase_spike_params.at[0, 'FWHM'] * 20
        box_vpad = rheobase_spike_params.at[0, 'v_amplitude'] * 0.2
        box_ymin = rheobase_spike_params.at[0, 'v_AHP'] - box_vpad
        box_xmin = rheobase_spike_params.at[0, 't_peaks'] - box_tpad_pre
        box_height = rheobase_spike_params.at[0, 'v_amplitude'] + (box_vpad*2)
        box_width = rheobase_spike_params.at[0, 'FWHM'] + box_tpad_pre + box_tpad_post
        
        axs_rheo_2[0].add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                                          width = box_width, 
                                          height = box_height,
                                          fill = False,
                                          color = colors_dict['primecolor'],
                                          linestyle = '--'))
        
        
        # add inset
        ## ([left, bottom, width, height]), percentages
        ax_inset = fig_rheo_2.add_axes([0.325, 0.4, 0.15, 0.50])
        
        # voltage trace of spike marked within box
        idc_fstAP_rheobase_v = np.arange(start = idx_rheobase_spike[0] - int(box_tpad_pre * SR_ms), stop = idx_rheobase_spike[0] + int(box_tpad_post * SR_ms))
        
        # plot marled voltage trace
        ax_inset.plot(t_rheobase[idc_fstAP_rheobase_v], 
                      v_rheobase[idc_fstAP_rheobase_v],
                      c = colors_dict['primecolor'],
                      lw = 1)
        
        # plot extracted spike shape 
        t_extracted_spike = calc_time_series(rheobase_spike_v, SR)
        t_extracted_spike = t_extracted_spike + rheobase_spike_params.at[0, 't_threshold']
        ax_inset.plot(t_extracted_spike, 
                      rheobase_spike_v,
                      c = colors_dict['color2'],
                      lw = 1)
    
        
        # x
        ax_inset.set_xticks(ticks = np.arange(0, 1500 + 1, 250), labels = [])
        ax_inset.set_xticks(ticks = np.arange(0, 1500 + 1, 50), labels = [], minor = True)
        ax_inset.set_xlim([box_xmin, box_xmin + box_width])
        # y
        ax_inset.set_yticks(ticks = np.arange(v_range[0], v_range[1] + 1, 20), labels = [])
        ax_inset.set_yticks(ticks = np.arange(v_range[0], v_range[1] + 1, 5), labels = [], minor = True)
        ax_inset.set_ylim([box_ymin, box_ymin + box_height])
        
        # dvdt v voltage
        axs_rheo_2[1].plot(v_rheobase, dvdt_rheobase, lw = 1, c = colors_dict['primecolor'])
        
        # plot dvdt v v of rheobase first spike
        axs_rheo_2[1].plot(rheobase_spike_v, calc_dvdt_padded(rheobase_spike_v, calc_time_series(rheobase_spike_v, SR)), 
                           lw = 1, 
                           c = colors_dict['color2'])
           
        # axs_rheo_2[1].set_box_aspect(1)
        #x
        axs_rheo_2[1].set_xlabel('Voltage [mV]')
        axs_rheo_2[1].set_xlim(v_range)
        axs_rheo_2[1].set_xticks(np.arange(v_range[0], v_range[1] + 1, 20))
        axs_rheo_2[1].set_xticks(np.arange(v_range[0], v_range[1] + 1, 5), minor = True)
        #y
        axs_rheo_2[1].set_ylabel('Rate of membrane potential change [mV/ms]')
        axs_rheo_2[1].set_ylim([-150, 250]) 
        axs_rheo_2[1].set_yticks(np.arange(-150, 250 + 1, 50))
        axs_rheo_2[1].set_yticks(np.arange(-150, 250 + 1, 10), minor = True)
        
        plt.show()
    
        vplots_path_rheobase_fstAP = os.path.join(vplot_dir, 'cc_IF', 'rheobase_1stAP')
        save_figures(fig_rheo_2, f'{cell_ID}-{PGF}-rheobase_1stAP', vplots_path_rheobase_fstAP, darkmode_bool,
                     figure_format = 'png')

    # %% max freq steps
    
    # beware of the difference between the index in the dataframe and the 
    # step index
    
    # get step indices of max firing frequencies
    
    # max frequency (number of APs)
    # get max freq    
    maxfreq = IF_df[cell_ID].max()
    # get index of max freq in df & step index
    maxfreq_step_idx = np.nanargmax(IF_df[cell_ID].dropna().to_numpy())
    # get input current for max freq
    maxfreq_iinput = IF_df.index.to_list()[np.nanargmax(IF_df[cell_ID].to_numpy())]
    
    # max inst. frequency
    # get max freq
    maxinstfreq = IF_inst_df[cell_ID].max()
    # get index of max freq in df
    maxinstfreq_df_idx = np.nanargmax(IF_inst_df[cell_ID].dropna().to_numpy())
    # get input current for max freq
    maxinstfreq_iinput = IF_inst_df.index.to_list()[np.nanargmax(IF_inst_df[cell_ID].to_numpy())]    
    # get step index
    maxinstfreq_step_idx = IF_df[cell_ID].dropna().index.to_list().index(maxinstfreq_iinput)
     
    # max inst. inital frequency
    # get max freq
    maxinstinitialfreq = IF_inst_initial_df[cell_ID].max()
    # get index of max freq in df
    maxinstinitialfreq_df_idx = np.nanargmax(IF_inst_initial_df[cell_ID].dropna().to_numpy())
    # get input current for max freq
    maxinstinitialfreq_iinput = IF_inst_initial_df.index.to_list()[np.nanargmax(IF_inst_initial_df[cell_ID].to_numpy())]    
    # get step index
    maxinstinitialfreq_step_idx = IF_df[cell_ID].dropna().index.to_list().index(maxinstinitialfreq_iinput)
    
    # write to dataframe
    step_idx_df.at[cell_ID, 'maxfreq_step_idx'] = maxfreq_step_idx
    step_idx_df.at[cell_ID, 'maxinstfreq_step_idx'] = maxinstfreq_step_idx
    step_idx_df.at[cell_ID, 'maxinstinitialfreq_step_idx'] = maxinstinitialfreq_step_idx
    
    if vplot_bool:
        fig_max, axs_max = plt.subplots(nrows = 3,
                                        ncols = 2,
                                        layout = 'constrained',
                                        figsize = get_figure_size(),
                                        width_ratios = [4,1],
                                        sharex = 'col',
                                        sharey = 'col'
                                        )
        
        set_font_sizes()
        
        fig_max.suptitle(f'{cell_ID} max firing frequencies')
        
        
        # max freq
        axs_max[0][0].plot(t[0], 
                           v[maxfreq_step_idx], 
                           lw = 1, 
                           c = colors_dict['primecolor'])
    
        axs_max[0][1].plot(v[maxfreq_step_idx], 
                           calc_dvdt_padded(v[maxfreq_step_idx], t[0]), 
                           lw = 1, 
                           c = colors_dict['primecolor'])
        
        axs_max[0][0].text(x = 50,
                           y = 50,
                           s = f'{maxfreq_iinput} pA\n{maxfreq} Hz',
                           ha = 'left', va = 'top',
                           size = 9)
        
        # max inst freq
        axs_max[1][0].plot(t[0], v[maxinstfreq_step_idx], lw = 1, c = colors_dict['primecolor'])
    
        axs_max[1][1].plot(v[maxinstfreq_step_idx], calc_dvdt_padded(v[maxinstfreq_step_idx], t[0]) , lw = 1, c = colors_dict['primecolor'])

        axs_max[1][0].text(x = 50,
                           y = 50,
                           s = f'{maxinstfreq_iinput} pA\n{round(maxinstfreq, 2)} Hz',
                           ha = 'left', va = 'top',
                           size = 9)

        # max inst initial freq
        axs_max[2][0].plot(t[0], v[maxinstinitialfreq_step_idx], lw = 1, c = colors_dict['primecolor'])
    
        axs_max[2][1].plot(v[maxinstinitialfreq_step_idx], calc_dvdt_padded(v[maxinstinitialfreq_step_idx], t[0]) , lw = 1, c = colors_dict['primecolor'])

        axs_max[2][0].text(x = 50,
                           y = 50,
                           s = f'{maxinstinitialfreq_iinput} pA\n{round(maxinstinitialfreq, 2)} Hz',
                           ha = 'left', va = 'top',
                           size = 9)

        #x v
        axs_max[-1][0].set_xlabel('Time [ms]')
        axs_max[-1][0].set_xlim([0, 1500])
        axs_max[-1][0].set_xticks(np.arange(0, 1500 + 1, 250))
        axs_max[-1][0].set_xticks(np.arange(0, 1500 + 1, 50), minor = True)

        #y v
        axs_max[1][0].set_ylabel('Voltage [mV]')
        axs_max[-1][0].set_ylim(v_range)
        axs_max[-1][0].set_yticks(np.arange(v_range[0], v_range[1] + 1, 50))
        axs_max[-1][0].set_yticks(np.arange(v_range[0], v_range[1] + 1, 10), minor = True)

        #x dvdt
        axs_max[-1][1].set_xlabel('Voltage [mV]')
        axs_max[-1][1].set_xlim(v_range)
        axs_max[-1][1].set_xticks(np.arange(v_range[0], v_range[1] + 1, 50))
        axs_max[-1][1].set_xticks(np.arange(v_range[0], v_range[1] + 1, 10), minor = True)
        #y dvdt
        axs_max[1][1].set_ylabel('Rate of membrane potential change [mV/ms]')
        axs_max[-1][1].set_ylim([-150, 250]) 
        axs_max[-1][1].set_yticks(np.arange(-100, 250 + 1, 100))
        axs_max[-1][1].set_yticks(np.arange(-150, 250 + 1, 50), minor = True)

        vplots_path_maxfreq = os.path.join(vplot_dir, 'cc_IF', 'max_freq')
        save_figures(fig_max, f'{cell_ID}-{PGF}-max_freq', vplots_path_maxfreq, darkmode_bool,
                     figure_format = 'png')
        
    
    # %%


    ### tau_mem & R_input ###
    useful_steps_bool = True
    useful_steps = 0
    
    # create list of indices for stimulation period
    idc_stim = np.arange(pre_points, pre_n_stim_points)
    
    # create list of indices for exponential fit window
    delta_points = int(t_expo_fit * SR_ms)
    delta_points_pre = int(50 * SR_ms)
    idc_expoFit = np.arange(pre_points-1, pre_points -1 + delta_points)
    
    # include time periode before step to compare volatges for R_input calculation
    idc_pre = np.arange(pre_points - delta_points_pre, pre_points)
    idc_post = np.arange(pre_points, pre_points + delta_points)
    idc_withbuffer = np.arange(pre_points - delta_points_pre, pre_points + delta_points)
    x_withbuffer = np.arange(-delta_points_pre, +delta_points)
    
    # initialise dataframe for r_input calculation
    r_input_calc_df = pd.DataFrame(columns = ['step_idx', 'mean_v_pre', 'v_post_fitted', 'r_squared', 'delta_v', 'delta_i', 'r_input'])
    tau_mem_calc_df = pd.DataFrame(columns = ['step_idx', 'delta_v', 'delta_v_63', 'v_tau', 'tau_mem'])
    
    
    #vplot
    if vplot_bool:
        fig_expfit, axs_expfit = plt.subplots(nrows = 2,
                                              ncols = 3,
                                              layout = 'constrained',
                                              dpi = 600,
                                              sharex = True,
                                              sharey = True,
                                              figsize = get_figure_size())
        
        # melt array of axis for easier handling
        axs_expfit = axs_expfit.flatten()
        
        # set figure title
        fig_expfit.suptitle(cell_ID)
    
    # while loop until at least specified number of steps were analyized
    while useful_steps_bool:
        
        # loop through steps
        for step_idx in np.arange(0, n_steps, 1):
            
            i_input_step = i_input[step_idx]
            
            v_step_withbuffer = v[step_idx][idc_withbuffer]
            
            v_step_fit = v[step_idx][idc_expoFit]
            
            v_step_min = np.min(v_step_fit)
            
            v_step_min_idx = np.argmin(v_step_fit)
            
            # limit trace of first minimum
            v_step_expFit = v_step_fit[:v_step_min_idx]
            
            # calc x
            x_expFit = np.arange(len(v_step_expFit))

            # vplot
            if vplot_bool:
                # set step as subplot title
                axs_expfit[step_idx].set_title(f'Step #: {step_idx} {i_input_step} pA')
                
                # plot
                axs_expfit[step_idx].plot(x_withbuffer,
                                          v_step_withbuffer,
                                          color = 'gray',
                                          alpha = 0.5)
                axs_expfit[step_idx].hlines(y = v_step_min,
                                            xmin = -delta_points_pre,
                                            xmax = delta_points,
                                            color = 'r',
                                            linestyle = '--',
                                            alpha = 0.5)
                axs_expfit[step_idx].scatter(v_step_min_idx, v_step_min,
                                             c = 'r',
                                             marker = 'x',
                                             alpha = 0.5)

            try:
                #delta I
                delta_i = i_input[step_idx]
                
                if delta_i > -10:
                    raise ValueError('Zero input current')
                    print(f'{cell_ID} close to zero input current')
            
                
                # calc delta v with mini
                v_pre = v[step_idx][idc_pre]
                v_pre_mean = np.mean(v_pre)
                delta_v = abs(v_step_min - v_pre_mean)
                
                # fit exponential curve
                popt, pcov = sc.optimize.curve_fit(exp_func, x_expFit, v_step_expFit, p0 = [delta_v, *popt_guess[1:]], maxfev = 5000,
                                                   bounds = ([delta_v-3, 0, -200], [delta_v+3, 1, -80]))
                
                # popt = [delta_v, *popt]
                
                r_squared = calc_rsquared_from_exp_fit(x_expFit, v_step_expFit, popt)
                
                ### R_INPUT ###
                #R_input = ∆U/∆I
                #∆U = a = popt[0], ∆I = 20 (for first step)
                # v_pre = v[step_idx][idc_pre]
                # v_pre_mean = np.mean(v_pre)
                # delta_v = (popt[2] - v_pre_mean)
                delta_v = -popt[0]
                    
                # calc r_input for current step
                r_input = ( delta_v / delta_i ) * 1e3 #in MOhm
                
                r_input_calc_df.at[step_idx, 'step_idx'] = step_idx
                r_input_calc_df.at[step_idx, 'mean_v_pre'] = v_pre_mean
                r_input_calc_df.at[step_idx, 'v_post_fitted'] = popt[2]
                r_input_calc_df.at[step_idx, 'r_squared'] = r_squared
                r_input_calc_df.at[step_idx, 'delta_v'] = delta_v
                r_input_calc_df.at[step_idx, 'delta_i'] = delta_i
                r_input_calc_df.at[step_idx, 'r_input'] = r_input
                
                
                ### MEMBRANE TIME CONSTANT
                # tau_mem
                # time it takes the potential to reach 1 - (1/e) (~63%) of the max voltage
            
                # calc 1 - 1/e
                tau_perc_value = 1-(1/np.exp(1))
                
                # calc max voltage delta
                delta_v_63 = delta_v * tau_perc_value
                
                # calc 63 % of max voltage
                v_tau = v_pre_mean + delta_v_63
                
                # calc time (indices first) it takes to reach v_tau
                idx_63 = - (np.log((v_tau - popt[2]) / (popt[0]))) / (popt[1])
                tau_mem = idx_63 / (SR / 1e3) 
                
                # append to list
                tau_mem_calc_df.at[step_idx, 'step_idx'] = step_idx
                tau_mem_calc_df.at[step_idx, 'delta_v'] = delta_v
                tau_mem_calc_df.at[step_idx, 'delta_v_63'] = delta_v_63
                tau_mem_calc_df.at[step_idx, 'v_tau'] = v_tau
                tau_mem_calc_df.at[step_idx, 'tau_mem'] = tau_mem
                
     
                # vplot
                if vplot_bool:               
                    # plot exponential fit
                    axs_expfit[step_idx].plot(x_expFit,
                                              exp_func(x_expFit, *popt),
                                              linestyle = '--',
                                              color = colors_dict['color2'])
                    
                    axs_expfit[step_idx].set_ylim([-150, -50])
                    
                    # text field with fit info
                    axs_expfit[step_idx].text(x = delta_points-250,
                                              y = -145,
                                              s = f'{popt[0]}\n{popt[1]}\n{popt[2]}\nr^2: {r_squared}\nr_input: {r_input}\ntau_mem: {tau_mem}',
                                              va = 'bottom',
                                              ha = 'right',
                                              fontsize = 8)
                    
                if r_squared > r_squared_thresh:
                    useful_steps += 1
            
            except RuntimeError:
                print(f'{cell_ID} step number {step_idx+1} has been omitted')
            
            except ValueError as m:
                print(f'{cell_ID} {str(m)}')
                useful_steps_bool = False
                break
            
            if useful_steps == 3 or step_idx == 6:
                useful_steps_bool = False
                break


    r_input_calc_df = r_input_calc_df.set_index('step_idx')
    tau_mem_calc_df = tau_mem_calc_df.set_index('step_idx')
        
    # save excel sheet    
    cell_path = os.path.join(quant_data_dir, 'cc_IF', cell_ID)
    
    if not os.path.exists(cell_path):
        os.mkdir(cell_path)
        
    r_input_calc_path = os.path.join(cell_path, f'{cell_ID}-{PGF}-R_input_calc.xlsx')
    r_input_calc_df.to_excel(r_input_calc_path, index_label='step_idx')
    
    tau_mem_calc_path = os.path.join(cell_path, f'{cell_ID}-{PGF}-tau_mem_calc.xlsx')
    tau_mem_calc_df.to_excel(tau_mem_calc_path, index_label='step_idx')

              
    if vplot_bool:
        # show plot
        [ax.grid(False) for ax in axs_expfit]
        plt.show()
        
        # save figure
        vplots_path_passive_properties = os.path.join(vplot_dir, 'cc_IF', 'passive_properties')
        
        save_figures(fig_expfit, f'{cell_ID}-ccIF-expon_fit', vplots_path_passive_properties, darkmode_bool,
                     figure_format = 'png')
        
    
    ### break out if fitting to first steps of IF is not successful ###
    n_steps_w_good_fit = len(r_input_calc_df[r_input_calc_df['r_squared'] > r_squared_thresh])

    if n_steps_w_good_fit < 3:
        print(f'{cell_ID} will try to use cc_sag')
        
        try:
            # call sag_analysis function
            r_input, tau_mem = get_rinput_n_taumem_from_cc_sag(cell_ID, vplot_bool, darkmode_bool)
            
        except ValueError:
            print(f'{cell_ID} will be disgarded')
            cells_todrop.append(cell_ID)
    
    else:
        # calc values as mean of 3 steps
        r_input = r_input_calc_df['r_input'].mean()   
        tau_mem = tau_mem_calc_df['tau_mem'].mean()


    passiv_properties_df.at[cell_ID, 'r_input'] = r_input
    passiv_properties_df.at[cell_ID, 'tau_mem'] = tau_mem


# %%

# calculate the membrane capacitance as
# c_mem = tau_mem / r_input
# times 1e3 for conversion to pF
passiv_properties_df['c_mem'] = (passiv_properties_df['tau_mem'] / passiv_properties_df['r_input']) * 1e3


# %% drop cells with unsuccessful calculations

passiv_properties_df.drop(index=cells_todrop, inplace = True)
active_properties_df.drop(index=cells_todrop, inplace = True)
IF_df.drop(columns=cells_todrop, inplace = True)
IF_inst_df.drop(columns=cells_todrop, inplace = True)
IF_inst_initial_df.drop(columns=cells_todrop, inplace = True)
fstAP_df.drop(index=cells_todrop, inplace = True)

# %% add max frequencies

active_properties_df['max_freq'] = IF_df.max(axis = 0)
active_properties_df['max_inst_freq'] = IF_inst_df.max(axis = 0)
active_properties_df['max_inst_initial_freq'] = IF_inst_initial_df.max(axis = 0)

# %% save measurements to excel file

if False:

    passiv_properties_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_label = 'cell_ID')    
    active_properties_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_label = 'cell_ID')   
    IF_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_label = 'i_input')       
    IF_inst_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_label = 'i_input')
    IF_inst_initial_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-IF_inst_initial.xlsx'), index_label = 'i_input')  
    
    fstAP_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_label = 'cell_ID')
    
    step_idx_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-step_indices.xlsx'), index_label = 'cell_ID')

# %%



# plt.figure()
# for cell_ID in cell_IDs:
#     plt.plot(IF_df[cell_ID])
    
# plt.show












