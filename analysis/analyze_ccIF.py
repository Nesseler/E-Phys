# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:16:29 2024

@author: nesseler
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn
import numpy as np
import scipy as sc

# custom directories & parameters
from parameters.directories_win import table_file, quant_data_dir, cell_descrip_dir, vplot_dir
from parameters.parameters import min_peak_prominence_ccIF, min_peak_distance_ccIF, min_max_peak_width_ccIF, dvdt_threshold, AP_parameters, t_expo_fit, popt_guess, r_squared_thresh
from parameters.PGFs import cc_IF_parameters
from getter.get_cell_IDs import get_cell_IDs_one_protocol, get_cell_IDs_all_ccAPfreqs

from functions.functions_constructors import construct_current_array
from functions.functions_ccIF import get_IF_data
from functions.functions_import import get_traceIndex_n_file
from functions.functions_useful import calc_time_series, butter_filter, calc_dvdt, calc_dvdt_padded, round_to_base, exp_func, calc_rsquared_from_exp_fit
from functions.functions_plotting import get_colors, get_figure_size, save_figures, plot_t_vs_v
from functions.functions_extractspike import get_AP_parameters

from analysis.analyze_ccsag_for_rinput_taumem import get_rinput_n_taumem_from_cc_sag

# %% settings

vplot_bool = False

darkmode_bool = False
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

# create dataframe for other parameters
active_properties_df = pd.DataFrame(columns=['rheobase_abs', 'rheobase_rel', 'v_thres_rheobase_spike'], index = cell_IDs)
passiv_properties_df = pd.DataFrame(columns=['r_input', 'tau_mem'], index = cell_IDs)
fstAP_df = pd.DataFrame(columns = AP_parameters, index = cell_IDs)

# test cell 
# cell_IDs = ['E-111']

# cell_IDs = ['E-082', 'E-137', 'E-138', 'E-140']

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
    
    # rheobase as voltage at threshold of first spike
    # get voltage trace with rheobase step
    v_rheobase = v[rheobase_idx]
    t_rheobase = t[rheobase_idx]
    dvdt_rheobase = calc_dvdt_padded(v_rheobase, t_rheobase)
    
    # get index of first spike
    idc_rheobase_spikes = peak_idc[rheobase_idx]
    
    # get only first indec
    if len(idc_rheobase_spikes) > 1:
        idx_rheobase_spike = [idc_rheobase_spikes[0]]
    elif len(idc_rheobase_spikes) == 1:
        idx_rheobase_spike = idc_rheobase_spikes
    else:
        raise ValueError('size of list for rheobase spike')

    # rheobase step verification plot
    if vplot_bool:
        plt.plot(v_rheobase, linewidth = 1, c = colors_dict['primecolor'])
        plt.eventplot(idx_rheobase_spike, lineoffsets=60, colors = 'r', linelengths=5)
        plt.title(f'{cell_ID}: rheobase step #: {rheobase_idx}')
        plt.ylim([-140, 75])
        plt.grid(False)
        plt.show()
    
    # get AP parameters of first spike
    rheobase_spike_params, rheobase_spike_v = get_AP_parameters(t_rheobase, v_rheobase, dvdt_rheobase, idx_rheobase_spike)
    
    # write to active properties dataframe
    active_properties_df.at[cell_ID, 'v_thres_rheobase_spike'] = rheobase_spike_params.at[0, 'v_threshold']
    
    # concatenate all AP parameters of first spike to dataframe
    fstAP_df.loc[cell_ID] = rheobase_spike_params.iloc[0]
    fstAP_df.at[cell_ID, 'SR_ms'] = SR_ms


    ### tau_mem & R_input ###
    useful_steps_bool = True
    useful_steps = 0
    
    # create list of indices for stimulation period
    idc_stim = np.arange(pre_points, pre_n_stim_points)
    
    # create list of indices for exponential fit window
    delta_points = int(t_expo_fit * SR_ms)
    idc_expoFit = np.arange(pre_points, pre_points + delta_points)
    
    # include time periode before step to compare volatges for R_input calculation
    idc_pre = np.arange(pre_points - delta_points, pre_points)
    idc_post = np.arange(pre_points, pre_points + delta_points)
    idc_withbuffer = np.arange(pre_points - delta_points, pre_points + delta_points)
    x_withbuffer = np.arange(-delta_points, +delta_points)
    
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
                                            xmin = -delta_points,
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
            
                # calc x
                x_expFit = np.arange(len(v_step_expFit))
                
                # fit exponential curve
                popt, pcov = sc.optimize.curve_fit(exp_func, x_expFit, v_step_expFit, p0=popt_guess, maxfev = 1000)
                
                r_squared = calc_rsquared_from_exp_fit(x_expFit, v_step_expFit, popt)
                
                ### R_INPUT ###
                #R_input = ∆U/∆I
                #∆U = a = popt[0], ∆I = 20 (for first step)
                v_pre = v[step_idx][idc_pre]
                v_pre_mean = np.mean(v_pre)
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
                    
                    axs_expfit[step_idx].set_ylim([-140, -75])
                    
                    # text field with fit info
                    axs_expfit[step_idx].text(x = -delta_points+100,
                                              y = -140,
                                              s = f'{popt[0]}\n{popt[1]}\n{popt[2]}\nr^2: {r_squared}',
                                              va = 'bottom',
                                              ha = 'left',
                                              fontsize = 8)
                    
                    axs_expfit[step_idx].text(x = 0+100,
                                              y = -140,
                                              s = f'r_input: {r_input}\ntau_mem: {tau_mem}',
                                              va = 'bottom',
                                              ha = 'left',
                                              fontsize = 8)
                    
                if r_squared > r_squared_thresh:
                    useful_steps += 1
            
            except RuntimeError:
                print(f'{cell_ID} step number {step_idx+1} has been omitted')
            
            except ValueError:
                useful_steps_bool = False
                break
            
            if useful_steps == 3 or step_idx ==6:
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
        
        # save excel sheet    
        cell_vplots_path = os.path.join(vplot_dir, 'cc_IF', cell_ID)
        
        if not os.path.exists(cell_vplots_path):
            os.mkdir(cell_vplots_path)
        
        save_figures(fig_expfit, f'{cell_ID}-ccIF-expon_fit', cell_vplots_path, darkmode_bool)
        
    
    ### break out if fitting to first steps of IF is not successful ###
    n_steps_w_good_fit = len(r_input_calc_df[r_input_calc_df['r_squared'] > 0.75])

    if n_steps_w_good_fit < 3:
        print(f'{cell_ID} will use cc_sag')
        
        # call sag_analysis function
        r_input, tau_mem = get_rinput_n_taumem_from_cc_sag(cell_ID, vplot_bool, darkmode_bool)
    
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


# %%

# save measurements to excel file
passiv_properties_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_label = 'cell_ID')    
active_properties_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_label = 'cell_ID')   
IF_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_label = 'i_input')       
IF_inst_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_label = 'i_input')   
 
fstAP_df.to_excel(os.path.join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_label = 'cell_ID')

# %%



# plt.figure()
# for cell_ID in cell_IDs:
#     plt.plot(IF_df[cell_ID])
    
# plt.show














