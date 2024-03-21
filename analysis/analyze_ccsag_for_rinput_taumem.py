# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:55:49 2024

@author: nesseler
"""

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
from parameters.directories_win import table_file, quant_data_dir, vplot_dir
from parameters.parameters import t_expo_fit, popt_guess, r_squared_thresh
from parameters.PGFs import cc_sag_parameters

from functions.functions_constructors import construct_current_array
from functions.functions_ccIF import get_IF_data
from functions.functions_import import get_traceIndex_n_file
from functions.functions_useful import calc_time_series, butter_filter, round_to_base, exp_func, calc_rsquared_from_exp_fit
from functions.functions_plotting import get_colors, get_figure_size, save_figures

# cell_ID = 'E-138'


def get_rinput_n_taumem_from_cc_sag(cell_ID, vplot_bool, darkmode_bool):
    
    colors_dict, region_colors = get_colors(darkmode_bool)
    
    # protocol 
    PGF = 'cc_sag'
    
    # get hold current as table
    I_hold_table = pd.read_excel(table_file, sheet_name="V_or_I_hold", index_col='cell_ID').loc[cell_ID, :]
    
    
    
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)
    
    # get IF data form file
    i, v, t, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')
    
    # sampling rate in ms
    SR_ms = SR / 1e3
    
    # concatenate individual steps
    
    v_concat = v.flatten() 
    
    t_ms = calc_time_series(v_concat, SR)
    
    # filter voltage (to vf)
    vf = butter_filter(v_concat, 
                       order = 3,
                       cutoff = 1e3,
                       sampling_rate = SR)
    
    ### construct current dataframe
    i_hold = I_hold_table[PGF]
    
    # calculate current steps relative to I_hold
    ## rounded to nearest 5
    i_hold_rounded = round_to_base(i_hold, 5)
    
    # rewrite I start for each cells since that varies in the protocols
    cc_sag_parameters['i_start'] = I_hold_table['cc_sag-i_start']
    cc_sag_parameters['i_delta'] = I_hold_table['cc_sag-i_delta']
    
    # get current arrays and list of input current relative to i_hold
    i, i_input = construct_current_array(i_hold = i_hold_rounded,
                                         n_steps = n_steps,
                                         parameters_dict = cc_sag_parameters,
                                         SR_ms = SR_ms)
    
    # split concatenate arrays back to steps wise 
    # needs to occurs after filtering because of the filtering artifact
    v = [None] * n_steps
    t = [None] * n_steps
    
    step_dur = cc_sag_parameters['t_pre'] + cc_sag_parameters['t_stim'] + cc_sag_parameters['t_post']
    step_points = step_dur * SR_ms
    
    pre_points = int(cc_sag_parameters['t_pre'] * SR_ms)
    pre_n_stim_points = int((cc_sag_parameters['t_pre'] + cc_sag_parameters['t_stim']) * SR_ms)
    
    
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
                # calc delta v with mini
                v_pre = v[step_idx][idc_pre]
                v_pre_mean = np.mean(v_pre)
                delta_v = abs(v_step_min - v_pre_mean)
                
                # fit exponential curve
                popt, pcov = sc.optimize.curve_fit(exp_func, x_expFit, v_step_expFit, p0 = [delta_v, *popt_guess[1:]], maxfev = 5000,
                                                   bounds = ([delta_v-3, 0, -200], [delta_v+3, 1, -80]))
                
                r_squared = calc_rsquared_from_exp_fit(x_expFit, v_step_expFit, popt)
                
                
                ### R_INPUT ###
                #R_input = ∆U/∆I
                #∆U = a = popt[0], ∆I = 20 (for first step)
                # calc delta v
                # v_pre = v[step_idx][idc_pre]
                # v_pre_mean = np.mean(v_pre)
                # delta_v = (popt[2] - v_pre_mean)
                
                # reverse sign for r_input calculation
                delta_v = -popt[0]
            
                #delta I
                delta_i = i_input[step_idx]
            
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
                    axs_expfit[step_idx].text(x = -delta_points+100,
                                              y = -145,
                                              s = f'{popt[0]}\n{popt[1]}\n{popt[2]}\nr^2: {r_squared}',
                                              va = 'bottom',
                                              ha = 'left',
                                              fontsize = 8)
                    
                    axs_expfit[step_idx].text(x = 0+100,
                                              y = -145,
                                              s = f'r_input: {r_input}\ntau_mem: {tau_mem}',
                                              va = 'bottom',
                                              ha = 'left',
                                              fontsize = 8)
            
                if r_squared > r_squared_thresh:
                    useful_steps += 1
            
            except (RuntimeError, ValueError):
                print(f'{cell_ID} step number {step_idx+1} has been omitted')
            
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
        vplot_path = os.path.join(vplot_dir, 'cc_IF', 'passive_properties')
        save_figures(fig_expfit, f'{cell_ID}-{PGF}-r_input_n_tau_mem', vplot_path, darkmode_bool)
        
    
    ### break out if fitting to first steps of IF is not successful ###
    n_steps_w_good_fit = len(r_input_calc_df[r_input_calc_df['r_squared'] > 0.75])
    
    if n_steps_w_good_fit < 3:
        print(f'{cell_ID} needs to be discarded')
    
    
    # calc values as mean of 3 steps
    r_input = r_input_calc_df['r_input'].mean()   
    tau_mem = tau_mem_calc_df['tau_mem'].mean()
    
    return r_input, tau_mem
















