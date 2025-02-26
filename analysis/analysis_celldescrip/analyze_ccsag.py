# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 09:42:29 2024

@author: nesseler
"""

import warnings
warnings.warn('Script not up-to-date!')
# %%

from os.path import join
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

from parameters.directories_win import table_file, quant_data_dir, cell_descrip_dir, vplot_dir

from parameters.PGFs import cc_IF_parameters, cc_sag_parameters
from parameters.parameters import min_peak_prominence_ccIF, min_peak_distance_ccIF, min_max_peak_width_ccIF

from functions.functions_useful import round_to_base, calc_time_series, butter_filter, calc_dvdt_padded
from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_constructors import construct_current_array
from functions.functions_extractspike import get_AP_parameters
from functions.functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes

from getter.get_cell_IDs import get_cell_IDs_one_protocol, get_cell_IDs_all_ccAPfreqs

vplot_bool = True
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)


# %% analyzable cells

# test cell
cell_ID = 'E-140'
# get cell IDs
cell_IDs = get_cell_IDs_one_protocol('cc_IF')

sag_potential = -130 #mV

# load passive properties
passive_prop = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col= 'cell_ID')

# excel sheet with PGF indices as lookup table
lookup_table = pd.read_excel(table_file, sheet_name="PGFs", index_col='cell_ID')

# get hold current as table
I_hold_table = pd.read_excel(table_file, sheet_name="V_or_I_hold", index_col='cell_ID').loc[cell_IDs, :]

cellIDs_toAnalyze = pd.DataFrame()

# output dataframe
sag_df = pd.DataFrame()

for cell_ID in cell_IDs:
    
    cell_toAnalyze = True
    
    # calc current to polarize cell to given potential
    # ΔI = ΔU / R_input
    necessary_deltaI = ( (cc_sag_parameters['v_hold_pre'] - sag_potential) / passive_prop.at[cell_ID, 'r_input'] ) *1e3 #pA
    
    # round to nearest 10 base, since cc_sag uses delta I of 10 per step
    necessary_deltaI_rounded = round_to_base(necessary_deltaI, 10)
    
    # check if necessary deltaI is reached in ccIF
    ## get holding current and substract current of first step
    max_deltaI_ccIF = I_hold_table.at[cell_ID, 'cc_IF'] - cc_IF_parameters['i_start']
       
    if max_deltaI_ccIF >= necessary_deltaI:
        PGF_toUse = 'cc_IF'
    else:
        PGF_toUse = 'cc_sag'
    
    # try loading the desired protocol
    try:
        # get indices of current cell with the dataframe containing all indices    
        group_idx = int(lookup_table.at[cell_ID, 'group'])-1
        series_idx = int(lookup_table.at[cell_ID, f'{PGF_toUse}'])-1
    
    # accept specific ValueError where a nan value cannot be converted to integer
    # otherwise raise same ValueError and message
    except ValueError as e:
        if str(e) == 'cannot convert float NaN to integer':
            print(f'{cell_ID} will be skipped since {PGF_toUse} is not available')
            cell_toAnalyze = False
        else:
            raise ValueError(e)
        
    # repeat for cc_sag protocol  
    max_deltaI_ccsag = I_hold_table.at[cell_ID, 'cc_sag'] - I_hold_table.at[cell_ID, 'cc_sag-i_start']
    
    if max_deltaI_ccsag <= necessary_deltaI:
        print(f'{cell_ID} cc_sag does not reach (needed: {str(necessary_deltaI)} pA, max: {str(max_deltaI_ccsag)} pA)')
        cell_toAnalyze = False
        
    if cell_toAnalyze:
        cellIDs_toAnalyze.at[cell_ID, 'toAnalyze'] = 1
        cellIDs_toAnalyze.at[cell_ID, 'PGF_toUse'] = PGF_toUse
        cellIDs_toAnalyze.at[cell_ID, 'deltaI_tosag'] = necessary_deltaI
        cellIDs_toAnalyze.at[cell_ID, 'deltaI_tosag_rounded'] = necessary_deltaI_rounded

    # insert deltaI to sag
    sag_df.at[cell_ID, 'deltaI_tosag'] = necessary_deltaI
    sag_df.at[cell_ID, 'deltaI_tosag_rounded'] = necessary_deltaI_rounded


# %% load data and analyze

for cell_ID in cellIDs_toAnalyze.index.to_list():

# cell_ID = 'E-187'

    print(f'Started: {cell_ID}')
    
    PGF_toUse = cellIDs_toAnalyze.at[cell_ID, 'PGF_toUse']
    
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF_toUse, cell_ID)
    
    # get IF data form file
    i, v, t, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')
    
    # sampling rate in ms
    SR_ms = SR / 1e3
    
    # concatenate individual steps
    n_points = int(np.shape(i)[0] * np.shape(i)[1])
    
    
    
    ### construct current dataframe
    if PGF_toUse == 'cc_sag':
        PGF_parameters = {'i_hold' : I_hold_table.at[cell_ID, PGF_toUse],
                          'i_start' : I_hold_table.at[cell_ID, 'cc_sag-i_start'],
                          'i_delta' : I_hold_table.at[cell_ID, 'cc_sag-i_delta'],
                          't_pre' : cc_sag_parameters['t_pre'],
                          't_stim' : cc_sag_parameters['t_stim'],
                          't_post' : cc_sag_parameters['t_post']
                          }
        
    elif PGF_toUse == 'cc_IF':
        cc_IF_parameters['i_hold'] = I_hold_table.at[cell_ID, PGF_toUse]
        PGF_parameters = cc_IF_parameters
         
     
    # filter voltage (to vf)
    v_concat = v.flatten()
    t_ms = calc_time_series(v_concat, SR)
    t_s = calc_time_series(v_concat, SR, scale = 's')
    
    vf = butter_filter(v_concat, 
                       order = 3,
                       cutoff = 1e3,
                       sampling_rate = SR)   
        
    # split voltage back to steps (after filtering)
    v = [None] * n_steps
    step_dur = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']
    step_points = step_dur * SR_ms
    
    # loop through steps to limit voltage trace
    for step_idx in np.arange(0, n_steps, 1):
        # calc start and stop indices for step
        start_idx = int(step_points * step_idx)
        stop_idx = int(start_idx + step_points)
        
        # set voltage trace of step
        v_step = vf[start_idx:stop_idx]
        
        # assign voltage and time traces to arrays
        v[step_idx] = v_step
        
    # time series
    t_step = calc_time_series(v_step, sampling_rate=SR, scale='ms')
    
    # get current arrays and list of input current relative to i_hold
    i, i_input = construct_current_array(i_hold = PGF_parameters['i_hold'],
                                          n_steps = n_steps,
                                          parameters_dict = PGF_parameters,
                                          SR_ms = SR_ms)
    
    # find closest i input to i_necessary
    deltaI_tosag = cellIDs_toAnalyze.at[cell_ID, 'deltaI_tosag']
    closest_steps_idx = np.argmin(np.abs(i_input + deltaI_tosag))
    
    # set v_step for sag potential calculation
    v_step_sag = v[closest_steps_idx]
    
    
    
    # calc indices for last 10 % of points for stim period 
    stim_idc = np.arange(start = PGF_parameters['t_pre'] * SR_ms, stop = (PGF_parameters['t_pre'] + PGF_parameters['t_stim']) * SR_ms,
                         dtype = int)
    
    # limit voltage trace to just the stimulation period
    v_stim = v_step_sag[stim_idc]
    t_stim = t_step[stim_idc]
    
    # find min in first 25 % 
    n_points_thirtypercent = int(PGF_parameters['t_stim'] * SR_ms * 0.3)
    n_points_pre = int(PGF_parameters['t_pre'] * SR_ms)
    
    v_first_thirtypercent = v_stim[:n_points_thirtypercent]
    
    # find local min and its index
    v_step_min = np.min(v_first_thirtypercent)
    v_step_min_idx = np.argmin(v_first_thirtypercent) + n_points_pre
    t_step_min = np.argmin(v_first_thirtypercent) / SR_ms
    
    # calc mean in last 10 % of stimulation period
    n_points_tenpercent = int(PGF_parameters['t_stim'] * SR_ms * 0.1)
    v_last_tenpercent = v_stim[-n_points_tenpercent:]
    v_mean_last_tenpercent = np.mean(v_last_tenpercent)
    
    # calc the sag potential delta
    sag_delta = np.abs(v_step_min - v_mean_last_tenpercent)
    
    
    # write to output dataframe
    sag_df.at[cell_ID, 'v_min'] = v_step_min
    sag_df.at[cell_ID, 't_min'] = t_step_min
    sag_df.at[cell_ID, 'v_mean_steadystate'] = v_mean_last_tenpercent
    sag_df.at[cell_ID, 'sag_delta'] = sag_delta
    
    
    
    ### rebound spikes ###
    
    # calc array for post stim indices
    idc_post = np.arange(start = (PGF_parameters['t_pre'] + PGF_parameters['t_stim']) * SR_ms, stop = (PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']) * SR_ms,
                         dtype = int)
    
    # limit voltage step
    v_post = v_step_sag[idc_post]
    t_post = t_step[idc_post]
    dvdt_post = calc_dvdt_padded(v_post, t_post)
    
    # find peaks
    idc_peaks, dict_peak = sc.signal.find_peaks(v_post, 
                                                prominence = min_peak_prominence_ccIF, 
                                                distance = min_peak_distance_ccIF * (SR_ms),
                                                width = np.multiply(min_max_peak_width_ccIF, SR_ms))
    
    t_peaks = (idc_peaks / SR_ms) + (PGF_parameters['t_pre'] + PGF_parameters['t_stim'])
    
    number_of_reboundspikes = len(idc_peaks)
    
    if number_of_reboundspikes > 0:
    
        # get AP parameter from rebound spikes
        rebound_spike_parameters, rebound_spike_v = get_AP_parameters(t_post, v_post, dvdt_post, idc_peaks)
        
        # save parameter from first rebound spike
        sag_df.at[cell_ID, 'n_reboundspikes'] = number_of_reboundspikes
        sag_df.at[cell_ID, 'reboundspike_t_peak'] = rebound_spike_parameters.at[0, 't_peaks'] - (PGF_parameters['t_pre'] + PGF_parameters['t_stim'])
        sag_df.at[cell_ID, 'reboundspike_t_threshold'] = rebound_spike_parameters.at[0, 't_threshold'] - (PGF_parameters['t_pre'] + PGF_parameters['t_stim'])
        sag_df.at[cell_ID, 'reboundspike_v_threshold'] = rebound_spike_parameters.at[0, 'v_threshold']
        sag_df.at[cell_ID, 'reboundspike_v_amplitude'] = rebound_spike_parameters.at[0, 'v_amplitude']
        sag_df.at[cell_ID, 'reboundspike_t_toPeak'] = rebound_spike_parameters.at[0, 't_toPeak']
        sag_df.at[cell_ID, 'reboundspike_t_rise'] = rebound_spike_parameters.at[0, 't_rise']
        sag_df.at[cell_ID, 'reboundspike_FWHM'] = rebound_spike_parameters.at[0, 'FWHM']
    
    else:
        sag_df.at[cell_ID, 'n_reboundspikes'] = number_of_reboundspikes
    
    
    if vplot_bool:
        from matplotlib.patches import Rectangle
    
        fig_sag, ax_sag = plt.subplots(figsize = get_figure_size(width = 247.75), 
                                       layout = 'constrained')
        
        set_font_sizes()
        
        ax_sag.set_title(cell_ID)
        
        # entire step
        ax_sag.plot(t_step, v_step_sag, c = colors_dict['primecolor'], lw = 1)
        
        # period to look for min
        ax_sag.plot(t_stim[:n_points_thirtypercent], v_first_thirtypercent, c = colors_dict['color2'], lw = 1)
        
        # steady state for mean
        ax_sag.plot(t_stim[-n_points_tenpercent:], v_last_tenpercent, c = colors_dict['color2'], lw = 1)
        
        # horizontal line
        ax_sag.hlines(y = v_step_min, 
                      xmin = PGF_parameters['t_pre'] + t_step_min, 
                      xmax = (PGF_parameters['t_pre'] + PGF_parameters['t_stim']), 
                      colors=colors_dict['color3'],
                      lw = 1)
        
        # steady state horizontal line
        ax_sag.hlines(y = v_mean_last_tenpercent, 
                      xmin = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] * 0.9, 
                      xmax = (PGF_parameters['t_pre'] + PGF_parameters['t_stim']), 
                      colors=colors_dict['color3'],
                      lw = 1)
        
        # vertical line
        ax_sag.vlines(x = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] * 0.95,
                      ymin = v_step_min,
                      ymax = v_mean_last_tenpercent,
                      colors=colors_dict['color3'],
                      lw = 1)
        
        # spike events
        ax_sag.eventplot(t_peaks, lineoffsets=50, colors = 'r', linelengths=5)
        
        # inset marker
        box_pad = np.abs(sag_delta)
        box_ymin = v_step_min - box_pad
        box_xmin = PGF_parameters['t_pre'] - 20
        box_height = np.abs(sag_delta) + 2*box_pad
        box_width = PGF_parameters['t_stim'] + 40
        
        ax_sag.add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                                   width = box_width, 
                                   height = box_height,
                                   fill = False,
                                   color = colors_dict['primecolor'],
                                   linestyle = '--'))
        
        ax_sag.set_ylabel('Voltage [mV]')
        ax_sag.set_yticks(np.arange(-140, 60+1, 20))
        ax_sag.set_yticks(np.arange(-140, 60+1, 5), minor = True)
        ax_sag.set_ylim([-140, 60])
        
        ax_sag.set_xlabel('Time [ms]')
        ax_sag.set_xlim([0, PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']])
        ax_sag.set_xticks(np.arange(0, PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']+1, 250))
        ax_sag.set_xticks(np.arange(0, PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']+1, 50), minor = True)
        
        # add inset
        ## ([left, bottom, width, height]), percentages
        ax_inset = fig_sag.add_axes([0.175, 0.4, 0.55, 0.45])
    
        # plot stimulation period
        ax_inset.plot(t_step, v_step_sag, c = colors_dict['primecolor'], lw = 1)
        
        # period to look for min
        ax_inset.plot(t_stim[:n_points_thirtypercent], v_first_thirtypercent, c = colors_dict['color2'], lw = 1)
        
        # steady state for mean
        ax_inset.plot(t_stim[-n_points_tenpercent:], v_last_tenpercent, c = colors_dict['color2'], lw = 1)
        
        # set y ticks and y limits
        ax_inset.set_yticks(ticks = np.arange(-140, -100+1, 5), labels = [])
        ax_inset.set_ylim([box_ymin, box_ymin + box_height])
    
        # set x ticks
        ax_inset.set_xticks(ticks = np.arange(PGF_parameters['t_pre'], (PGF_parameters['t_pre'] + PGF_parameters['t_stim'])+1, 250), labels = [])
        ax_inset.set_xticks(ticks = np.arange(PGF_parameters['t_pre'], (PGF_parameters['t_pre'] + PGF_parameters['t_stim'])+1, 50), labels = [], minor = True)
        ax_inset.set_xlim([box_xmin, box_xmin + box_width])
    
        # sag delta vertical line
        ax_inset.vlines(x = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] * 0.95,
                        ymin = v_step_min,
                        ymax = v_mean_last_tenpercent,
                        colors=colors_dict['color3'], 
                        lw = 1)
        
        # horizontal line
        ax_inset.hlines(y = v_step_min, 
                        xmin = PGF_parameters['t_pre'] + t_step_min, 
                        xmax = (PGF_parameters['t_pre'] + PGF_parameters['t_stim']), 
                        colors=colors_dict['color3'],
                        lw = 1)
        
        # steady state horizontal line
        ax_inset.hlines(y = v_mean_last_tenpercent, 
                        xmin = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] * 0.9, 
                        xmax = (PGF_parameters['t_pre'] + PGF_parameters['t_stim']), 
                        colors=colors_dict['color3'],
                        lw = 1)
        
        # text
        ax_inset.text(x = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] * 0.90, 
                      y = v_step_min + sag_delta / 2, 
                      s = f'Δsag = {round(sag_delta, 3)} mV', 
                      ha = 'right', va = 'center',
                      size = 9)
        
        plt.show()
            
        save_figures(fig_sag, f'{cell_ID}-sag_step', join(vplot_dir, 'cc_sag'), darkmode_bool)


# save activity dataframe to quant data folder
sag_df.to_excel(join(cell_descrip_dir, 'cc_sag-sagdelta.xlsx'), index_label='cell_ID')

        
# %%


import seaborn as sbn 

sbn.jointplot(data = sag_df,
              x = 'n_reboundspikes',
              y = 'sag_delta')







