# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:52:47 2024

@author: nesseler
"""


import os
os.chdir('C:/Users/nesseler/E-Phys')
from os.path import join
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sbn
import numpy as np
import scipy as sc

# custom directories & parameters
from parameters.directories_win import table_file, quant_data_dir, cell_descrip_dir, vplot_dir, figure_dir
from parameters.parameters import min_peak_prominence_ccIF, min_peak_distance_ccIF, min_max_peak_width_ccIF, dvdt_threshold, AP_parameters, t_expo_fit, popt_guess, r_squared_thresh
from parameters.PGFs import cc_IF_parameters
# from getter.get_cell_IDs import get_cell_IDs_one_protocol, get_cell_IDs_all_ccAPfreqs

from functions.functions_constructors import construct_current_array
from functions.functions_ccIF import get_IF_data
from functions.functions_import import get_traceIndex_n_file
from functions.functions_useful import calc_time_series, butter_filter, calc_dvdt, calc_dvdt_padded, round_to_base, exp_func, calc_rsquared_from_exp_fit
from functions.functions_plotting import get_colors, get_figure_size, save_figures, plot_t_vs_v, set_font_sizes
from functions.functions_extractspike import get_AP_parameters

IF_step_idc = pd.read_excel(join(cell_descrip_dir, 'ccIF-step_indices.xlsx'), index_col='cell_ID')

# %% settings

vplot_bool = True

darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)


# cell_ID = 'E-122'
PGF = 'cc_IF'


#BAOT
cell_IDs = ['E-130', 'E-089', 'E-152', 'E-109', 'E-112', 'E-186']

# load data
# active membrane properties for rheobase step idx
IF_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_col = 'i_input')
active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')

# import meta data sheet
MetaData = pd.read_excel(table_file,
                      sheet_name="MetaData",
                      index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]


supfig = plt.figure(layout = 'constrained', 
                    figsize=get_figure_size())

subfigs = supfig.subfigures(nrows = 2, 
                            ncols = 3,
                            width_ratios = [1, 1, 1])

subfigs = subfigs.flatten()


for idx_cell, cell_ID in enumerate(cell_IDs):
    
    rheobase_step_idx = IF_step_idc.at[cell_ID, 'rheobase_step_idx']

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
    t_step = calc_time_series(v[0], SR)
    
    # filter voltage (to vf)
    vf = butter_filter(v_concat, 
                       order = 3,
                       cutoff = 1e3,
                       sampling_rate = SR)
    
    # %% i hold and input current
    
    # get hold current as table
    i_hold = pd.read_excel(table_file, sheet_name="V_or_I_hold", index_col='cell_ID').at[cell_ID, PGF]
    
    # calculate current steps relative to I_hold
    ## rounded to nearest 5
    i_hold_rounded = round_to_base(i_hold, 5)
    
    # get current arrays and list of input current relative to i_hold
    i, i_input = construct_current_array(i_hold = i_hold_rounded,
                                         n_steps = n_steps,
                                         parameters_dict = cc_IF_parameters,
                                         SR_ms = SR_ms)
    
    
    # %% start plotting
    
    # zero_input_step_idx = np.where(i_input == 0)[0][0]
    max_freq_step_idx = IF_df[cell_ID].dropna().argmax()
    
    above_max_freq_step_idx = max_freq_step_idx + (max_freq_step_idx - rheobase_step_idx)
    
    if above_max_freq_step_idx >= len(i_input):
        above_max_freq_step_idx = len(i_input) -1
        
    first_step_idx = 0
    last_step_idx = len(i_input) -1
    
    steps_to_plot = [first_step_idx, rheobase_step_idx, max_freq_step_idx, last_step_idx]
    
    
    ax_IF = subfigs[idx_cell].subplots(nrows = 2,
                                       ncols = 1,
                                       gridspec_kw={'height_ratios': [1, 4]})
    
    subfigs[idx_cell].suptitle(cell_ID)
    
    set_font_sizes()
    
    # x axis at 0
    ax_IF[0].axhline(y = 0, xmin = 0, xmax = 1500,
                     linestyle = '--',
                     lw = 1,
                     color = 'gray')
    
    # x axis at -85
    ax_IF[1].axhline(y = -85, xmin = 0, xmax = 1500,
                     linestyle = '--',
                     lw = 1,
                     color = 'gray')
    
    
    
    # initialise color code
    ## same normalisation for all cells
    norm_min = IF_df.index.to_list()[0]
    norm_max = IF_df.index.to_list()[-1]
    
    ## cell specific normalisation
    # norm_min = i_input[0]
    # norm_max = i_input[above_max_freq_step_idx]
    cmap_str = 'rainbow'
    
    
    norm = mpl.colors.Normalize(norm_min, norm_max)
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
    cmap.set_array(list(i_input))
    
    

    
    
    step_dur = cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim'] + cc_IF_parameters['t_post']
    step_points = step_dur * SR_ms
    
    pre_points = int(cc_IF_parameters['t_pre'] * SR_ms)
    pre_n_stim_points = int((cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim']) * SR_ms)
        
    n_plotted_steps = 0
        
    # loop through steps to limit voltage trace
    for step_idx in np.arange(0, n_steps, 1):
         
        if step_idx in steps_to_plot:
            # calc start and stop indices for step
            start_idx = int(step_points * step_idx)
            stop_idx = int(start_idx + step_points)
        
            # set voltage trace of step
            v_step = vf[start_idx:stop_idx]
            i_step = i[step_idx]
    
    
            # plot each step
            ax_IF[0].plot(t_step, i_step, 
                          c=cmap.to_rgba(i_input[step_idx]),
                          lw = 1)
            
            ax_IF[1].plot(t_step, v_step, 
                          c=cmap.to_rgba(i_input[step_idx]),
                          lw = 1)
            
            # find peaks
            idc_peaks, dict_peak = sc.signal.find_peaks(v_step, 
                                                        prominence = min_peak_prominence_ccIF, 
                                                        distance = min_peak_distance_ccIF * (SR_ms),
                                                        width = np.multiply(min_max_peak_width_ccIF, SR_ms))
            # get times of spikes
            t_spikes = np.divide(idc_peaks, SR_ms)
            
            # eventplot
            ax_IF[1].eventplot(t_spikes, 
                               lineoffsets=55+(n_plotted_steps*5), 
                               colors = cmap.to_rgba(i_input[step_idx]), 
                               linelengths=4.5,
                               linewidth = 1)
            
            n_plotted_steps += 1
    
    
    [ax.grid(False) for ax in ax_IF]
    
    
    # if idx_cell in [2, 5]:
    #     # colorbar
    #     subfigs[idx_cell].colorbar(cmap, ax = ax_IF[1])
    
    
    ## current
    # y
    ax_IF[0].set_ylabel('')
    ax_IF[0].set_ylim([-100-5, 400+5])
    ax_IF[0].set_yticks(ticks = np.arange(-100, 400 + 1, 100),
                        labels = [])
    ax_IF[0].set_yticks(ticks = np.arange(-100, 400 + 1, 50),
                        labels = [], minor = True)
    
    ## voltage
    # y
    ax_IF[1].set_ylabel('')
    ax_IF[1].set_ylim([-150-2, 75+2])
    ax_IF[1].set_yticks(ticks = np.arange(-150, 75 + 1, 50),
                        labels = [])
    ax_IF[1].set_yticks(ticks = np.arange(-150, 75 + 1, 10),
                        labels = [], minor = True)
    
    
    ax_IF[0].set_xlabel('')
    ax_IF[0].set_xlim([-10, 1500+10])
    ax_IF[0].set_xticks(ticks = [],
                        labels = [])
    
    
    
    if idx_cell in [3, 4, 5]:
        # x
        ax_IF[1].set_xlabel('Time [ms]')
        ax_IF[1].set_xlim([-10, 1500+10])
        ax_IF[1].set_xticks(ticks = np.arange(0, 1500 + 1, 500))
        ax_IF[1].set_xticks(ticks = np.arange(0, 1500 + 1, 50), minor = True)
        
    else:
        # x
        ax_IF[1].set_xlabel('')
        ax_IF[1].set_xlim([-10, 1500+10])
        ax_IF[1].set_xticks(ticks = np.arange(0, 1500 + 1, 500),
                            labels = [])
        ax_IF[1].set_xticks(ticks = np.arange(0, 1500 + 1, 50),
                            labels = [], minor = True)
    
    # ## current
    # # y
    # ax_IF[0].set_ylabel('Current\n[pA]')
    # ax_IF[0].set_ylim([-100-5, 400+5])
    # ax_IF[0].set_yticks(ticks = np.arange(-100, 400 + 1, 100),
    #                     labels = [None, 0, None, 200, None, 400])
    # ax_IF[0].set_yticks(ticks = np.arange(-100, 400 + 1, 50), minor = True)
    
    # ## voltage
    # # y
    # ax_IF[1].set_ylabel('Voltage [mV]')
    # ax_IF[1].set_ylim([-150-2, 75+2])
    # ax_IF[1].set_yticks(ticks = np.arange(-150, 75 + 1, 50))
    # ax_IF[1].set_yticks(ticks = np.arange(-150, 75 + 1, 10), minor = True)
    

    
    
    ax_IF[0].spines['bottom'].set_visible(False)
    
    # despine
    [ax.spines[spine].set_visible(False) for ax in ax_IF for spine in ['top', 'right']]
    [ax.spines['bottom'].set_bounds([0, 1500]) for ax in ax_IF]
    
    
# supfig.colorbar(cmap)
   
    
plt.show()






# save figure
save_figures(supfig, 'overview-ccIF-example_traces-color_coded', figure_dir, darkmode_bool)

























