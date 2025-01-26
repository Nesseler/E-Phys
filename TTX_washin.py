# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:07:52 2024

@author: nesseler
"""

import scipy as sc
import numpy as np
import pandas as pd
from os.path import join
from tqdm import tqdm


# define protocol to analyze
PGF = 'vc_TTX_washin'

# get cell IDs
from functions.get_cell_IDs import get_cell_IDs_one_protocol
cell_IDs = get_cell_IDs_one_protocol(PGF, sheet_name = 'PGFs_Syn')
cell_IDs_leak = get_cell_IDs_one_protocol(PGF + '_leak', sheet_name = 'PGFs_Syn')

# load parameters
from parameters.PGFs import vc_TTX_washin_parameters, vc_TTX_washin_leak_parameters

from functions.initialize_plotting import *

# verification plots
vplots = False


# %% set peak detection parameters

peak_prominence = 40 #pA
peak_width = 50 # points 


# %% load and process

from functions.functions_import import get_traceIndex_n_file, get_vc_data, get_vc_leak_data
from functions.functions_useful import calc_time_series
from parameters.directories_win import table_file, vplot_dir
import re

# init dataframe for measurements
peakCurrents = pd.DataFrame(columns = cell_IDs + cell_IDs_leak)

# cell_ID = 'E-268' #E-268

# load lookup table
lookup_table = pd.read_excel(table_file,
                             sheet_name = 'PGFs_Syn',
                             index_col = 'cell_ID')


for cell_ID in tqdm(cell_IDs + cell_IDs_leak):
           
    # check if cell has leak subtraction
    if cell_ID in cell_IDs_leak:
        
        PGF = 'vc_TTX_washin_leak'
        PGF_params = vc_TTX_washin_leak_parameters
        
        # get the traceIndex and the file path string for data import functions
        traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name='PGFs_Syn')
        
        # load data
        i, v, ileak, t, SR, n_steps = get_vc_leak_data(file_path, traceIndex, scale = 'ms')
        
    else: 
        
        PGF = 'vc_TTX_washin'
        PGF_params = vc_TTX_washin_parameters
        
        # get the traceIndex and the file path string for data import functions
        traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name='PGFs_Syn')
        
        # load data
        i, v, t, SR, n_steps = get_vc_data(file_path, traceIndex, scale = 'ms')
        


    # convert sampling rate
    SR_ms = SR/1000
    
    # load which steps to use
    ## get steps from lookup table
    steps_str = lookup_table.at[cell_ID, PGF+'_steps']

    # find integers in string
    steps_start_stop = [int(i) - 1 for i in steps_str.split('-') if i.isdigit()]
    
    # convert to indices
    steps_start_stop[0] = steps_start_stop[0]-1
    
    # create iteratable list
    steps_to_use = np.arange(steps_start_stop[0], steps_start_stop[1])
    
    # redefine number of steps
    n_steps = len(steps_to_use)
    
    # iterate through steps
    for step_idx in range(n_steps):
        
        # limit data to stimulation period
        step = i[step_idx][int(PGF_params['t_pre'] * SR_ms):int(PGF_params['t_post'] * SR_ms + PGF_params['t_stim'] * SR_ms)]
    
        # find peaks
        peaks, properties = sc.signal.find_peaks(-step, prominence = peak_prominence, width = peak_width)
    
        # get peak measurements
        peak_idx = peaks[0]
        peak_current = step[peak_idx]
        t_peak = peak_idx / SR_ms
        
        # get measurements of peak base
        leftbase_idx = properties['left_bases'][0]
        leftbase = step[leftbase_idx]
        t_leftbase = leftbase_idx / SR_ms
        
        # calculate peak current data
        delta = leftbase - peak_current
        
        # write to dataframe
        peakCurrents.at[step_idx, cell_ID] = delta
    
        # %
        # create verification plot
        if vplots:
        
            vfig, vax = plt.subplots(nrows = 1,
                                     ncols = 1,
                                     figsize = get_figure_size(width = 150, height = 120),
                                     dpi = 100,
                                     layout = 'constrained')
            
            vfig.suptitle(f'{cell_ID} step: {step_idx}')
            
            # get trace
            current_trace = i[step_idx][int((PGF_params['t_pre']-5) * SR_ms):int((PGF_params['t_post']+5) * SR_ms + PGF_params['t_stim'] * SR_ms)]
        
            # calc time series
            t_plot = np.arange(-5, 25, step = 1/SR_ms)
        
            vax.plot(t_plot,
                     current_trace,
                     c = colors_dict['primecolor'],
                     lw = 1)
            
            # plot peak
            vax.scatter(t_peak, peak_current,
                        marker = 'x',
                        c = 'r',
                        s = 20,
                        zorder = 2)
            
            # plot leftbase
            vax.scatter(t_leftbase, leftbase,
                        marker = 'x',
                        c = 'g',
                        s = 20,
                        zorder = 2)
            
            # y axis
            apply_axis_settings(ax = vax, 
                                axis = 'y', 
                                ax_min = -6000, 
                                ax_max = 2000, 
                                pad = 50, 
                                step = 2000, 
                                stepminor = 250, 
                                label = 'Current [pA]')

            # x axis
            apply_axis_settings(ax = vax, 
                                axis = 'x', 
                                ax_min = -5, 
                                ax_max = 25, 
                                pad = 1, 
                                step = 5, 
                                stepminor = 1, 
                                label = 'Time [ms]' )
            
            # despine
            [vax.spines[spine].set_visible(False) for spine in ['top', 'right']]
            
            # display plot
            plt.show()
            
            # set directory for figure
            vfig_dir = join(vplot_dir, 'TTX_washin')

            # save figure
            save_figures(vfig, f'figure-TTXwashin-{cell_ID}-step{step_idx}', 
                          save_dir = vfig_dir,
                          darkmode_bool = darkmode_bool,
                          figure_format = 'both')

            
            

    
# %%

from parameters.directories_win import quant_data_dir, figure_dir

peakCurrents.to_excel(join(quant_data_dir, 'TTX_washin', 'TTX_washin-peakCurrents.xlsx'), index_label = 'step_idc')

# %%

# min max normalize currents
mm_peakCurrents = (peakCurrents - peakCurrents.min()) / (peakCurrents.max() - peakCurrents.min())

plt_cell_IDs = cell_IDs + cell_IDs_leak
# plt_cell_IDs = cell_IDs_leak

#  plot graph
fig, ax = plt.subplots(nrows = 1,
                       ncols = 1,
                       figsize = get_figure_size(width = 328.67/2),
                       dpi = 300,
                       layout = 'constrained')

# create x dimension
x = np.arange(0, 120) * 5 + 2.5

for cell_ID in plt_cell_IDs:
    
    
    if cell_ID in cell_IDs_leak:
        linestyle = 'dashed'
    else:
        linestyle = 'solid'

    # plot
    ax.plot(x, mm_peakCurrents[cell_ID],
            c = 'grey',
            lw = 1,
            ls = linestyle)
    
# calc mean
mean_mm_peakcurrent = mm_peakCurrents[plt_cell_IDs].mean(axis = 1)
std_mm_peakcurrent = mm_peakCurrents[plt_cell_IDs].std(axis = 1)

# plot mean and std
# ax.fill_between(x = x,
#                 y1 = list(mean_mm_peakcurrent+std_mm_peakcurrent),
#                 y2 = list(mean_mm_peakcurrent-std_mm_peakcurrent),
#                 facecolor = 'grey',
#                 alpha = 0.5)

# plot mean and std
ax.plot(x, mean_mm_peakcurrent,
        c = colors_dict['primecolor'],
        lw = 2)

# y axis
apply_axis_settings(ax = ax, 
                    axis = 'y', 
                    ax_min = 0, 
                    ax_max = 1, 
                    pad = 0.05, 
                    step = 0.2, 
                    stepminor = 0.1, 
                    label = 'Current [pA]' )

# x axis
apply_axis_settings(ax = ax, 
                    axis = 'x', 
                    ax_min = 0, 
                    ax_max = 600, 
                    pad = 30, 
                    step = 60, 
                    stepminor = 30, 
                    label = 'Time [s]' )
# despine
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# show figure
plt.show()

# set directory for figure
fig_dir = join(figure_dir, 'temp_figs')

# save figure
save_figures(fig, 'figure-TTXwashin-leakonly-minmax', 
              save_dir = fig_dir,
              darkmode_bool = darkmode_bool,
              figure_format = 'both')

