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


from plotting_functions import get_figure_size, get_colors


# %%

## n/3 by 3 subplots with each a 30 sec resting activity plot


from patchview.HekaIO.HekaHelpers import HekaBundleInfo

lookup_table = pd.read_csv('/Users/moritznesseler/local E-Phys/cc_IF.csv', 
                           delimiter=';',
                           index_col='cell_ID')


data_folder = '/Users/moritznesseler/local E-Phys'

figure_dir = '/Users/moritznesseler/local E-Phys'

darkmode_bool = True

# convert indices of dataframe to list to loop through
all_cell_IDs = lookup_table.index.to_list()

for cell_idx, cell_ID in enumerate(all_cell_IDs):
    print(cell_ID)
    
    
    
    group_idx = int(lookup_table.at[cell_ID, 'group_idx'])
    series_idx = int(lookup_table.at[cell_ID, 'series_idx'])

    # construct traceIndex with indices
    traceIndex = [group_idx, series_idx, 0, 0]

    current_file = lookup_table.at[cell_ID, 'file']
    
    data_file_path = os.path.join(data_folder, current_file + '.dat')

    data_file_path_str = fr"{data_file_path}"



    bundleTester = HekaBundleInfo(data_file_path_str)

    ## Get data from a single sweep and single channel
    data = bundleTester.getSeriesData(traceIndex)

    SR = bundleTester.getSeriesSamplingRate(traceIndex)

    v = data[:,0,:] * 1e3
    i = data[:,1,:] * 1e12
    t = calc_time_series(v, SR, scale='ms')
    
    fig_IF, axs_IF = plt.subplots(2, 1, 
                                  layout = 'constrained',
                                  figsize=get_figure_size(),
                                  gridspec_kw={'height_ratios': [2,8]},
                                  sharex='col')


    n_steps = np.shape(v)[1]
    step_idx = np.arange(0,n_steps)

    # axs_IF[1][0].plot(t, v, c = prime_color)


    vs = np.transpose(v)
    i_s = np.transpose(i)

    segs = [np.column_stack([t, vi]) for vi in vs]


    # line_segments = mtl.collections.LineCollection(segs, array=step_idx, cmap="plasma", linewidth=0.5)

    # axs_IF[1][0].add_collection(line_segments)
    #cbar = fig_IF.colorbar(line_segments)




    norm = mtl.colors.Normalize(vmin=0, vmax=n_steps)
    cmap = mtl.cm.ScalarMappable(norm=norm, cmap="plasma")
    cmap.set_array(step_idx)


    fig_IF.colorbar(cmap, label = 'Number of steps', ax=axs_IF[1])


    n_peaks = []

    #example_steps = [0, 10, 25, 60]

    example_steps = step_idx

    for i, vi in enumerate(vs):
        
        if i in example_steps:
            axs_IF[1].plot(t, vi, c=cmap.to_rgba(i + 1))
            axs_IF[0].plot(t, i_s[i], c=cmap.to_rgba(i + 1))
        
        v_pulse = vi[(250*50):(1250*50)]
        i_pulse = i_s[i][(250*50):(1250*50)]
        
        idx_peaks, dict_peak = sc.signal.find_peaks(v_pulse, prominence = 20, width = 50)
        
        n_peaks.append(len(idx_peaks))
        
        #axs_IF[1].plot(n_peaks, step_idx[:i+1])
        

    axs_IF[0].set_ylim([-100,300])

    axs_IF[0].set_yticks(np.arange(-50,251,50))
    axs_IF[0].set_ylabel('Inj. current\n[pA]')



    axs_IF[1].set_ylim([-150,50])
    axs_IF[1].set_ylabel('Membrane potential\n[mV]')

    axs_IF[1].set_xlabel('Time [ms]')

    axs_IF[1].set_xlim([0,1500])
    axs_IF[1].set_xticks(np.arange(0,1501,250))
    
    plt.pause(0.4)
    
    plt.show()
    





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
    
    
    
    
    
    
    
    
    