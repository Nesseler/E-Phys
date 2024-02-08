# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:51:40 2024

@author: nesseler
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from os.path import join, exists
from os import mkdir

from parameters.directories_win import table_dir, quant_data_dir, vplot_dir
from parameters.PGFs import cc_cntrest_parameters


from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_plotting import get_colors, save_figures
from functions.functions_useful import calc_time_series, butter_filter


table = pd.read_excel(table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')

# loop to create string to include all frequencies in query
PGF = 'cc_cnt_rest'  

# limit lookup table
lookup_table = table.query(f'{PGF}.notnull()')

# cell IDs 
cell_IDs = lookup_table.query('cc_cnt_rest.notnull()').index.to_list()

autocorr_dir = join(quant_data_dir, 'cnt_rest', 'autocorrelations')

# autocorrelation dataframe will be filled for all cells
autocorr_df = pd.DataFrame(columns=cell_IDs)


# plotting specifications
darkmode_bool = False

vplot_bool = True

save_bool = True

autocorr_bool = False

colors_dict = get_colors(darkmode_bool)



for cell_ID in cell_IDs:

    vplot_dir_cell = join(vplot_dir, PGF, cell_ID)
    if not exists(vplot_dir_cell):
        mkdir(vplot_dir_cell)

    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)
    
    # get IF data form file
    _, v, _, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')
    
    # sampling rate in ms
    SR_ms = SR / 1e3
    
    # concatenate individual steps
    n_points = int(np.shape(v)[0] * np.shape(v)[1])
    
    # reconstruct i array from parameters
    i = np.multiply(np.ones(n_points), cc_cntrest_parameters['i_hold'])
    
    v_concat = v.flatten() 

    print(f'Started: {cell_ID}')
    
    t_ms = calc_time_series(v_concat, SR)
    t_s = calc_time_series(v_concat, SR, scale = 's')
    t_total = len(t_s) / SR
    
    
    # filter voltage (to vf)
    vf = butter_filter(v_concat, 
                        order = 3,
                        cutoff = 1e3,
                        sampling_rate = SR)


    # %% autocorrelation
    
    # copy filtered version of voltage trace to circumvent pass by reference
    # create pandas series for its autocorrelation function
    vf_ds = copy(pd.Series(vf))
    
    # define total time frame of autocorrelation
    lag_s = 100
    # convert to datapoints
    lag_points = lag_s * SR
    
    # define time steps to shift series by, in seconds
    lag_steps_s = 1
    # convert to datapoints
    lag_steps_points = lag_steps_s * SR
    
    # create array of lag intervals, datapoints
    lag = np.arange(-lag_points, lag_points + lag_steps_points, lag_steps_points, dtype = int)
    autocorr_values = np.zeros((len(lag)))
    
    for idx, ilag in enumerate(lag):
        autocorr_values[idx] = vf_ds.autocorr(ilag)
        print(f'{idx+1}/{len(lag)}')
    
    # convert lag from datapoints to seconds
    lag_s = [l / SR for l in lag]
    
    # write to autocorrelation dataframe
    autocorr_df['lag_s'] = lag_s
    autocorr_df = autocorr_df.set_index('lag_s')
    autocorr_df[cell_ID] = autocorr_values
    
    
    # %%

    if vplot_bool:
        fig_autocor, ax_autocor = plt.subplots(nrows = 1,
                                                ncols =1,
                                                layout = 'constrained')
        
        fig_autocor.suptitle(f'autocorrelogram {cell_ID}')
        
        ax_autocor.axhline(c='grey', lw=1)
        ax_autocor.axvline(c='grey', lw=1)
        
        ax_autocor.plot(lag_s, autocorr_values, 
                        marker = '',
                        markersize = 5,
                        linestyle = '-',
                        c = colors_dict['primecolor'])
        
        
        ax_autocor.set_xlim([-t_total/3, t_total/3])
        ax_autocor.set_xlabel('Time [s]')
        
        ax_autocor.set_ylim([-0.25, 1.05])
        ax_autocor.set_ylabel('Autocorrelation (Pearson)')
        
        ax_autocor.grid(False)
    
        if save_bool:
            save_figures(fig_autocor, f'autocorrelogram-{cell_ID}', vplot_dir_cell, darkmode_bool)
            
        plt.show()
