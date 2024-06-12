# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:01:03 2024

@author: nesseler
"""

from os.path import join
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import seaborn as sbn

from parameters.directories_win import table_file, quant_data_dir, cell_descrip_dir, vplot_dir

from parameters.PGFs import cc_IF_parameters

from functions.functions_useful import calc_time_series, calc_dvdt_padded, butter_filter, round_to_base
from functions.functions_constructors import construct_current_array
from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from getter.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes

vplot_bool = True
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)


# %% load passive properties

# cc_rest
activity_df = pd.read_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col = 'cell_ID')

# cc_IF
passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col = 'cell_ID')
active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')
fstAP_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_col = 'cell_ID')
# IF_step_idc = pd.read_excel(join(cell_descrip_dir, 'ccIF-step_indices.xlsx'), index_col = 'cell_ID')

IF_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_col = 'i_input')
IF_inst_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_col = 'i_input')
IF_inst_initial_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst_initial.xlsx'), index_col = 'i_input')


# calc delta v_thres to v_rest
delta_vrest_to_vthres_df = fstAP_df['v_threshold'] - activity_df['v_rest']

# %% analyzable cells

# protocol 
PGF = 'cc_IF'

# get cell IDs
cell_IDs = get_cell_IDs_one_protocol(PGF)

cell_IDs = ['E-122']

for cell_ID in cell_IDs:
    # test if desired cells in passive properties
    if cell_ID not in passiv_properties_df.index.to_list():
        raise ValueError('Passiv properties not calculated for provided cell_ID!')
    
    # test if desired cells in passive properties
    if cell_ID not in activity_df.index.to_list():
        raise ValueError('Resting membrane properties not calculated for provided cell_ID!')

# get hold current as table
I_hold_table = pd.read_excel(table_file, sheet_name="V_or_I_hold", index_col='cell_ID').loc[cell_IDs, :]

# %% load steps

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

























# %%




plt.plot(IF_df[cell_IDs], label = 'max')
plt.plot(IF_inst_df[cell_IDs], label = 'max_inst')
plt.plot(IF_inst_initial_df[cell_IDs], label = 'max_inst_initial')

plt.legend()











    
