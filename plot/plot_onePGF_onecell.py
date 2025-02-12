# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 11:24:17 2024

@author: nesseler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from parameters.directories_win import table_file, figure_dir

# from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size
from functions.functions_constructors import construct_current_array
# from functions.functions_ccIF import get_IF_data
from functions.functions_import import get_traceIndex_n_file, get_cc_data
from functions.functions_useful import calc_time_series, calc_dvdt, calc_dvdt_padded, round_to_base
from functions.functions_filter import butter_filter


cell_ID = 'E-122'
PGF = 'cc_IF'

from parameters.PGFs import cc_th1Ap_parameters, cc_IF_parameters, cc_sag_parameters#, vc_rest_EPSC_parameters

PGFs_dict = {'cc_IF' : cc_IF_parameters,
             'cc_th1AP' : cc_th1Ap_parameters,
             'cc_sag' : cc_sag_parameters}

PGF_parameters = PGFs_dict[PGF]


# init plotting
from functions.initialize_plotting import *

# get hold current as table
I_hold_table = pd.read_excel(table_file, sheet_name="V_or_I_hold", index_col='cell_ID').loc[cell_ID, :]

# %% load data

# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)

# get IF data form file
i, v, t, SR, n_steps = get_cc_data(file_path, traceIndex, 'ms')

# sampling rate in ms
SR_ms = SR / 1e3

# concatenate individual steps
n_points = int(np.shape(i)[0] * np.shape(i)[1])
v_concat = v.flatten() 

# filter voltage (to vf)
vf = butter_filter(v_concat, 
                   order = 3,
                   cutoff = 1e3,
                   sampling_rate = SR)

# split concatenate arrays back to steps wise 
# needs to occurs after filtering because of the filtering artifact
v = [None] * n_steps

step_dur = PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']
step_points = step_dur * SR_ms

# loop through steps to limit voltage trace
for step_idx in np.arange(0, n_steps, 1):
    
    # calc start and stop indices for step
    start_idx = int(step_points * step_idx)
    stop_idx = int(start_idx + step_points)
    
    # set voltage trace of step
    v[step_idx] = vf[start_idx:stop_idx]

# calc time series for one step
t = calc_time_series(v[0], SR_ms)


# ### construct current dataframe
i_hold = I_hold_table[PGF]

# # calculate current steps relative to I_hold
# ## rounded to nearest 5
# i_hold_rounded = round_to_base(i_hold, 5)

# get current arrays and list of input current relative to i_hold
i_cons, i_input = construct_current_array(i_hold = i_hold,
                                          n_steps = n_steps,
                                          parameters_dict = PGF_parameters,
                                          SR_ms = SR_ms)


# %% figure

fig, axs = plt.subplots(nrows = 2, 
                        ncols = 1,
                        # layout = 'tight',
                        sharex = 'col',
                        figsize = get_figure_size(328.67/4, 165.5/2))

fig.subplots_adjust(hspace = 0.06)

for step_idx in np.arange(0, n_steps, 5):
    axs[0].plot(t, i_cons[step_idx], lw = 1, c = 'gray')
    axs[1].plot(t, v[step_idx], lw = 1, c = 'gray')


add_stepidx = 22

axs[0].plot(t, i_cons[add_stepidx], lw = 1, c = colors_dict['primecolor'])
axs[1].plot(t, v[add_stepidx], lw = 1, c = colors_dict['primecolor'])


# #y i
# i_range = [-150, 200]
# axs[0].set_ylabel('Current [pA]')
# axs[0].set_ylim(i_range) 
# axs[0].set_yticks(np.arange(-100, i_range[1] + 1, 100))
# axs[0].set_yticks(np.arange(i_range[0], i_range[1] + 1, 25), minor = True)



# #x
# axs[1].set_xlabel('Time [ms]')
# axs[1].set_xlim([0, step_dur]) 
# axs[1].set_xticks(np.arange(0, step_dur + 1, step_dur / 4))
# axs[1].set_xticks(np.arange(0, step_dur + 1, step_dur / 8), minor = True)    

# #y
# v_range = [-160, 40]
# axs[1].set_ylabel('Voltage [mV]')
# axs[1].set_ylim(v_range) 
# axs[1].set_yticks(np.arange(-150, v_range[1] + 1, 50))
# axs[1].set_yticks(np.arange(v_range[0], v_range[1] + 1, 10), minor = True)


# y
ydict_v = {'ax_min' : -100,
           'ax_max' : 60,
           'pad' : None,
           'step' : 50,
           'stepminor' : 5,
           'label' : '',
           'limits_n_0' : True}

apply_axis_settings(axs[1], axis = 'y', **ydict_v)



# plt.show()

save_figures(fig, f'example_trace-{cell_ID}_{PGF}', figure_dir, darkmode_bool,
             figure_format='both')












