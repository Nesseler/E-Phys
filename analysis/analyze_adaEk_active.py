# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 17:29:10 2025

@author: nesseler
"""

import pandas as pd
import numpy as np
from os.path import join
import scipy as sc
import warnings
from tqdm import tqdm


# custom directories & parameters
from parameters.directories_win import vplot_dir, figure_dir
from parameters.parameters import min_peak_prominence, min_peak_distance

# custom functions
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_useful import butter_filter, calc_dvdt_padded

# define protocol
PGF = 'cc_sup1AP'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF + '_pre', sheet_name = 'PGFs_Syn')

# define conditions
conditions = ['pre', 'post']

# set parameters
cc_sup1APk_parameters = {'t'             : 30,
                         'i_hold'        : 0 ,
                         'SR'            : 50000}

SR_ms = cc_sup1APk_parameters['SR'] / 1e3



# init plotting
from functions.initialize_plotting import * # analysis:ignore


cell_ID = 'E-316'

condition = 'post'

idx_start = int(491*SR_ms)
idx_stop = int(515*SR_ms)
dur = (idx_stop - idx_start) / SR_ms
t = np.arange(0, dur, 1/SR_ms)

# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')

# get data with file path & trace index
_, v, _, _, n_steps = get_cc_data(file_path, traceIndex, scale='s')

# 
v2 = np.empty((n_steps, len(t)))

for step in range(n_steps):
    # filter all data with 1kHz cutoff
    vf = butter_filter(v[step], order=3, cutoff=1e3, sampling_rate=cc_sup1APk_parameters['SR'])
    
    # pads with nans
    vf = np.pad(vf,
                pad_width = (0, 1500000),
                mode = 'constant',
                constant_values = (np.nan,np.nan))

    # write to v2
    v2[step] = vf[idx_start:idx_stop]
    



for step in range(n_steps): 
    plt.plot(v2[step])
    
plt.ylim([-90, 50])
plt.show()

# %%

from functions.functions_extractspike import extract_spike

v3 = v2[1]



# find peaks
idc_spikes, dict_peak = sc.signal.find_peaks(v3, 
                                             prominence = min_peak_prominence, 
                                             distance = min_peak_distance * SR_ms)
    
# calculate spike times in seconds
t_spikes = np.divide(idc_spikes, SR_ms)
n_spikes = len(t_spikes)

    
# calc first derivative
dvdt = calc_dvdt_padded(v3, t)

# define negative dvdt threshold, i.e.: detection of fast AHP
dvdt_n_threshold = -1

# # copy voltage trace to adjustable variable 
# vf_wo_spikes = np.copy(vf)


# plt.plot(t, v3)


# iterate through spikes
# for spike_idx in idc_spikes:

spike_idc, spike_t, spike_v, spike_dvdt = extract_spike(t = t, 
                                                        v = v3, 
                                                        dvdt = dvdt, 
                                                        idx_peak = idc_spikes[0],
                                                        dvdt_n_threshold=dvdt_n_threshold)

plt.plot(spike_v, spike_dvdt)
plt.ylim([-50, 150])
plt.show()





