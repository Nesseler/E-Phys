# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:27:54 2024

@author: nesseler
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt


# custom directories & parameters
from PGFs import cc_APs_parameters, cc_APs_n_stims

from directories_win import cell_descrip_dir, quant_data_dir, figure_dir


frequencies = list(cc_APs_parameters.keys())

# initialise dataframe to populate in loop
tAPs_df = pd.DataFrame()

# set directory to APs subfolder to get parameters of action potentials
quant_dir = os.path.join(quant_data_dir, 'APs')

# create a list of cell-IDs that have subfolder in the APs folder i.e. that have been analysed
cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]



def import_AP_measurement_all_freqs(frequencies, parameter, cell_ID):

    # create file path
    cell_data_path = os.path.join(quant_dir, cell_ID)
    
    # get number of APs for each stim
    measurement_df = pd.DataFrame()
    
    for frequency in frequencies:
    # frequency = '1Hz'
            
        # create file path
        file_path = os.path.join(cell_data_path, f'{cell_ID}_{frequency}.xlsx')
        
        # read excel file
        AP_all_params = pd.read_excel(file_path, index_col = 'idx_step')
        
        # write to dictionary
        measurement_df[frequency] = AP_all_params[parameter]

    return measurement_df



# %% get data

# start with one cell
cell_ID = 'E-092'
t_APs = import_AP_measurement_all_freqs(frequencies, 't_peaks', cell_ID)



# %% create dataframe with times that

# create list with interstimulus interval per stimulation frequency
ISIs = [sum(cc_APs_parameters[f].values()) for f in frequencies]

# create dataframe with all stimuation time of every frequency
t_stims_df = pd.DataFrame()

for idx_freq, freq in enumerate(frequencies):
    
    t_pre = cc_APs_parameters[freq]['t_pre']
    
    t_stims = [t_pre + ISIs[idx_freq] * i for i in np.arange(cc_APs_n_stims)]

    t_stims_df[freq] = t_stims

# add both dataframes and create spike times in continues timeseries
t_APs_cont = t_APs.add(t_stims_df)


# %%

fig_event, axs_event = plt.subplots(6,1,
                                    layout = 'constrained')

for idx_freq, freq in enumerate(frequencies):
    axs_event[idx_freq].eventplot(t_APs_cont[freq].dropna())










