# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 11:24:22 2024

@author: nesseler
"""

import pandas as pd
import os

# custom directories & parameters
from directories_win import quant_data_dir, cell_descrip_dir
from PGFs import cc_APs_parameters

from functions_ccAPs import import_AP_measurement_all_freqs
from functions_export import set_df_to_cell_descrips


# %%

frequencies = list(cc_APs_parameters.keys())

# set directory to APs subfolder to get parameters of action potentials
quant_dir = os.path.join(quant_data_dir, 'APs')

# create a list of cell-IDs that have subfolder in the APs folder i.e. that have been analysed
cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]



# %%

parameters = ['t_peaks', 'FWHM', 'v_amplitude']

# cell_ID = 'E-092'

for parameter in parameters:

    # initialise dataframe to populate in loop
    mean_cols = ['mean' + '-' + parameter + '-' + freq for freq in frequencies]
    std_cols = ['std' + '-' + parameter + '-' + freq for freq in frequencies]
    APparam_df = pd.DataFrame(columns = mean_cols + std_cols)
    APparam_df.index.name = 'cell_ID'
    
    for cell_ID in cell_IDs:
        # get parameter for all frequencies as dataframe
        param_df = import_AP_measurement_all_freqs(frequencies, parameter, cell_ID)
    
        # calc mean and std and set to the dataframe at appropriate locations
        APparam_df.loc[cell_ID, mean_cols] = param_df.mean(axis = 0).to_list()
        APparam_df.loc[cell_ID, std_cols] = param_df.std(axis = 0).to_list()
    
    
    
    set_df_to_cell_descrips(APparam_df)
    
    APparam_df.to_excel(os.path.join(cell_descrip_dir, f'mean_std-{parameter}-ccAPs.xlsx'), index=cell_ID)

























