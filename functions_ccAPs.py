# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 12:06:54 2024

@author: nesseler
"""

import os
import pandas as pd

from directories_win import quant_data_dir


def import_AP_measurement_all_freqs(frequencies, parameter, cell_ID):
    '''
    Functions imports excel sheet with AP parameters for one cell and the given 
    stimulation frequencies. Returns the specified parameters as a Dataframe of 
    all frequencies.
    Parameters:
        frequencies : List of frequencies to import.
        parameter : Specified parameter to be imported
        cell_ID : Cell-ID
    '''

    # set directory to APs subfolder to get parameters of action potentials
    quant_dir = os.path.join(quant_data_dir, 'APs')

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