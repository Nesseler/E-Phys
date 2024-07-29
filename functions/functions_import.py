# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 18:57:17 2024

@author: nesseler
"""
import pandas as pd
import os

from parameters.directories_win import table_dir, table_file, raw_data_dir

def get_traceIndex_n_file(PGF = 'ccth1AP', cell_ID = 'E-092'):
    # excel sheet with PGF indices as lookup table
    lookup_table = pd.read_excel(table_file,
                                 sheet_name="PGFs",
                                 index_col='cell_ID')
    
    # get indices of current cell with the dataframe containing all indices    
    group_idx = int(lookup_table.at[cell_ID, 'group'])-1
    series_idx = int(lookup_table.at[cell_ID, f'{PGF}'])-1

    # construct traceIndex with indices
    traceIndex = [group_idx, series_idx, 0, 0]

    # call on data file with indices from dataframe above
    current_file = lookup_table.at[cell_ID, 'file']

    data_file_path = os.path.join(raw_data_dir, current_file + '.dat')

    data_file_path_str = fr"{data_file_path}"
    
    return traceIndex, data_file_path_str




def get_onlyfiles_list(dir_path):

    from os import listdir
    from os.path import isfile, join  

    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]