# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:41:05 2024

@author: nesseler
"""
 
import pandas as pd
import os

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir, cell_descrip_file


# test variables
# n_APs_path = os.path.join(cell_descrip_dir, 'n_APs.xlsx')
# add_df = pd.read_excel(n_APs_path, index_col = 'frequency')
# add_header_ext = 'nAPs'

# function to add additional parameter to cell descriptive values
def set_df_to_cell_descrips(add_df, add_header_ext = ''):
    '''
    Function gets dataframe that is to be added to the excel sheet containing
    cell descriptive values. Additional extension to the header can be parsed 
    if header alone is not self-explanatory.
    Parameters:
        add_df : Dataframe that is to be added.
        add_header_ext : String extension that will be added to the exsisting header.
                          Default is '' (empty).
    Returns:
        nothing, excel file is changed.
    '''
    
    # get cell descriptor table
    cells_df = pd.read_excel(cell_descrip_file, index_col = 'cell_ID')
    
    # test if index or header are cell_IDs
    if add_df.columns[0][:2] == 'E-':
        add_df = add_df.transpose()
    
    if add_header_ext != '':
        # rename header elements when add_header_ext is parsed
        rename_dict = {col_str: str(col_str) + '-' + add_header_ext for i, col_str in enumerate(add_df.columns)}
        add_df = add_df.rename(columns = rename_dict)       
        
    # test if data is already contained in dataframe
    if all(x in cells_df.columns.to_list() for x in add_df.columns.to_list()):
        cells_df.update(add_df)
        print('cell_descrips.xlsx has been updated')
    else:
        cells_df = pd.concat([cells_df, add_df], axis = 1)
        print('cell_descrips.xlsx has been extended')
    
    # save dataframe to excel file 
    cells_df.to_excel(cell_descrip_file, index_label='cell_ID')




