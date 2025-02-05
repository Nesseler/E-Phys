# -*- coding: utf-8 -*-
"""
Created on Wed Feb 02 14:40:13 2025

@author: nesseler
"""
 
import pandas as pd
from os.path import join


def write_exportvars_to_excel(export_vars, export_prefix):
    '''
    This function takes export variables in the provided dictionary and
    writes them to the corresponding excel files. If the file does not exist, 
    a new version will be created. If that file does exist, the new values will
    be concatnated if the values are not already included or written in the same
    position if the values are alreay included.
    Parameters:
        export_vars: dict, containing the pandas dataframes for export and the
                     name of their file as key.
        export_prefix: str, prefix of filenames
    '''

    # set filename extension to excel fileformat
    export_extension = '.xlsx'
    
    from parameters.directories_win import cell_descrip_syn_dir
    
    for export_name, export_var in export_vars.items():
        
        # get index label
        index_label = export_var.index.name
        
        # check for non unique cols
        cols = export_var.nunique()[export_var.nunique() > 0].index.to_list()
        
        # use these columns to check for rows
        if index_label == 'cell_ID':
            rows = export_var.nunique(axis = 1)[export_var.nunique(axis = 1) > 0].index.to_list()
            concat_axis = 0
        else:
            rows = export_var.index.to_list()
            concat_axis = 1
    
        # try loading and writing or create new file
        try:
            loaded_export_var = pd.read_excel(join(cell_descrip_syn_dir, export_prefix + export_name + export_extension),
                                              index_col = index_label)
          
            # break
            try:
                # combine both dataframes
                loaded_export_var.loc[rows, cols] = export_var.loc[rows, cols].values
    
            # if values not yet in file, append/concat new values
            except KeyError:
                loaded_export_var = pd.concat([loaded_export_var, export_var.loc[rows, cols]], axis = concat_axis)
                    
            # save activity dataframe
            loaded_export_var.to_excel(join(cell_descrip_syn_dir, export_prefix + export_name + export_extension), 
                                        index_label=index_label)
        
        except FileNotFoundError:
            # save activity dataframe
            export_var.to_excel(join(cell_descrip_syn_dir, export_prefix + export_name + export_extension), 
                                index_label=index_label)



