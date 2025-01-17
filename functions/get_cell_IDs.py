# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:32:07 2024

@author: nesseler
"""

import pandas as pd

from parameters.directories_win import table_file
from parameters.PGFs import cc_APs_parameters


def get_cell_IDs_all_ccAPfreqs():
    '''
    This function returns a list of cell IDs that were stimulated with all 
    frequencies (1 - 75 Hz) in the cc_APs protocols.
    Parameter:
        -
    Returns:
        cell_IDs : List of cell_IDs
    '''
    # get table
    table = pd.read_excel(table_file, sheet_name="PGFs", index_col='cell_ID')
    
    # loop to create string to include all frequencies in query
    query_str = ''
    
    frequencies = list(cc_APs_parameters.keys())
    
    for idx, frequency in enumerate(frequencies):
        PGF = 'cc_APs_' + frequency
        
        if idx > 0:
            query_str = query_str + ' and '
            
        query_str = query_str + f'{PGF}.notnull()'
    
    query_str = query_str + ' and cc_th1AP.notnull()'    
    
    # limit lookup table
    lookup_table = table.query(query_str)
    
    # cell IDs 
    cell_IDs = list(lookup_table.index)

    return cell_IDs



def get_cell_IDs_one_protocol(PGF = 'cc_IF', sheet_name = 'PGFs'):
    '''
    Functions returns list of cell_IDs that contain the provided protocol.
    Parameters:
        PGF : str of protocol. Default is cc_IF.
    Returns:
        cell_IDs : List of cell_IDs
    '''
    
    # database as table
    table = pd.read_excel(table_file, sheet_name=sheet_name, index_col='cell_ID')

    # limit 
    lookup_table = table.query(f'{PGF}.notnull()')

    # cell IDs 
    cell_IDs = list(lookup_table.index)
    
    return cell_IDs
    
    
    
    
    
    
    
    