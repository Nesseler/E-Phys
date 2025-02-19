# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:37:36 2025

@author: nesseler
"""

# initialize packages for script
import pandas as pd


def get_cells_list(cell_IDs = None):
    '''
    This function loads the "Cell_List_2Photons" sheet. The sheet can be 
    limited it to only the cell_IDs specified in the cell_IDs list.
    Parameters:
        cell_IDs: list of str, list of cell_ID strings like 'Exp-313', default
                  is None
    Returns:
        MetaData: pandas Dataframe, containing the metadata
    '''
    
    from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_cellsList
    
    MetaData = pd.read_excel(AMCs_cellsList,
                             sheet_name="MetaData",
                             index_col='cell_ID')

    if cell_IDs:
        MetaData = MetaData.loc[cell_IDs, :]
    
    return MetaData
