# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 19:16:44 2024

@author: nesseler
"""
import pandas as pd
from cc_APs_onecell import export_all_freqs_and_AP_parameters
import directories_win as directories
from PGFs import cc_APs_parameters

# %%

table = pd.read_excel(directories.table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')


frequencies = list(cc_APs_parameters.keys())

# loop to create string to include all frequencies in query
query_str = ''

for idx, frequency in enumerate(frequencies):
    PGF = 'cc_APs_' + frequency
    
    if idx > 0:
        query_str = query_str + ' and '
        
    query_str = query_str + f'{PGF}.notnull()'
    

# limit lookup table
lookup_table = table.query(query_str)


# %% 

for cell_ID in lookup_table.index:
    print(f'Started: {cell_ID}')
    export_all_freqs_and_AP_parameters(cell_ID, lookup_table)
    print(f'Done: {cell_ID}')
