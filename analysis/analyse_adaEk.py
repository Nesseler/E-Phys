#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 14:47:38 2025

@author: moritznesseler
"""

import pandas as pd
import scipy as sc
import numpy as np
from os.path import join

from parameters.directories_win import raw_data_dir, table_file

from functions.functions_import import get_cc_data, get_traceIndex_n_file

from functions.get_cell_IDs import get_cell_IDs_one_protocol



PGF = 'cc_rest_adaEk'

# get cell_IDs
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF,
                                     sheet_name = 'PGFs_Syn')

cell_ID = 'E-311'

# define sampling rate
SR = 50000

# define dataframes
pre_df = pd.DataFrame(columns=cell_IDs, index=np.arange(0, ))

conditions = ['pre', 'washin', 'test']  

# get trace indeex and file path
traceIndex, file_path = get_traceIndex_n_file(PGF = PGF + '_' + conditions[0], 
                                              cell_ID = cell_ID,
                                              sheet_name = 'PGFs_Syn')



if traceIndex is list:
    
    for traceI in traceIndex:
        
        # get data
        i, v, t, SR = get_cc_data(file_path, traceIndex, 's')
    
    
    
# print(len(tr1aceIndex[1]))




## load cc_rest

# pre

# washin

# post



# vplot timeline 

