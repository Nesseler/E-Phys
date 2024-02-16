# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 14:55:42 2024

@author: nesseler
"""


import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtl
import seaborn as sbn

from os.path import join

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir, figure_dir, quant_data_dir
from parameters.PGFs import cc_APs_parameters
from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size


frequencies = cc_APs_parameters.keys()
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

resul_freq_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccAPs-resul_freq.xlsx'), index_col = 'frequencies')
cell_IDs = list(resul_freq_df.columns)

# %%



AP_parameter = 'v_amplitude'


#for cell_ID in ['E-092']:
for cell_ID in cell_IDs:
    
    test = []
    
    cell_APs_path = join(quant_data_dir, 'APs', cell_ID)
    
    
    for frequency in frequencies:
    
        AP_params = pd.read_excel(join(cell_APs_path, f'{cell_ID}_{frequency}.xlsx'))
        
        step_idx_last_AP = AP_params.query('v_amplitude.notnull()').iloc[-4:-1].index.to_list()

        mean_last_x = AP_params.loc[step_idx_last_AP, AP_parameter].mean()

        print(cell_ID, frequency, ((mean_last_x - AP_params.at[0, AP_parameter]) / AP_params.at[0, AP_parameter]))

        
        test.append(((mean_last_x - AP_params.at[0, AP_parameter]) / AP_params.at[0, AP_parameter]))
        

    plt.plot(freqs_int, test, 'x-k', alpha = 0.5)
plt.show()

