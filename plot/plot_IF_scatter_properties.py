# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 18:06:06 2024

@author: nesseler
"""

import pandas as pd
from os.path import join
import seaborn as sbn
import matplotlib.pyplot as plt
import numpy as np

from functions.functions_plotting import set_font_sizes, get_colors, get_figure_size, save_figures

from parameters.directories_win import quant_data_dir, cell_descrip_dir, table_file, figure_dir


active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')
passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col = 'cell_ID')

sbn.jointplot(data = passiv_properties_df,
                x = 'r_input',
                y = 'tau_mem')












