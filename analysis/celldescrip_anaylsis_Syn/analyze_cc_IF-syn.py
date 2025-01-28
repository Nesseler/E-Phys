# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:27:42 2025

@author: nesseler
"""
# import standard packages
from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_syn_dir
from parameters.parameters import min_peak_prominence, min_peak_distance

# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_import import get_cc_data, get_traceIndex_n_file

# PGF specific
from parameters.PGFs import cc_IF_syn_parameters
t = cc_IF_syn_parameters['t']
SR = cc_IF_syn_parameters['SR']
    
# define protocol
PGF = 'cc_IF'
sheet_name = 'PGFs_Syn'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = sheet_name)

# get number of cells
n_cells = len(cell_IDs)

# init plotting
from functions.initialize_plotting import *


cell_ID = 'E-303'


# load IF protocol
# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name = sheet_name)

# get data with file path & trace index
i, v, _, _, n_steps = get_cc_data(file_path, traceIndex, scale='s')


# %%

from functions.functions_constructors import construct_I_array

# construct current array
i_calc, i_steps = construct_I_array(cell_ID, n_steps,
                                    PGF = PGF,
                                    sheet_name = sheet_name,
                                    parameters = cc_IF_syn_parameters)



for step in range(n_steps):
    plt.plot(i[step], 'grey', lw=0.5)
    
    plt.plot(i_calc[step], 'w', lw=1)