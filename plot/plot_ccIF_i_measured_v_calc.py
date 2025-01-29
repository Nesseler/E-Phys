# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 09:19:13 2025

@author: nesseler
"""
# import standard packages
from functions.initialize_packages import *

# custom functions
from functions.functions_import import get_cc_data, get_traceIndex_n_file

# PGF specific
from parameters.PGFs import cc_IF_syn_parameters
t = cc_IF_syn_parameters['t']
SR = cc_IF_syn_parameters['SR']


# set cell_ID
cell_ID = 'E-309'

# load IF protocol
# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF = 'cc_IF', cell_ID = cell_ID, sheet_name = 'PGFs_Syn')

# get data with file path & trace index
i, v, _, _, n_steps = get_cc_data(file_path, traceIndex, scale='s')

# construct current array
from functions.functions_constructors import construct_I_array
i_calc, i_steps = construct_I_array(cell_ID, n_steps,
                                    PGF = 'cc_IF',
                                    sheet_name = 'PGFs_Syn',
                                    parameters = cc_IF_syn_parameters)

# %% figure

# init plotting
from functions.initialize_plotting import *

fig, ax = plt.subplots(nrows = 1,
                       ncols = 1,
                       figsize = get_figure_size(width = 120, height = 100),
                       layout = 'constrained')

# set axis title
ax.set_title(f'{cell_ID} cc_IF',
             fontsize=9, 
             loc='left',
             x = 0.015)   

i_collection = LineCollection([np.column_stack([t, i[step]]) for step in range(n_steps)], 
                              colors = 'grey', 
                              linestyle = 'solid',
                              lw = 1, 
                              label = 'i measured') 
ax.add_collection(i_collection)

icalc_collection = LineCollection([np.column_stack([t, i_calc[step]]) for step in range(n_steps)], 
                                  colors= 'w',
                                  lw = 0.5,
                                  linestyle = 'dashed', 
                                  label = 'i calculated') 
ax.add_collection(icalc_collection)

fig.legend(fontsize = 9, frameon = False)
    
# y
ydict = {'ax_min' : -100,
         'ax_max' : 1000,
         'pad' : 10,
         'step' : 100,
         'stepminor' : 50,
         'label' : 'Input current [pA]',
         'start_at_0' : True}

apply_axis_settings(ax, axis = 'y', **ydict)

# x
xdict = {'ax_min' : 0,
         'ax_max' : 1500,
         'pad' : 10,
         'step' : 250,
         'stepminor' : 50,
         'label' : 'Time [ms]'}

apply_axis_settings(ax, axis = 'x', **xdict)

# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# create saving path and save
from parameters.directories_win import tempfigs_dir
save_figures(fig, f'{cell_ID}-current_input-record_v_calc', tempfigs_dir, darkmode_bool, figure_format='both')

plt.show()


