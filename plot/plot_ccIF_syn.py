# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 16:21:03 2025

@author: nesseler
"""

from functions.initialize_plotting import *

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_syn_dir

# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol

# get cell_IDs
cell_IDs = get_cell_IDs_one_protocol('cc_IF', sheet_name = 'PGFs_Syn')

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# load IF
IF = pd.read_excel(join(cell_descrip_syn_dir, 'cc_IF-syn-IF.xlsx'), index_col = 'i_input')


# %% figure

# init figure and axes
fig, ax = plt.subplots(nrows = 1,
                       ncols = 1,
                       layout = 'constrained',
                       figsize = get_figure_size(width = 130, height = 100),
                       dpi = 300)

# set axis title
ax.set_title(f'sagdelta',
             fontsize=9, 
             loc='left',
             x = 0.02)

for cell_ID in cell_IDs:
    ax.plot(IF[cell_ID].dropna().index.to_list(), IF[cell_ID].dropna(),
             lw = 0.5,
             color = colors_dict['primecolor'])
    
    
    
# y
ydict = {'ax_min' : 0,
         'ax_max' : 200,
         'pad' : None,
         'step' : 10,
         'stepminor' : 2,
         'label' : 'Firing frequency [Hz]'}

apply_axis_settings(ax, axis = 'y', **ydict)

# x
xdict = {'ax_min' : -100,
         'ax_max' : 1000,
         'pad' : None,
         'step' : 100,
         'stepminor' : 25,
         'label' : 'Input current [pA]'}

apply_axis_settings(ax, axis = 'x', **xdict)


# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# align labels
fig.align_labels()

# create saving path and save
# from parameters.directories_win import tempfigs_dir
# save_figures(fig, 'cc_sag-sagdelta', tempfigs_dir, darkmode_bool, figure_format='png')

# display figure
plt.show()