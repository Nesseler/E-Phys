# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:43:42 2025

@author: nesseler
"""

from functions.initialize_plotting import *

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import quant_data_dir

# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol

# get cell_IDs
cell_IDs = get_cell_IDs_one_protocol('cc_sag', sheet_name = 'PGFs_Syn')

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# limit cell_IDs
cell_IDs = MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list()

# %% figure

# init figure and axes
fig, ax = plt.subplots(nrows = 1,
                        ncols = 1,
                        layout = 'constrained',
                        figsize = get_figure_size(width = 100, height = 100),
                        dpi = 300)

   
# set axis title
ax.set_title(f'sagdelta',
             fontsize=9, 
             loc='left',
             x = 0.02)

# load files
for cell_ID in cell_IDs:
    
    # load
    sagdeltas_cell = pd.read_excel(join(quant_data_dir, 'cc_sag-sagdeltas', f'{cell_ID}-sagdeltas.xlsx'))
    
    # get sagdeltas
    v_mins = sagdeltas_cell['v_min']
    sagdeltas = sagdeltas_cell['sagdelta']
    
    if sagdeltas.tail(1).values > 30:
        v_mins.iat[-1] = np.nan
        sagdeltas.iat[-1] = np.nan
        
    if v_mins.tail(1).values > -120:
        print(cell_ID)
    
    plt.plot(sagdeltas, v_mins,
             lw = 0.5,
             color = region_colors[MetaData.at[cell_ID, 'Region']])
    
# y
ydict = {'ax_min' : -160,
         'ax_max' : -80,
         'pad' : None,
         'step' : 20,
         'stepminor' : 5,
         'label' : 'Membrane potential [mV]'}

apply_axis_settings(ax, axis = 'y', **ydict)

# x
xdict = {'ax_min' : 0,
         'ax_max' : 25,
         'pad' : None,
         'step' : 5,
         'stepminor' : 1,
         'label' : 'Sag delta [mV]'}

apply_axis_settings(ax, axis = 'x', **xdict)


# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# align labels
fig.align_labels()

# create saving path and save
from parameters.directories_win import tempfigs_dir
save_figures(fig, 'cc_sag-sagdelta', tempfigs_dir, darkmode_bool, figure_format='png')

# display figure
plt.show()

    

    
    
