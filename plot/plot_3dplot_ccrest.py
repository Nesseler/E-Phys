# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 14:51:30 2025

@author: nesseler
"""

import pandas as pd
import numpy as np
from os.path import join
import scipy as sc
import warnings
from tqdm import tqdm

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir, figure_dir

# custom functions
from functions.functions_useful import butter_filter, calc_time_series, calc_dvdt_padded
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol


# define protocol
PGF = 'cc_rest'

# get all cell_IDs for cc_rest
# cell_IDs_original = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = 'PGFs')
# cell_IDs_new = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = 'PGFs_Syn')
# cell_IDs = cell_IDs_original + cell_IDs_new

cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = 'PGFs_Syn')

# get number of cells
n_cells = len(cell_IDs)

# init plotting
from functions.initialize_plotting import * 


# verification plots
vplots = True


# %% data loading

# initialize dataframes to populate in the loop
v_df = pd.DataFrame(columns = cell_IDs)
vf_df = pd.DataFrame(columns = cell_IDs)
SR_df = pd.DataFrame(index = cell_IDs, columns = ['SR'])
v_mem_df = pd.DataFrame(index = cell_IDs, columns = ['v_mem'])


print('loading ...')

for cell_idx, cell_ID in enumerate(tqdm(cell_IDs)):
    
    # if cell_ID in cell_IDs_original:
    #     sheet_name = 'PGFs'
    # elif cell_ID in cell_IDs_new:
    sheet_name = 'PGFs_Syn'
            
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name = sheet_name)
    
    # get data with file path & trace index
    i, v, t, SR, n_step = get_cc_data(file_path, traceIndex, scale='s')
    
    # get first and only step of protocol
    v = v[0]
    
    # edge case when file exceeds the 30 sec recording (E-069)
    if len(v) > (30 * SR):
        warnings.warn(str(cell_ID) + ' exceeds 30 sec and will be cut down.')
        v = v[0:(30*SR)]
        
    # filter all data with 1kHz cutoff
    vf = butter_filter(v, order=3, cutoff=1e3, sampling_rate=SR)
    
    # replace first values with nans to eliminate filter artifact
    vf[:100] = np.nan

    # populate the dataframes & lists
    v_df[cell_ID] = v
    vf_df[cell_ID] = vf
    SR_df.at[cell_ID, 'SR'] = SR
    v_mem_df.at[cell_ID, 'v_mem'] = np.nanmean(vf)
    
    
    
# %% plot testing

print('generating figure ...')

t = calc_time_series(v, sampling_rate=SR, scale = 's')


# subsample dataframe
subsample_facter = 1
vf_sub = vf_df.iloc[::subsample_facter]
t_sub = t[::subsample_facter]

# sort cell_IDs by resting membrane potential
sorted_cell_IDs = v_mem_df.sort_values('v_mem').index.to_list()[::2]


fig = plt.figure(figsize= get_figure_size(),
                 layout = 'constrained')

ax = plt.axes(projection='3d')

# specify color map
cmap_str = 'flare_r'

# min max normalize time for color-code
norm = mtl.colors.Normalize(0, len(sorted_cell_IDs))

# create mappable colormap object for colorbar
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

print('plotting cells ...')

for c_idx, cell_ID in tqdm(enumerate(sorted_cell_IDs)):
    ax.plot(t_sub, vf_sub[cell_ID], c_idx,
            c = cmap.to_rgba(c_idx),
            lw = 0.25)



for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
    axis.set_pane_color((1.0, 1.0, 1.0), 0.0) 
    axis._axinfo['grid']['linewidth'] = 0.0
    axis._axinfo['grid']['color'] = "#d1d1d1"
    axis._axinfo['tick']['inward_factor'] = 0.2
    axis._axinfo['tick']['outward_factor'] = 0.0

 
ax.view_init(vertical_axis = 'y', elev = 5, azim = 7)
# ax.view_init(vertical_axis = 'y', elev = 5, azim = 25)

# (zxy)
ax.set_box_aspect(aspect=(4, 1.5, 1.5), zoom=1.5)

# x
xdict = {'ax_min' : 0,
         'ax_max' : 30,
         'pad' : 0,
         'step' : 10,
         'stepminor' : 5,
         'label' : 'Time [s]'}

# edit axis
apply_axis_settings(ax, axis = 'x', **xdict)

# y
ydict = {'ax_min' : -100,
         'ax_max' : 60,
         'pad' : 0,
         'step' : 25,
         'stepminor' : 5,
         'label' : 'Membrane potential [mV]'}

# edit axis
apply_axis_settings(ax, axis = 'y', **ydict)

# x
zdict = {'ax_min' : 0,
         'ax_max' : len(sorted_cell_IDs),
         'pad' : 0,
         'step' : 20,
         'stepminor' : 5,
         'label' : 'Cell #'}

# edit axis
apply_axis_settings(ax, axis = 'z', **zdict)

ax.zaxis.set_inverted(True)


# create saving path and save
print('saving ...')
path_figure = join(figure_dir, 'temp_figs')
save_figures(fig, f'figure-ccrest_3d-syn', path_figure, darkmode_bool, figure_format='png')

plt.show()
