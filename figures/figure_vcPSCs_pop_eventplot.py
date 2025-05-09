# -*- coding: utf-8 -*-
"""
Created on Fri May  9 13:41:22 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, vplot_dir, synaptic_dir, quant_data_dir

# custom functions
from functions.functions_import import get_vc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_filter import butter_filter

import pickle
import gc

plot = False

PGF = 'vc-Erest-3min'

conditions = ['ctrl', 'adaEk']

SR = 100_000

t = np.arange(0, (6*30), 1/SR)

# get all cell_IDs
cell_IDs = get_cell_IDs_one_protocol('vc-Erest-3min-ctrl', 'PGFs_Syn')


# %% set dataframe / dicts

event_peak_locs = dict.fromkeys(cell_IDs)
amplitudes = dict.fromkeys(cell_IDs)
halfdecays = dict.fromkeys(cell_IDs)


# %% load

for cell_ID in tqdm(cell_IDs):
    # set filename
    filename = f'miniMLdetect_{cell_ID}_Erest_3min_ctrl' 
    
    # open a file, where you stored the pickled data
    file = open((synaptic_dir + '/miniML_dtc-Erest-ctrl/' + filename + '.pickle'), 'rb')
    
    # dump information to that file
    detection = pickle.load(file)
    
    # close and remove (from memory) the file
    file.close()
    del file 
    gc.collect()
    
    # write to dicts
    amplitudes[cell_ID] = detection['individual_values']['amplitudes']
    halfdecays[cell_ID] = detection['individual_values']['half_decaytimes'] *1e3
    event_peak_locs[cell_ID] = detection['event_location_parameters']['event_peak_locations'] / SR


# %%

n_events = pd.DataFrame(index = cell_IDs)
for cell_ID in cell_IDs:
    n_events.at[cell_ID, 'n_events'] = event_peak_locs[cell_ID].shape[0]

cell_IDs = n_events.sort_values(by = ['n_events']).index.to_list()

# %% plot

# init plotting
from functions.initialize_plotting import *


fig, ax = plt.subplots(nrows = 1, ncols = 1,
                       figsize = get_figure_size(width = 100, height = 80),
                       dpi = 600,
                       layout = 'constrained')

### time color coding ###

# specify color map
cmap_str = 'inferno_r'

# min max normalize time for color-code
norm = mtl.colors.Normalize(-40, 0)

# create mappable colormap object for colorbar
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)



for ci, cell_ID in enumerate(cell_IDs):
    p = event_peak_locs[cell_ID]
    ll = [np.abs(amplitudes[cell_ID])/60]
    lo = [[ci-(ll[0][i]/2) for i in range(len(p))]]
    lw = halfdecays[cell_ID]/30
    
    # print(np.max(l))
    
    ax.eventplot(positions = p,
                 lineoffsets = ci,
                 lw = 0.25,
                 linelengths = 0.8,
                 color = [cmap.to_rgba(amplitudes[cell_ID])])
    
    
plt.yticks(ticks = np.arange(0, len(cell_IDs)),
           labels = cell_IDs)
plt.show()


