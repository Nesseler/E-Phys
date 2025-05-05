# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 19:02:06 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir

# custom functions
from functions.functions_import import get_vc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_filter import butter_filter


PGF = 'vc-Erest-3min'

conditions = ['ctrl', 'adaEk']

SR = 100000

t = np.arange(0, (6*30), 1/SR)

# get cell_IDs
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF + '-' + 'adaEk', sheet_name= 'PGFs_Syn')


# %% 

# initialize dicts and dataframes
i_s = {'ctrl'  : pd.DataFrame(columns=cell_IDs, index = t),
       'adaEk' : pd.DataFrame(columns=cell_IDs, index = t)}

print('loading ...')

for cell_ID in tqdm(cell_IDs):
    
    for condition in conditions:

        # get the traceIndex and the file path string for data import functions
        traceIndex, file_path = get_traceIndex_n_file(PGF + '-' + condition, cell_ID, sheet_name = 'PGFs_Syn')
        
        # get data with file path & trace index
        i, _, _, _, n_steps = get_vc_data(file_path, traceIndex, scale='ms')
        
        # concatenate steps
        i = i.flatten('C')
        
        # filter
        # filter all data with 1kHz cutoff
        i = butter_filter(i, order=3, cutoff=1e3, sampling_rate=SR)
        
        # replace first values with nans to eliminate filter artifact
        i[:100] = np.nan
        
        # write to dataframe
        i_s[condition][cell_ID] = i


# %%

# init plotting
from functions.initialize_plotting import *


fig, axs = plt.subplots(nrows = len(cell_IDs),
                        ncols = 2,
                        figsize = get_figure_size(),
                        dpi = 300,
                        layout = 'constrained',
                        sharex=True,
                        sharey=True)

fig.suptitle('synaptic currents')

    
axs[0][0].set_title('$E_k$ = -98 mV')
axs[0][1].set_title('$E_k$ = -85 mV')


for cell_i, cell_ID in enumerate(cell_IDs):
    
    for con_i, condition in enumerate(conditions):
        
        axs[cell_i][con_i].plot(t, 
                                i_s[condition][cell_ID],
                                c = colors_dict['primecolor'],
                                lw = 0.25)

# y
ydict = {'ax_min' : -100,
         'ax_max' : 10,
         'pad' : 1.1,
         'step' : 100,
         'stepminor' : 5,
         'label' : ''}

# edit axis
for col in axs:
    for ax in col:
        apply_axis_settings(ax, axis = 'y', **ydict)

# x
xdict = {'ax_min' : 0,
         'ax_max' : 180,
         'pad' : 3,
         'step' : 60,
         'stepminor' : 5,
         'label' : ''}

# edit axis
for col in axs:
    for ax in col:
        apply_axis_settings(ax, axis = 'x', **xdict)

# set sup labels
fig.supylabel('Current [pA]')
fig.supxlabel('Time [s]')

# remove spines
[ax.spines[spine].set_visible(False) for col in axs for ax in col for spine in ['top', 'right']]

# align labels
fig.align_labels()

# create saving path and save
fig_path = join(figure_dir, 'temp_figs')
save_figures(fig, f'adaEK-pre_post_PSCs', fig_path, darkmode_bool, figure_format='both')

# display figure
plt.show()










# %%

import pickle
import gc


# import directories 
from parameters.directories_win import synaptic_dir
# miniML_path = synaptic_dir + '/miniML_dtc-validation'

# set cell_ID
cell_ID = 'E-315'
treatments = ['ctrl', 'adaEk']

# set dicts
traces = dict.fromkeys(['ctrl', 'adaEk'])
amplitudes = dict.fromkeys(['ctrl', 'adaEk'])
amplitudes_hist = dict.fromkeys(['ctrl', 'adaEk'])
events = dict.fromkeys(['ctrl', 'adaEk'])
event_locations = dict.fromkeys(['ctrl', 'adaEk'])

# set histogram bins
bins = np.arange(-50, 0+1, 1)

for treatment in ['ctrl', 'adaEk']:

    # set filename
    filename = f'miniMLdetect_{cell_ID}_Erest_3min_{treatment}' 
    
    # open a file, where you stored the pickled data
    file = open((synaptic_dir + f'/miniML_dtc-Erest-{treatment}/' + filename + '.pickle'), 'rb')
    
    # dump information to that file
    detection = pickle.load(file)
    
    # close and remove (from memory) the file
    file.close()
    del(file)
    gc.collect()
    
    # write to dicts
    traces[treatment] = detection['mini_trace']
    events[treatment] = detection['events']
    amplitudes[treatment] = detection['individual_values']['amplitudes']
    event_locations[treatment] = detection['event_location_parameters']['event_peak_locations']
    
    # calc occurrances
    amplitudes_hist[treatment], _ = np.histogram(a = detection['individual_values']['amplitudes'], bins = bins)


# %%

# # set data and bins
# data = amplitudes[treatment]

# hist, _ = np.histogram(a = amplitudes['ctrl'], bins = bins)

# init figure
fig, ax = plt.subplots(1,1, dpi = 300, layout = 'constrained')

# set title
fig.suptitle(f'{cell_ID} amplitude histogram')
# ax.hist(amplitudes['ctrl'], bins=bins, 
#         histtype='step', 
#         color='k',
#         lw = 1)

# for treatment in treatments:
# plot ctrl
ax.stairs(amplitudes_hist['ctrl'], bins, 
          fill = False,
          lw = 1,
          color = colors_dict['primecolor'],
          label = 'ctrl')

# plot treatment
ax.stairs(amplitudes_hist['adaEk'], bins, 
          fill = False,
          lw = 1,
          color = 'goldenrod',
          label = 'adaEk')
    
ax.legend(title = 'Treatment', 
          frameon = False,
          loc = 'upper left',
          fontsize = 9)

# edit axis
ax.set_ylabel('Event count [#]')
ax.set_ylim([-2, 132])
ax.set_yticks(ticks = np.arange(0, 130+1, 20))
ax.set_yticks(ticks = np.arange(0, 130+1, 2), minor = True)
ax.spines['left'].set_bounds([0, 130])

ax.set_xlabel('Amplitude [pA]')
ax.set_xticks(ticks = np.arange(-50, 0+1, 10))
ax.set_xticks(ticks = np.arange(-50, 0, 1), minor = True)
ax.set_xlim([-51, 1])
ax.spines['bottom'].set_bounds([-50, 0])

# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# display figure
plt.show()


# %%

# init figure
fig, axs = plt.subplot_mosaic(mosaic = [[0, 0, 2, 2], [1, 1, 3, 3], [4, 4, 4, 5], [4, 4, 4, 6]],#'AACC;BBDD;EEEF;EEEG', 
                             dpi = 300, 
                             layout = 'constrained',
                             figsize = get_figure_size(width = 159.2, height = 159.2))




# set histogram axis
ax = axs[4]

# title
ax.set_title('E: Amplitude histogram',
             loc = 'left')

ax.stairs(amplitudes_hist['ctrl'], bins, 
          fill = False,
          lw = 1,
          color = colors_dict['primecolor'],
          label = 'ctrl')

# plot treatment
ax.stairs(amplitudes_hist['adaEk'], bins, 
          fill = False,
          lw = 1,
          color = 'goldenrod',
          label = 'adaEk')
    
ax.legend(title = 'Treatment', 
          frameon = False,
          loc = 'upper left',
          fontsize = 9)

# edit axis
ax.set_ylabel('Event count [#]')
ax.set_ylim([-2, 132])
ax.set_yticks(ticks = np.arange(0, 130+1, 20))
ax.set_yticks(ticks = np.arange(0, 130+1, 5), minor = True)
ax.spines['left'].set_bounds([0, 130])

ax.set_xlabel('Amplitude [pA]')
ax.set_xticks(ticks = np.arange(-50, 0+1, 10))
ax.set_xticks(ticks = np.arange(-50, 0, 1), minor = True)
ax.set_xlim([-51, 1])
ax.spines['bottom'].set_bounds([-50, 0])





# remove spines
[axs[i].spines[spine].set_visible(False) for i in range(7) for spine in ['top', 'right']]

# display figure
plt.show()
