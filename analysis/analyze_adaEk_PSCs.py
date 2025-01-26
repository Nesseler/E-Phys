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
from functions.functions_useful import butter_filter


PGF = 'vc_Erest_3min'

conditions = ['ctrl', 'adaEk']

SR = 100000

t = np.arange(0, (6*30), 1/SR)

# get cell_IDs
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF + '_' + 'adaEk', sheet_name= 'PGFs_Syn')


# %% 

# initialize dicts and dataframes
i_s = {'ctrl'  : pd.DataFrame(columns=cell_IDs, index = t),
       'adaEk' : pd.DataFrame(columns=cell_IDs, index = t)}

print('loading ...')

for cell_ID in tqdm(cell_IDs):
    
    for condition in conditions:

        # get the traceIndex and the file path string for data import functions
        traceIndex, file_path = get_traceIndex_n_file(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')
        
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






















# # define min PSC width for detection
# PSC_min_width = 1 * SR/1e3 # ms

# # find peaks
# idc_spikes, dict_peak = sc.signal.find_peaks(-i, 
#                                              prominence = 10, 
#                                              width = PSC_min_width)

# plt.plot(i, lw = 0.25)
# plt.scatter(idc_spikes, i[idc_spikes],
#             marker = 'x',
#             c = 'r')
# plt.show()


# # %%

# for ip, idx_peak in enumerate(idc_spikes):
    
#     plt.plot(i[idx_peak-750:idx_peak+2250], lw = 0.25)
#     plt.title(dict_peak['prominences'][ip])
#     # plt.show()