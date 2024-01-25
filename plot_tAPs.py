# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:27:54 2024

@author: nesseler
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtl
import seaborn as sbn


# custom directories & parameters
from PGFs import cc_APs_parameters, cc_APs_n_stims

from directories_win import cell_descrip_dir, quant_data_dir, figure_dir, vplot_dir

from functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size
from functions_export import set_df_to_cell_descrips


# %%

frequencies = list(cc_APs_parameters.keys())
frequencies.reverse()
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

# initialise dataframe to populate in loop
tAPs_df = pd.DataFrame()

# set directory to APs subfolder to get parameters of action potentials
quant_dir = os.path.join(quant_data_dir, 'APs')

# create a list of cell-IDs that have subfolder in the APs folder i.e. that have been analysed
cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]



def import_AP_measurement_all_freqs(frequencies, parameter, cell_ID):   

    # create file path
    cell_data_path = os.path.join(quant_dir, cell_ID)
    
    # get number of APs for each stim
    measurement_df = pd.DataFrame()
    
    for frequency in frequencies:
    # frequency = '1Hz'
            
        # create file path
        file_path = os.path.join(cell_data_path, f'{cell_ID}_{frequency}.xlsx')
        
        # read excel file
        AP_all_params = pd.read_excel(file_path, index_col = 'idx_step')
        
        # write to dictionary
        measurement_df[frequency] = AP_all_params[parameter]

    return measurement_df



# %% load additional data

nAPs_df = pd.read_excel(os.path.join(cell_descrip_dir, 'nAPs.xlsx'), index_col = 'frequencies')
mISIs_df = pd.read_excel(os.path.join(cell_descrip_dir, 'mean_ISIs.xlsx'), index_col = 'frequencies')
rfreq_df = pd.read_excel(os.path.join(cell_descrip_dir, 'resul_freq.xlsx'), index_col = 'frequencies')


# %% get data

# start with one cell
cell_ID = 'E-117'


mean_ISIs_df = pd.DataFrame(index = frequencies)

# for cell_ID in cell_IDs:
t_APs = import_AP_measurement_all_freqs(frequencies, 't_peaks', cell_ID)



# %% create dataframe with times that

# create list with interstimulus interval per stimulation frequency
ISIs = [sum(cc_APs_parameters[f].values()) for f in frequencies]
total_dur = [ISI * cc_APs_n_stims for ISI in ISIs]

# create dataframe with all stimuation time of every frequency
t_stims_df = pd.DataFrame()

for idx_freq, freq in enumerate(frequencies):
    
    t_pre = cc_APs_parameters[freq]['t_pre']
    
    t_stims = [t_pre + ISIs[idx_freq] * i for i in np.arange(cc_APs_n_stims)]

    t_stims_df[freq] = t_stims

# add both dataframes and create spike times in continues timeseries
t_APs_cont = t_APs.add(t_stims_df)



# %% verification plot

darkmode_bool = True

colors_dict = get_colors(darkmode_bool)

set_font_sizes()

ax_keys = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

fig_event, axs_event = plt.subplot_mosaic('AGHI;BGHI;CGHI;DGHI;EGHI;FGHI', 
                                          layout = 'constrained',
                                          figsize = get_figure_size(),
                                          gridspec_kw = {'width_ratios': [6,2,2,2]})


for idx_freq, freq in enumerate(frequencies):
    axs_event[ax_keys[idx_freq]].eventplot(t_APs_cont[freq].dropna(), colors = colors_dict['primecolor'])
    axs_event[ax_keys[idx_freq]].set_xlim([0, total_dur[idx_freq]])
    axs_event[ax_keys[idx_freq]].set_xticks(np.linspace(0, total_dur[idx_freq], 6))
    axs_event[ax_keys[idx_freq]].set_ylabel(freq)
    
    # eventplot for stimulations
    axs_event[ax_keys[idx_freq]].eventplot(t_stims_df[freq] + 5, lineoffsets = 2.25, linelength = 0.75, colors = 'grey')

fig_event.suptitle(cell_ID)

[axs_event[ax_i].set_yticks([]) for ax_i in ax_keys[:6]]

[axs_event[ax_i].spines[spine].set_visible(False) for ax_i in ax_keys[:6] for spine in ['left', 'top', 'right']]

[axs_event[ax_i].grid(False) for ax_i in ax_keys[:]]

fig_event.supxlabel('Time [ms]')


# number of APs
axs_event['G'].plot(nAPs_df[cell_ID], nAPs_df.index)

axs_event['G'].set_xlim([0, 100])
axs_event['G'].set_xlabel('Number of APs\n[#]')

# mean ISI
axs_event['H'].plot(mISIs_df[cell_ID], mISIs_df.index)

axs_event['H'].set_xlim([0, 1000])
axs_event['H'].set_xlabel('Mean ISI\n[ms]')


# resulting frequency
axs_event['I'].plot(rfreq_df[cell_ID], rfreq_df.index)

axs_event['I'].set_xlim([0, 75])
axs_event['I'].set_xlabel('Resulting freq.\n[Hz]')

# save figure as verification plot
fig_path = os.path.join(vplot_dir, 'APs', cell_ID)
if not os.path.exists(fig_path):
    os.mkdir(fig_path)
save_figures(fig_event, f'{cell_ID}_ccAPs', fig_path, darkmode_bool)

plt.show()

plt.pause(0.1)













