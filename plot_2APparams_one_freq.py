# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 21:44:19 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import directories_win as directories
from PGFs import cc_APs_parameters, cc_th1Ap_parameters
import pandas as pd
from useful_functions import calc_time_series, butter_filter, calc_dvdt
import os
from cc_IF_functions import get_IF_data
from plotting_functions import get_colors, save_figures, get_figure_size, set_font_sizes, return_segments
import scipy as sc
import parameters

from matplotlib.collections import LineCollection

from spiketrains_functions import get_colorcode


# %%

frequencies = list(cc_APs_parameters.keys())

freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

n_APs_df = pd.DataFrame()

quant_dir = os.path.join(directories.quant_data_dir, 'APs')

cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]



# %%

frequency = '30Hz'

freq_int = int(frequency.replace('Hz', ''))
freq_idx = frequencies.index(frequency)


def import_AP_measurement_one_cell(frequency, parameter, cell_ID):

    # get number of APs for each stim
    measurement_df = pd.DataFrame()
    
    # create file path
    file_path = os.path.join(cell_data_path, f'{cell_ID}_{frequency}.xlsx')
    
    # read excel file
    AP_all_params = pd.read_excel(file_path, index_col = 'idx_step')
    
    # write to dictionary
    measurement_df[frequency] = AP_all_params[parameter]

    return measurement_df


FWHM_df = pd.DataFrame()
ampl_df = pd.DataFrame()

for cell_ID in cell_IDs:
    
    # create file path
    cell_data_path = os.path.join(os.path.join(directories.quant_data_dir, 'APs'), cell_ID)

    cur_FWHM_df = import_AP_measurement_one_cell(frequency, 'FWHM', cell_ID)
    cur_ampl_df = import_AP_measurement_one_cell(frequency, 'v_amplitude', cell_ID)

    FWHM_df[cell_ID] = cur_FWHM_df
    ampl_df[cell_ID] = cur_ampl_df

# FWHM_df = import_AP_measurement_all_cells(frequencies, 'FWHM')
# ampl_df = import_AP_measurement_all_freqs(frequencies, 'v_amplitude')


# %% interpolate nan values


FWHM_df_i = FWHM_df.interpolate()
ampl_df_i = ampl_df.interpolate()

# %%

darkmode_bool = True

colors_dict = get_colors(darkmode_bool)

fig_APps, axs_APps = plt.subplots(nrows = 2, 
                                  ncols = 1,
                                  layout = 'constrained',
                                  dpi = 600,
                                  figsize = get_figure_size(),
                                  sharex = 'col')

set_font_sizes()


norm_min = 1
norm_max = 75
cmap_str = 'plasma'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
cmap.set_array(freqs_int)


fig_APps.colorbar(cmap, 
                  ax = axs_APps,
                  ticks = freqs_int,
                  label = 'Stimulation frequency [Hz]')


for idx, cell_ID in enumerate(cell_IDs):
    axs_APps[0].plot(FWHM_df_i[cell_ID].index, FWHM_df_i[cell_ID],
                     c = cmap.to_rgba(freqs_int[freq_idx]))
    
    axs_APps[0].scatter(FWHM_df[cell_ID].index, FWHM_df[cell_ID],
                        marker = '.',
                        s = 75,
                        color = cmap.to_rgba(freqs_int[freq_idx]))
    
    axs_APps[1].plot(ampl_df_i[cell_ID].index, ampl_df_i[cell_ID],
                      c = cmap.to_rgba(freqs_int[freq_idx]))
    
    axs_APps[1].scatter(ampl_df[cell_ID].index, ampl_df[cell_ID],
                        marker = '.',
                        s = 75,
                        color = cmap.to_rgba(freqs_int[freq_idx]))



# FWHM y axis 
axs_APps[0].set_ylim([0, 6])
axs_APps[0].set_ylabel('FWHM [ms]')

axs_APps[0].set_yticks(ticks = np.linspace(0, 6, 7))


# ampli y axis 
axs_APps[1].set_ylim([15, 125])
axs_APps[1].set_ylabel('Amplitude [mV]')

axs_APps[1].set_yticks(ticks = np.linspace(20, 120, 11))


# x axis
axs_APps[1].set_xlim([0, 99])
axs_APps[1].set_xlabel('n-th stimulation [#]')
axs_APps[1].set_xticks(ticks = np.arange(19, 99 + 1, 20),
                       labels = np.arange(20, 100 + 1, 20))
axs_APps[1].set_xticks(ticks = np.arange(0, 99 + 1, 1),
                       minor = True)


[ax.grid(False) for ax in axs_APps]


plt.show()


save_figures(fig_APps, f'FWHM_ampli_all_cells_{frequency}', directories.figure_dir, darkmode_bool)


