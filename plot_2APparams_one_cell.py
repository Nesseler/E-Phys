# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 20:36:45 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import pandas as pd
import os

# custom directories & parameters
import directories_win as directories
from PGFs import cc_APs_parameters

from functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes


# %%

frequencies = list(cc_APs_parameters.keys())

freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

n_APs_df = pd.DataFrame()

quant_dir = os.path.join(directories.quant_data_dir, 'APs')

cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]



# %%

cell_ID = 'E-092'

    
# create file path
cell_data_path = os.path.join(os.path.join(directories.quant_data_dir, 'APs'), cell_ID)


def import_AP_measurement_all_freqs(frequencies, parameter):

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

FWHM_df = import_AP_measurement_all_freqs(frequencies, 't_peaks')
ampl_df = import_AP_measurement_all_freqs(frequencies, 'v_amplitude')


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


for idx, frequency in enumerate(frequencies):
    axs_APps[0].plot(FWHM_df_i[frequency].index, FWHM_df_i[frequency],
                     c = cmap.to_rgba(freqs_int[idx]))
    
    axs_APps[0].scatter(FWHM_df[frequency].index, FWHM_df[frequency],
                        marker = '.',
                        s = 75,
                        color = cmap.to_rgba(freqs_int[idx]))
    
    axs_APps[1].plot(ampl_df_i[frequency].index, ampl_df_i[frequency],
                     c = cmap.to_rgba(freqs_int[idx]))
    
    axs_APps[1].scatter(ampl_df[frequency].index, ampl_df[frequency],
                        marker = '.',
                        s = 75,
                        color = cmap.to_rgba(freqs_int[idx]))



# FWHM y axis 
# axs_APps[0].set_ylim([1.0, 2.5])
# axs_APps[0].set_ylabel('FWHM [ms]')

# axs_APps[0].set_yticks(ticks = np.linspace(1,2.5, 4))


# ampli y axis 
axs_APps[1].set_ylim([50, 110])
axs_APps[1].set_ylabel('Amplitude [mV]')

axs_APps[1].set_yticks(ticks = np.linspace(50, 110, 7))


# x axis
axs_APps[1].set_xlim([0, 99])
axs_APps[1].set_xlabel('n-th stimulation [#]')
axs_APps[1].set_xticks(ticks = np.arange(19, 99 + 1, 20),
                       labels = np.arange(20, 100 + 1, 20))
axs_APps[1].set_xticks(ticks = np.arange(0, 99 + 1, 1),
                       minor = True)


[ax.grid(False) for ax in axs_APps]


plt.show()


# save_figures(fig_APps, f'FWHM_ampli_one_cell_{cell_ID}', directories.figure_dir, darkmode_bool)


