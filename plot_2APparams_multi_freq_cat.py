# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 22:27:24 2024

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
import seaborn as sbn


# %%

frequencies = list(cc_APs_parameters.keys())

freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

n_APs_df = pd.DataFrame()

quant_dir = os.path.join(directories.quant_data_dir, 'APs')

cell_IDs = [name for name in os.listdir(quant_dir) if os.path.isdir(os.path.join(quant_dir, name))]



# %%

FWHM_freq = pd.DataFrame()
ampl_freq = pd.DataFrame()

for frequency in frequencies:
 
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
    
    
    
    FWHM_melt = pd.melt(FWHM_df, var_name='cell_ID', value_name='FWHM') 
    ampl_melt = pd.melt(ampl_df, var_name='cell_ID', value_name='v_amplitude')
    
    FWHM_freq[frequency] = FWHM_melt['FWHM']
    ampl_freq[frequency] = ampl_melt['v_amplitude']


    


 # %%
 
FWHM_freq_melt = pd.melt(FWHM_freq, var_name='frequency', value_name='FWHM')
ampl_freq_melt = pd.melt(ampl_freq, var_name='frequency', value_name='v_amplitude')


# %% color code

norm_min = 1
norm_max = 75
cmap_str = 'plasma'

norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
cmap.set_array(freqs_int)


# creating specific color pallet for seaborn plotting functions
color_pal = {'1Hz'  : cmap.to_rgba(freqs_int[0]), 
             '5Hz'  : cmap.to_rgba(freqs_int[1]),
             '10Hz' : cmap.to_rgba(freqs_int[2]),
             '30Hz' : cmap.to_rgba(freqs_int[3]),
             '50Hz' : cmap.to_rgba(freqs_int[4]),
             '75Hz' : cmap.to_rgba(freqs_int[5])}

# %% plotting


darkmode_bool = True

colors_dict = get_colors(darkmode_bool)

fig_APcats, axs_APcats = plt.subplots(nrows = 2, 
                                      ncols = 1,
                                      layout = 'constrained',
                                      dpi = 600,
                                      figsize = get_figure_size(),
                                      sharex = 'col')

set_font_sizes()


fig_APcats.colorbar(cmap, 
                  ax = axs_APcats,
                  ticks = freqs_int,
                  label = 'Stimulation frequency [Hz]')


violins1 = sbn.violinplot(data = FWHM_freq_melt, 
                          x = "frequency", 
                          y = "FWHM", 
                          inner = 'quart', 
                          ax = axs_APcats[0], 
                          linewidth = 1,
                          palette = color_pal, 
                          width = 1)

for l in violins1.lines:
    l.set_color(colors_dict['primecolor'])

[violin.set_edgecolor(colors_dict['primecolor']) for violin in violins1.collections]

axs_APcats[0].set(xlabel=None)

# swarms = sbn.swarmplot(data = FWHM_freq_melt, 
#                         x = "frequency", 
#                         y = "FWHM", 
#                         size = 0.5, 
#                         ax = axs_APcats[0], 
#                         color = colors_dict['primecolor'])

violins2 = sbn.violinplot(data = ampl_freq_melt, 
                          x = "frequency", 
                          y = "v_amplitude", 
                          inner = 'quart', 
                          ax = axs_APcats[1], 
                          linewidth = 1,
                          palette = color_pal)

for l in violins2.lines:
    l.set_color(colors_dict['primecolor'])

[violin.set_edgecolor(colors_dict['primecolor']) for violin in violins2.collections]


# FWHM y axis 
axs_APcats[0].set_ylim([0, 6])
axs_APcats[0].set_ylabel('FWHM [ms]')

axs_APcats[0].set_yticks(ticks = np.linspace(0, 6, 7))


# ampli y axis 
axs_APcats[1].set_ylim([15, 125])
axs_APcats[1].set_ylabel('Amplitude [mV]')

axs_APcats[1].set_yticks(ticks = np.linspace(20, 120, 6))

[ax.grid(False) for ax in axs_APcats]


plt.show()



save_figures(fig_APcats, 'FWHM_ampli_all_cells_all_freq_cats', directories.figure_dir, darkmode_bool)
