# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 16:59:43 2024

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

def import_AP_measurement_one_cell(cell_data_path, frequency, parameter, cell_ID):

    # get number of APs for each stim
    measurement_df = pd.DataFrame()
    
    # create file path
    file_path = os.path.join(cell_data_path, cell_ID, f'{cell_ID}_{frequency}.xlsx')
    
    # read excel file
    AP_all_params = pd.read_excel(file_path, index_col = 'idx_step')
    
    # write to dictionary
    measurement_df[frequency] = AP_all_params[parameter]

    return measurement_df


# %%

parameter = 'FWHM'

mean_df = pd.DataFrame()

# loop through cells and then frequencies to get mean value of desired parameter

for cell_ID in cell_IDs: 
    
    mean_p_freqs = []
    
    for frequency in frequencies:
     
        freq_int = int(frequency.replace('Hz', ''))
        freq_idx = frequencies.index(frequency)
        
        cur_df = import_AP_measurement_one_cell(cell_data_path = quant_dir,
                                                frequency = frequency,
                                                parameter = parameter,
                                                cell_ID = cell_ID)
        
        mean_p_freqs.append(cur_df[frequency][:].mean())
    
    # create dataframe with index to add to the mean_df
    mean_p_freqs_df = pd.DataFrame({cell_ID : mean_p_freqs}, 
                                    index = freqs_int)
    
    # mean_p_freqs_df = pd.DataFrame({cell_ID : mean_p_freqs,
    #                                 'freqs_str' : frequencies}, 
    #                                index = freqs_int)    
    
    mean_df[cell_ID] = mean_p_freqs_df
        
mean_df = mean_df.transpose()

mean_df_melt = pd.melt(mean_df, var_name='frequency', value_name=parameter, ignore_index= False)

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

# %% 

mean_df = mean_df.transpose()
mean_df['freq_str'] = frequencies

mean_df_melt['colFromIndex'] = mean_df_melt.index
mean_df_melt['freq_str'] = mean_df_melt['frequency'].apply(lambda x: str(x) + 'Hz')
mean_df_melt = mean_df_melt.sort_values(by = ['colFromIndex', 'frequency'], ascending = [True, True])



# %% plotting

darkmode_bool = True

colors_dict = get_colors(darkmode_bool)

fig_APcats, axs_APcats = plt.subplots(nrows = 1, 
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


violins1 = sbn.violinplot(data = mean_df_melt, 
                          x = 'freq_str', 
                          y = parameter, 
                          inner = 'quart', 
                          ax = axs_APcats, 
                          linewidth = 1,
                          palette = color_pal, 
                          width = 0.9)

swarms = sbn.swarmplot(data = mean_df_melt, 
                        x = 'freq_str', 
                        y = parameter, 
                        ax = axs_APcats,
                        size = 6,
                        color = colors_dict['primecolor'])


# get positions of all points in swarmplots to plot the connecting lines
positions_df = pd.DataFrame()

for idx, frequency in enumerate(frequencies):
    positions = np.array(axs_APcats.collections[6 + idx].get_offsets())
    
    cur_df = pd.DataFrame({'x' : positions[:, 0],
                           'y' : positions[:, 1],
                           'frequency' : [freqs_int[idx]] * len(positions),
                           'freq_str' : [frequency] * len(positions),
                           'cell_ID' : cell_IDs
                            })

    positions_df = pd.concat([positions_df, cur_df])
    
     

for idx, cell_ID in enumerate(cell_IDs):
    lines = sbn.lineplot(data = positions_df[positions_df.cell_ID == cell_ID], 
                            x = 'x', 
                            y = 'y',
                            ax = axs_APcats,
                            estimator=None,
                            color = colors_dict['primecolor'],
                            alpha = 0.3)

for l in violins1.lines:
    l.set_color(colors_dict['primecolor'])

[violin.set_edgecolor(colors_dict['primecolor']) for violin in violins1.collections]
        


axs_APcats.grid(False)

plt.show()



save_figures(fig_APcats, f'mean_{parameter}_cell_all_freq_cats', directories.figure_dir, darkmode_bool)























