# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 17:25:56 2024

@author: nesseler
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn

from os.path import join

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir, figure_dir, quant_data_dir, table_file
from parameters.PGFs import cc_APs_parameters, cc_APs_t_stims_df
from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size


frequencies = cc_APs_parameters.keys()
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

resul_freq_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccAPs-resul_freq.xlsx'), index_col = 'frequencies')
cell_IDs = list(resul_freq_df.columns)

darkmode_bool = True

# import meta data sheet
MetaData = pd.read_excel(table_file,
                      sheet_name="MetaData",
                      index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

# %%


nAPs = 3

perc_of_fst_df = pd.DataFrame(index=frequencies, columns=cell_IDs)
perc_gain_loss_df = pd.DataFrame(index=frequencies, columns=cell_IDs)
abs_gain_loss_df = pd.DataFrame(index=frequencies, columns=cell_IDs)


# cell_IDs = ['E-087']


for cell_ID in cell_IDs:
    
    cell_APs_path = join(quant_data_dir, 'APs', cell_ID)
    
    # ISI_df = pd.DataFrame(columns = frequencies)
    
    for frequency in frequencies:
    
        # load dataframe with AP parameters for cell
        AP_params = pd.read_excel(join(cell_APs_path, f'{cell_ID}_{frequency}.xlsx'))
    
        # # filter to dataframe with APs
        # AP_params_onlyAPs = AP_params.query('v_amplitude.notnull()')
        
        # get AP times
        AP_times = cc_APs_t_stims_df[frequency] + AP_params['t_peaks']
        
        # calc ISIs
        ISIs = np.diff(AP_times.dropna())
        
        # get first ISI
        ISI_1stAP = ISIs[0]
        
        # get mean of last n APs
        ISIs_lstAPs = ISIs[-(1 + nAPs):-1]
        mean_ISI_lstAPs = np.mean(ISIs_lstAPs)
        
        # calc percentage of first, loss/gain (absolute and relative)
        ## percentage of first
        perc_of_fst = mean_ISI_lstAPs / ISI_1stAP
        perc_of_fst_df.at[frequency, cell_ID] = perc_of_fst

        ## absolute gain/loss
        abs_gain_loss = mean_ISI_lstAPs - ISI_1stAP
        abs_gain_loss_df.at[frequency, cell_ID] = abs_gain_loss
        
        ## percentage gain/loss
        perc_gain_loss = abs_gain_loss / ISI_1stAP
        perc_gain_loss_df.at[frequency, cell_ID] = perc_gain_loss
    
        
        
    print(f'Finished: {cell_ID}')




# %% set plot settings

# define plotting dataframe and get mean and std
plt_df = perc_of_fst_df
mean = plt_df.mean(axis = 1)
std = plt_df.std(axis = 1)

# define colors
colors_dict, region_c = get_colors(darkmode_bool)


# %% plot percentage of first spike

fig_1st_v_lst, axs_1st_v_lst = plt.subplots(nrows = 1,
                                            ncols = 1,
                                            layout = 'constrained',
                                            dpi = 600)

for cell_ID in cell_IDs:
    axs_1st_v_lst.plot(freqs_int, plt_df[cell_ID],
                       color = 'gray',
                       alpha = 0.5)


axs_1st_v_lst.errorbar(x = freqs_int,
                       y = mean,
                       yerr = std,
                       color = colors_dict['primecolor'],
                       ecolor = colors_dict['primecolor'],
                       linestyle = '--',
                       capsize = 3,
                       marker = '_',
                       markersize = 10)


axs_1st_v_lst.set_xlim([-1, 77])
axs_1st_v_lst.set_xlabel('Stimulation frequency [Hz]')
axs_1st_v_lst.set_xticks(ticks = freqs_int)

# x axis at zero
axs_1st_v_lst.axhline(y = 1.0, xmin = 0, xmax = 1,
                      lw = 1,
                      color = colors_dict['primecolor'])


# despine
[axs_1st_v_lst.spines[spine].set_visible(False) for spine in ['top', 'right']]

# xlimit axis spines to their limits
axs_1st_v_lst.spines['bottom'].set_bounds([1, 75])

axs_1st_v_lst.grid(False)





# %% same plot with color coded regions


fig_ccmeta, axs_ccmeta = plt.subplots(nrows = 1,
                                      ncols = 1,
                                      layout = 'constrained',
                                      dpi = 600)

for cell_ID in cell_IDs:
    
    cell_region = MetaData.at[cell_ID, 'Region']
    
    axs_ccmeta.plot(freqs_int, plt_df[cell_ID],
                    color = region_c[cell_region],
                    lw = 1,
                    marker = 'x')


axs_ccmeta.set_xlim([-1, 77])
axs_ccmeta.set_xlabel('Stimulation frequency [Hz]')
axs_ccmeta.set_xticks(ticks = freqs_int)

# x axis at zero
axs_ccmeta.axhline(y = 1.0, xmin = 0, xmax = 1,
                      lw = 1,
                      color = colors_dict['primecolor'])


# despine
[axs_ccmeta.spines[spine].set_visible(False) for spine in ['top', 'right']]

# xlimit axis spines to their limits
axs_ccmeta.spines['bottom'].set_bounds([1, 75])
axs_ccmeta.grid(False)


# %% same plot with color coded regions + separate planes

regions = ['BAOT/MeA', 'MeA', 'BAOT']

fig_regions, axs_regions = plt.subplots(nrows = 3,
                                        ncols = 1,
                                        layout = 'constrained',
                                        dpi = 600,
                                        figsize = get_figure_size(width = 161.835))


set_font_sizes()

# separate panels for regions

for idx_region, region in enumerate(regions):
    
    # get cell_IDs of dataframe per region
    cell_IDs_region = MetaData[MetaData['Region'] == region].index.to_list()
    
    # get number of cells in region
    n_cells_region = len(cell_IDs_region)
    
    ax = axs_regions[idx_region]
    
    for idx_cell, cell_ID in enumerate(cell_IDs_region):

        cell_region = MetaData.at[cell_ID, 'Region']
        
        ax.plot(freqs_int, plt_df[cell_ID],
                            color = region_c[cell_region],
                            lw = 1,
                            marker = 'x')


    ax.set_ylim([-1, 7])

    ax.set_xlim([-1, 77])

    ax.set_xticks(ticks = freqs_int, labels = [])
    # ax.tick_params('both', width = 5, length = 5, color = 'k', direction = 'in')

    # x axis at zero
    ax.axhline(y = 1.0, xmin = 0, xmax = 1,
               lw = 1,
               color = colors_dict['primecolor'],
               linestyle = '--')

    # xlimit axis spines to their limits
    ax.spines['bottom'].set_bounds([1, 75])

axs_regions[-1].set_xlabel('Stimulation frequency [Hz]')
axs_regions[-1].set_xticks(ticks = freqs_int, labels = freqs_int)
fig_regions.supylabel(f'first ISI / mean of last {nAPs} ISIs')


# despine
[ax.spines[spine].set_visible(False) for spine in ['top', 'right'] for ax in axs_regions]

[ax.grid(False) for ax in axs_regions]

save_figures(fig_regions, 'ccAPs_ISI_adaptation-sep_regions', figure_dir, darkmode_bool)


