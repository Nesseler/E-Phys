# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 15:30:01 2024

@author: nesseler
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sbn
import numpy as np

from parameters.directories_win import cell_descrip_dir, figure_dir,table_file

# custom functions
from functions.functions_plotting import get_colors, get_figure_size, save_figures, set_font_sizes


fstAP_parameter_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccth1AP-fst_AP_parameters.xlsx'), index_col = 'cell_ID')

fstAP_i_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccth1AP-fst_AP_i.xlsx'), index_col = 'cell_ID')

cell_IDs = fstAP_parameter_df.index.to_list()

todrop = ['v_peaks', 't_peaks', 't_threshold', 'idx_threshold', 'v_AHP', 't_AHP', 'idx_AHP', 'v_AHP_amplitude', 't_to_AHP', 'v_HM', 't1_HM', 't2_HM', 'SR_ms']

fstAP_df = pd.concat([fstAP_parameter_df, fstAP_i_df], axis = 1)


fstAP_df.drop(columns = todrop, inplace = True)





# %% plot

set_font_sizes()

darkmode_bool = True

colors_dict, _ = get_colors(darkmode_bool)

fig_fstAP, axs_fstAP = plt.subplots(nrows = 2,
                                    ncols = 4,
                                    figsize = get_figure_size(),
                                    layout = 'constrained')

axs_fstAP = axs_fstAP.flatten()
fig_fstAP.delaxes(axs_fstAP[-1])

for idx, parameter in enumerate(fstAP_df.columns):

    violin = sbn.violinplot(data = fstAP_df,
                            y = parameter, 
                            inner = 'quart',
                            linewidth = 1,
                            ax = axs_fstAP[idx])

    for l in violin.lines:
        l.set_color(colors_dict['primecolor'])

    for violin in violin.collections:
        violin.set_edgecolor(colors_dict['primecolor'])
        violin.set_facecolor('None')

    sbn.swarmplot(data = fstAP_df, 
                  y = parameter, 
                  ax = axs_fstAP[idx],
                  size = 7,
                  color = colors_dict['primecolor'])


[ax.grid(False) for ax in axs_fstAP]

# set axes 

# v_threshold
axs_fstAP[0].set_ylim([-85, -40])
axs_fstAP[0].set_ylabel('Voltage at threshold [mV]')
axs_fstAP[0].set_yticks(np.arange(-85, -40+1, 10))
axs_fstAP[0].set_yticks(np.arange(-85, -40+1, 5), minor = True)

# amplitude
axs_fstAP[1].set_ylim([40, 140])
axs_fstAP[1].set_ylabel('Spike amplitude [mV]')
axs_fstAP[1].set_yticks(np.arange(40, 140+1, 20))
axs_fstAP[1].set_yticks(np.arange(40, 140+1, 5), minor = True)

# t_topeak
axs_fstAP[2].set_ylim([0.5, 1.5])
axs_fstAP[2].set_ylabel('Time to peak [ms]')
axs_fstAP[2].set_yticks(np.arange(0.5, 1.5 + 0.1, 0.5))
axs_fstAP[2].set_yticks(np.arange(0.5, 1.5 + 0.1, 0.25), minor = True)

# t_rise
axs_fstAP[3].set_ylim([0, 0.75])
axs_fstAP[3].set_ylabel('Rise time [ms]')
axs_fstAP[3].set_yticks(np.arange(0, 0.75 + 0.1, 0.25))
axs_fstAP[3].set_yticks(np.arange(0., 0.75 + 0.05, 0.05), minor = True)

# FWHM
axs_fstAP[4].set_ylim([0.5, 2.0])
axs_fstAP[4].set_ylabel('FWHM [ms]')
axs_fstAP[4].set_yticks(np.arange(0.5, 2. + 0.1, 0.25))
axs_fstAP[4].set_yticks(np.arange(0.5, 2. + 0.05, 0.05), minor = True)

# i_th_abs
axs_fstAP[5].set_ylim([0, 400])
axs_fstAP[5].set_ylabel('Abs. input current [pA]')
axs_fstAP[5].set_yticks(np.arange(0, 400 + 1, 100))
axs_fstAP[5].set_yticks(np.arange(0, 400 + 1, 10), minor = True)

# i_th_rel
axs_fstAP[6].set_ylim([0, 400])
axs_fstAP[6].set_ylabel('Rel. input current [pA]')
axs_fstAP[6].set_yticks(np.arange(0, 400 + 1, 100))
axs_fstAP[6].set_yticks(np.arange(0, 400 + 1, 10), minor = True)

# despine
[ax.spines[spine].set_visible(False) for ax in axs_fstAP for spine in ['bottom', 'top', 'right']]


save_figures(fig_fstAP, 'cc_th1AP_cat', figure_dir, darkmode_bool)




# %% import meta data sheet

MetaData = pd.read_excel(table_file,
                      sheet_name="MetaData",
                      index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

plt_df = pd.concat([fstAP_df, MetaData['Region']], axis = 1)


# %% plot with color coded regions

set_font_sizes()

darkmode_bool = True

colors_dict, regions_c = get_colors(darkmode_bool)


fig_fstAP, axs_fstAP = plt.subplots(nrows = 2,
                                    ncols = 4,
                                    figsize = get_figure_size(),
                                    layout = 'constrained')

axs_fstAP = axs_fstAP.flatten()
fig_fstAP.delaxes(axs_fstAP[-1])

for idx, parameter in enumerate(fstAP_df.columns):

    # identified_regions_only_df = plt_df.query('Region != "BAOT / MeA"')

    violin = sbn.violinplot(data = plt_df,
                            y = parameter,
                            inner = 'quart',
                            linewidth = 1,
                            ax = axs_fstAP[idx])

    for l in violin.lines:
        l.set_color(colors_dict['primecolor'])

    for violin in violin.collections:
        violin.set_edgecolor(colors_dict['primecolor'])
        violin.set_facecolor('None')

    swarm = sbn.swarmplot(data = plt_df,
                          y = parameter, 
                          ax = axs_fstAP[idx],
                          size = 7,
                          color = colors_dict['primecolor'],
                          hue = 'Region',
                          palette = regions_c)

    axs_fstAP[idx].legend().set_visible(False)


    
[ax.grid(False) for ax in axs_fstAP]

# set axes 

# v_threshold
axs_fstAP[0].set_ylim([-85, -40])
axs_fstAP[0].set_ylabel('Voltage at threshold [mV]')
axs_fstAP[0].set_yticks(np.arange(-85, -40+1, 10))
axs_fstAP[0].set_yticks(np.arange(-85, -40+1, 5), minor = True)

# amplitude
axs_fstAP[1].set_ylim([40, 140])
axs_fstAP[1].set_ylabel('Spike amplitude [mV]')
axs_fstAP[1].set_yticks(np.arange(40, 140+1, 20))
axs_fstAP[1].set_yticks(np.arange(40, 140+1, 5), minor = True)

# t_topeak
axs_fstAP[2].set_ylim([0.5, 1.5])
axs_fstAP[2].set_ylabel('Time to peak [ms]')
axs_fstAP[2].set_yticks(np.arange(0.5, 1.5 + 0.1, 0.5))
axs_fstAP[2].set_yticks(np.arange(0.5, 1.5 + 0.1, 0.25), minor = True)

# t_rise
axs_fstAP[3].set_ylim([0, 0.75])
axs_fstAP[3].set_ylabel('Rise time [ms]')
axs_fstAP[3].set_yticks(np.arange(0, 0.75 + 0.1, 0.25))
axs_fstAP[3].set_yticks(np.arange(0., 0.75 + 0.05, 0.05), minor = True)

# FWHM
axs_fstAP[4].set_ylim([0.5, 2.0])
axs_fstAP[4].set_ylabel('FWHM [ms]')
axs_fstAP[4].set_yticks(np.arange(0.5, 2. + 0.1, 0.25))
axs_fstAP[4].set_yticks(np.arange(0.5, 2. + 0.05, 0.05), minor = True)

# i_th_abs
axs_fstAP[5].set_ylim([0, 400])
axs_fstAP[5].set_ylabel('Abs. input current [pA]')
axs_fstAP[5].set_yticks(np.arange(0, 400 + 1, 100))
axs_fstAP[5].set_yticks(np.arange(0, 400 + 1, 10), minor = True)

# i_th_rel
axs_fstAP[6].set_ylim([0, 400])
axs_fstAP[6].set_ylabel('Rel. input current [pA]')
axs_fstAP[6].set_yticks(np.arange(0, 400 + 1, 100))
axs_fstAP[6].set_yticks(np.arange(0, 400 + 1, 10), minor = True)

# despine
[ax.spines[spine].set_visible(False) for ax in axs_fstAP for spine in ['bottom', 'top', 'right']]


save_figures(fig_fstAP, 'cc_th1AP_cat-cc_region', figure_dir, darkmode_bool)



# %% plot with color coded regions and separated violins

set_font_sizes()

darkmode_bool = True

colors_dict, regions_c = get_colors(darkmode_bool)


fig_fstAP, axs_fstAP = plt.subplots(nrows = 2,
                                    ncols = 4,
                                    figsize = get_figure_size(),
                                    layout = 'constrained')

axs_fstAP = axs_fstAP.flatten()
fig_fstAP.delaxes(axs_fstAP[-1])

for idx, parameter in enumerate(fstAP_df.columns):

    # identified_regions_only_df = plt_df.query('Region != "BAOT / MeA"')

    violin = sbn.violinplot(data = plt_df,
                            x = 'Region',
                            y = parameter,
                            inner = 'quart',
                            linewidth = 1,
                            ax = axs_fstAP[idx],
                            size = 0.9)

    for l in violin.lines:
        l.set_color(colors_dict['primecolor'])

    for violin in violin.collections:
        violin.set_edgecolor(colors_dict['primecolor'])
        violin.set_facecolor('None')

    swarm = sbn.swarmplot(data = plt_df,
                          x = 'Region',
                          y = parameter, 
                          ax = axs_fstAP[idx],
                          size = 5,
                          color = colors_dict['primecolor'],
                          hue = 'Region',
                          palette = regions_c)

    axs_fstAP[idx].legend().set_visible(False)


    
[ax.grid(False) for ax in axs_fstAP]
[ax.set_xlabel('') for ax in axs_fstAP]
[ax.set_xticklabels(['BAOT/\nMeA', 'BAOT', 'MeA'], rotation = 45) for ax in axs_fstAP]

# set axes 

# v_threshold
axs_fstAP[0].set_ylim([-85, -40])
axs_fstAP[0].set_ylabel('Voltage at threshold [mV]')
axs_fstAP[0].set_yticks(np.arange(-85, -40+1, 10))
axs_fstAP[0].set_yticks(np.arange(-85, -40+1, 5), minor = True)

# amplitude
axs_fstAP[1].set_ylim([40, 140])
axs_fstAP[1].set_ylabel('Spike amplitude [mV]')
axs_fstAP[1].set_yticks(np.arange(40, 140+1, 20))
axs_fstAP[1].set_yticks(np.arange(40, 140+1, 5), minor = True)

# t_topeak
axs_fstAP[2].set_ylim([0.5, 1.5])
axs_fstAP[2].set_ylabel('Time to peak [ms]')
axs_fstAP[2].set_yticks(np.arange(0.5, 1.5 + 0.1, 0.5))
axs_fstAP[2].set_yticks(np.arange(0.5, 1.5 + 0.1, 0.25), minor = True)

# t_rise
axs_fstAP[3].set_ylim([0, 0.75])
axs_fstAP[3].set_ylabel('Rise time [ms]')
axs_fstAP[3].set_yticks(np.arange(0, 0.75 + 0.1, 0.25))
axs_fstAP[3].set_yticks(np.arange(0., 0.75 + 0.05, 0.05), minor = True)

# FWHM
axs_fstAP[4].set_ylim([0.5, 2.0])
axs_fstAP[4].set_ylabel('FWHM [ms]')
axs_fstAP[4].set_yticks(np.arange(0.5, 2. + 0.1, 0.25))
axs_fstAP[4].set_yticks(np.arange(0.5, 2. + 0.05, 0.05), minor = True)

# i_th_abs
axs_fstAP[5].set_ylim([0, 400])
axs_fstAP[5].set_ylabel('Abs. input current [pA]')
axs_fstAP[5].set_yticks(np.arange(0, 400 + 1, 100))
axs_fstAP[5].set_yticks(np.arange(0, 400 + 1, 10), minor = True)

# i_th_rel
axs_fstAP[6].set_ylim([0, 400])
axs_fstAP[6].set_ylabel('Rel. input current [pA]')
axs_fstAP[6].set_yticks(np.arange(0, 400 + 1, 100))
axs_fstAP[6].set_yticks(np.arange(0, 400 + 1, 10), minor = True)

# despine
[ax.spines[spine].set_visible(False) for ax in axs_fstAP for spine in ['bottom', 'top', 'right']]


save_figures(fig_fstAP, 'cc_th1AP_cat-cc_region+separate_violins', figure_dir, darkmode_bool)