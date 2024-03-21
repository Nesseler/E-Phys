# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 14:55:42 2024

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
from parameters.PGFs import cc_APs_parameters
from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size


frequencies = cc_APs_parameters.keys()
freqs_int = [int(f_str.replace('Hz', '')) for f_str in frequencies]

resul_freq_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccAPs-resul_freq.xlsx'), index_col = 'frequencies')
cell_IDs = list(resul_freq_df.columns)

darkmode_bool = True

# %%



AP_parameter = 'FWHM'

nAPs = 3

perc_of_fst_df = pd.DataFrame(index=frequencies, columns=cell_IDs)
perc_gain_loss_df = pd.DataFrame(index=frequencies, columns=cell_IDs)
abs_gain_loss_df = pd.DataFrame(index=frequencies, columns=cell_IDs)

# cell_IDs = ['E-087']


for cell_ID in cell_IDs:
    
    cell_APs_path = join(quant_data_dir, 'APs', cell_ID)
    
    for frequency in frequencies:
    
        # load dataframe with AP parameters for cell
        AP_params = pd.read_excel(join(cell_APs_path, f'{cell_ID}_{frequency}.xlsx'))
    
        # filter to dataframe with APs
        AP_params_onlyAPs = AP_params.query('v_amplitude.notnull()')
        
        # get first AP
        AP_params_1stAP = AP_params_onlyAPs.iloc[0]
        parameter_1stAP = AP_params_1stAP[AP_parameter]
        
        # get mean of last n APs
        ## get their indices
        idc_AP_params_lstAPs = AP_params_onlyAPs.iloc[-(1 + nAPs):-1].index.to_list()
        
        ## limit dataframe to only last APs
        AP_params_lstAPs = AP_params.iloc[idc_AP_params_lstAPs]
        
        # get mean of given parameter
        parameter_mean_lstAPs =AP_params_lstAPs[AP_parameter].mean()
        
        # calc percentage of first, loss/gain (absolute and relative)
        ## percentage of first
        perc_of_fst = parameter_mean_lstAPs / parameter_1stAP
        perc_of_fst_df.at[frequency, cell_ID] = perc_of_fst

        ## absolute gain/loss
        abs_gain_loss = parameter_mean_lstAPs - parameter_1stAP
        abs_gain_loss_df.at[frequency, cell_ID] = abs_gain_loss
        
        ## percentage gain/loss
        perc_gain_loss = abs_gain_loss / parameter_1stAP
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
                                            dpi = 600,
                                            figsize = get_figure_size(width = 165.5))

set_font_sizes()

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
                      color = colors_dict['primecolor'],
                      linestyle = '--')


# despine
[axs_1st_v_lst.spines[spine].set_visible(False) for spine in ['top', 'right']]

# xlimit axis spines to their limits
axs_1st_v_lst.spines['bottom'].set_bounds([1, 75])


# y axis settings per parameter
axs_1st_v_lst.set_ylabel(f'Percentage of first spike {AP_parameter}')

if AP_parameter == 'v_amplitude':
    axs_1st_v_lst.set_ylim([0.45, 1.15])
    
elif AP_parameter == 'FWHM':
    axs_1st_v_lst.set_ylim([0.5, 3])
    
    
axs_1st_v_lst.grid(False)

save_figures(fig_1st_v_lst, f'ccAPs-{AP_parameter}-perc_of_fst', figure_dir, darkmode_bool)



# %% import meta data sheet

MetaData = pd.read_excel(table_file,
                      sheet_name="MetaData",
                      index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]


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



# axs_ccmeta.errorbar(x = freqs_int,
#                        y = mean,
#                        yerr = std,
#                        color = colors_dict['primecolor'],
#                        ecolor = colors_dict['primecolor'],
#                        linestyle = '--',
#                        capsize = 3,
#                        marker = '_',
#                        markersize = 10)


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

# y axis settings per parameter
axs_ccmeta.set_ylabel(f'Percentage of first spike {AP_parameter}')

if AP_parameter == 'v_amplitude':
    axs_ccmeta.set_ylim([0.45, 1.15])
    
elif AP_parameter == 'FWHM':
    axs_ccmeta.set_ylim([0.5, 3])
    
    


save_figures(fig_ccmeta, f'ccAPs-{AP_parameter}-perc_of_fst-cc_region', figure_dir, darkmode_bool)


# %% separate subplots for the regions

fig_ccmeta_sep, axs_ccmeta_sep = plt.subplots(nrows = 3,
                                              ncols = 1,
                                              layout = 'constrained',
                                              dpi = 600,
                                              sharey = True,
                                              sharex = True,
                                              figsize = get_figure_size(width = 165.5))


regions = ['BAOT/MeA', 'MeA', 'BAOT']

for region_idx, region in enumerate(regions):
    
    for cell_ID in cell_IDs:
        
        cell_region = MetaData.at[cell_ID, 'Region']
        
        if cell_region == region:
            
            axs_ccmeta_sep[region_idx].plot(freqs_int, plt_df[cell_ID],
                                            color = region_c[cell_region],
                                            lw = 1,
                                            marker = 'x')


axs_ccmeta_sep[-1].set_xlim([-1, 77])
axs_ccmeta_sep[-1].set_xlabel('Stimulation frequency [Hz]')
axs_ccmeta_sep[-1].set_xticks(ticks = freqs_int)

# y axis settings per parameter
fig_ccmeta_sep.supylabel(f'Percentage of first spike {AP_parameter}')

for ax in axs_ccmeta_sep:
    # x axis at zero
    ax.axhline(y = 1.0, xmin = 0, xmax = 1,
                          lw = 1,
                          color = colors_dict['primecolor'])
    
    
    # despine
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # xlimit axis spines to their limits
    ax.spines['bottom'].set_bounds([1, 75])
    
    

    
    if AP_parameter == 'v_amplitude':
        ax.set_ylim([0.45, 1.15])
        
    elif AP_parameter == 'FWHM':
        ax.set_ylim([0.5, 3])


[ax.grid(False) for ax in axs_ccmeta_sep]


save_figures(fig_ccmeta_sep, f'ccAPs-{AP_parameter}-perc_of_fst-cc_region+sep_panels', figure_dir, darkmode_bool)




# %% plot again resulting frequency not stimulated

resul_freq_df = pd.read_excel(os.path.join(cell_descrip_dir, 'ccAPs-resul_freq.xlsx'), index_col = 'frequencies')



fig_rfreq, axs_rfreq = plt.subplots(nrows = 1,
                                            ncols = 1,
                                            layout = 'constrained',
                                            dpi = 600)

for cell_ID in cell_IDs:
    
    cell_region = MetaData.at[cell_ID, 'Region']
    
    # sort x and y value after resulting firing frequency
    x = resul_freq_df[cell_ID]
    x = x.rename('resul_freq')
    y = plt_df[cell_ID]
    y = y.rename('prec_of_first')
    
    # concate x and y
    xy = pd.concat([x, y], axis = 1)
    
    # sort for resulting frequency
    xy = xy.sort_values('resul_freq')
    
    
    axs_rfreq.plot(xy['resul_freq'], xy['prec_of_first'],
                       color = region_c[cell_region] ,#regions_c[cell_region],
                       alpha = 1)


# axs_ccmeta.errorbar(x = freqs_int,
#                        y = mean,
#                        yerr = std,
#                        color = colors_dict['primecolor'],
#                        ecolor = colors_dict['primecolor'],
#                        linestyle = '--',
#                        capsize = 3,
#                        marker = '_',
#                        markersize = 10)


axs_rfreq.set_xlim([-1, 77])
axs_rfreq.set_xlabel('Resulting inst. firing frequency [Hz]')
axs_rfreq.set_xticks(ticks = np.arange(0, 75, 10))

# x axis at zero
axs_rfreq.axhline(y = 1.0, xmin = 0, xmax = 1,
                      lw = 1,
                      color = colors_dict['primecolor'])


# despine
[axs_rfreq.spines[spine].set_visible(False) for spine in ['top', 'right']]

# xlimit axis spines to their limits
axs_rfreq.spines['bottom'].set_bounds([0, 75])


# y axis settings per parameter
axs_rfreq.set_ylabel(f'Percentage of first spike {AP_parameter}')

if AP_parameter == 'v_amplitude':
    axs_rfreq.set_ylim([0.45, 1.15])
    
elif AP_parameter == 'FWHM':
    axs_rfreq.set_ylim([0.5, 3])
    
    
axs_rfreq.grid(False)

save_figures(fig_rfreq, f'ccAPs-{AP_parameter}-perc_of_fst-rfreq', figure_dir, darkmode_bool)



# %%

region_df = MetaData['Region']

plt_df = perc_of_fst_df

plt_df = plt_df.transpose()



plt_df_melted = pd.melt(plt_df, var_name='frequency', value_name=AP_parameter, ignore_index= False) 

plt_df_melted['Region'] = region_df
plt_df_melted.index.name = 'cell_ID'

# %%

fig_ccmeta, axs_ccmeta = plt.subplots(nrows = 1,
                                            ncols = 1,
                                            layout = 'constrained',
                                            dpi = 600,
                                            figsize = get_figure_size())

violins = sbn.violinplot(data = plt_df_melted,
                         y = plt_df_melted[AP_parameter].to_list(),
                         x = 'frequency',
                         hue = 'Region',
                         inner = 'quart',
                         palette = ['k', 'k', 'k'],
                         hue_order = regions,
                         linewidth = 1.5)


# 6 stimulation frequencies, 3 regions and 3 lines (for quartiles) in violins
region_per_violin = []
[region_per_violin.append(region) for freq in frequencies for region in regions for i in range(3)]

for idx_l, l in enumerate(violins.lines):
    l.set_color(region_c[region_per_violin[idx_l]])

# 6 stimulation frequencies, 3 regions and 3 lines (for quartiles) in violins
region_per_violin = []
[region_per_violin.append(region) for freq in frequencies for region in regions]

for idx_violin, violin in enumerate(violins.collections):
    violin.set_edgecolor(region_c[region_per_violin[idx_violin]])
    


swarms = sbn.swarmplot(plt_df_melted,
                       x = 'frequency',
                       y = AP_parameter,
                       hue = 'Region',
                       palette = region_c,
                       dodge = True,
                       size = 5,
                       hue_order = regions)


for region_idx, region in enumerate(regions):
    
    print(region_idx, region)
    
    print(MetaData[MetaData['Region'] == region].index.to_list())
    
    region_cell_IDs = MetaData[MetaData['Region'] == region].index.to_list()

    # get positions of all points in swarmplots to plot the connecting lines
    positions_df = pd.DataFrame()
    
    for idx, frequency in enumerate(frequencies):
        positions = np.array(axs_ccmeta.collections[(idx*3)+region_idx+18].get_offsets())
        
        cur_df = pd.DataFrame({'x' : positions[:, 0],
                                'y' : positions[:, 1],
                                'frequency' : [freqs_int[idx]] * len(positions),
                                'freq_str' : [frequency] * len(positions),
                                'cell_ID' : region_cell_IDs
                                })
    
        positions_df = pd.concat([positions_df, cur_df])
        


    for idx, cell_ID in enumerate(region_cell_IDs):
        lines = sbn.lineplot(data = positions_df[positions_df.cell_ID == cell_ID], 
                              x = 'x', 
                              y = 'y',
                              ax = axs_ccmeta,
                              estimator=None,
                              color = region_c[region],
                              alpha = 0.3)



# axs_ccmeta.errorbar(x = freqs_int,
#                         y = mean,
#                         yerr = std,
#                         color = colors_dict['primecolor'],
#                         ecolor = colors_dict['primecolor'],
#                         linestyle = '--',
#                         capsize = 3,
#                         marker = '_',
#                         markersize = 10)


# axs_ccmeta.set_xlim([-1, 77])
# axs_ccmeta.set_xlabel('Stimulation frequency [Hz]')
# axs_ccmeta.set_xticks(ticks = freqs_int)

# x axis at one
axs_ccmeta.axhline(y = 1.0, xmin = 0, xmax = 1,
                      lw = 1,
                      color = colors_dict['primecolor'],
                      linestyle = '--')


# # despine
# [axs_ccmeta.spines[spine].set_visible(False) for spine in ['top', 'right']]

# # xlimit axis spines to their limits
# axs_ccmeta.spines['bottom'].set_bounds([1, 75])


# y axis settings per parameter
axs_ccmeta.set_ylabel(f'Percentage of first spike {AP_parameter}')

if AP_parameter == 'v_amplitude':
    axs_ccmeta.set_ylim([0.45, 1.15])
    
elif AP_parameter == 'FWHM':
    axs_ccmeta.set_ylim([0.5, 3])
    
    
axs_ccmeta.grid(False)

save_figures(fig_ccmeta, f'ccAPs-{AP_parameter}-perc_of_fst-violins-regions', figure_dir, darkmode_bool)
