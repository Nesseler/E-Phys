# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 19:21:26 2024

@author: nesseler
"""

import pandas as pd
import seaborn as sbn
from os.path import join
import matplotlib.pyplot as plt
import numpy as np


from parameters.directories_win import cell_descrip_dir, table_file, figure_dir, cell_morph_descrip_dir

from parameters.PGFs import cc_IF_parameters

from functions.functions_plotting import get_colors, get_figure_size, set_font_sizes, save_figures


# %% load dataframes to describe cells

# cc_rest
activity_df = pd.read_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col = 'cell_ID')

# cc_IF
active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')
passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col = 'cell_ID')
fstAP_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_col = 'cell_ID')
IF_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_col = 'i_input')
IF_inst_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_col = 'i_input')
IF_inst_initial_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst_initial.xlsx'), index_col = 'i_input')

# cc_sag
sag_df = pd.read_excel(join(cell_descrip_dir, 'cc_sag-sagdelta.xlsx'), index_col = 'cell_ID')

# cc_th1AP
th1AP_iinput = pd.read_excel(join(cell_descrip_dir, 'ccth1AP-fst_AP_i.xlsx'), index_col = 'cell_ID')
th1AP_parameters = pd.read_excel(join(cell_descrip_dir, 'ccth1AP-fst_AP_parameters.xlsx'), index_col = 'cell_ID')

# cc_APs
result_freqs_df = pd.read_excel(join(cell_descrip_dir, 'ccAPs-resul_freq.xlsx'), index_col = 'frequencies')

# sholl metrics
sholl_metrics_df = pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics.xlsx'), index_col = 'cell_ID')
sholl_metrics_df = sholl_metrics_df.loc[[cell_ID for cell_ID in activity_df.index.tolist() if cell_ID in sholl_metrics_df.index.to_list()], :]

# MetaData
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# add MetaData
MetaData = MetaData.loc[passiv_properties_df.index.to_list(), :]
activity_df = activity_df.loc[passiv_properties_df.index.to_list(), :]
sag_df = sag_df.loc[passiv_properties_df.index.to_list(), :]
# celldescriptors_df = pd.concat([celldescriptors_df, MetaData[['Region', 'Hemisphere', 'Sex', 'Stage', 'reconstructed']]], axis = 1)

fstAP_df['t_peaks'] = fstAP_df['t_peaks'] - 250 #ms (start of step)

# concatenate all dataframes
plt_df = pd.concat([active_properties_df, 
                    passiv_properties_df, 
                    activity_df,
                    fstAP_df,
                    sag_df,
                    MetaData['Region']], axis = 1)

# drop cell-IDs
from hierarchical_clustering.ePhys_hierarchical_parameters import cell_IDs_toDrop

plt_df.drop(index = cell_IDs_toDrop, inplace = True)

# params_toplot = ['r_input', 'tau_mem', 'c_mem',
#                  'rheobase_rel']

# params_toplot['']

params_toplot = {'1':'r_input', 
                  '2' : 'tau_mem', 
                  '3' : 'c_mem',
                 
                  '5' : 'rheobase_rel',
                  '6' : 't_peaks',
                  '7' : 'n_rheobasespikes',
                 
                 
                  # '8' : 'v_threshold',                 
                  # '9' : 'delta_vrest_to_vthres',
                  # '10': 'v_AHP_amplitude', 
                  # '11': 'FWHM',
                 
                  '8' : 'delta_vrest_to_vthres',                 
                  '9' : 'FWHM',
                  '10': 't_toPeak', 
                  '11': 'v_AHP_amplitude',
                 
                  '13': 'sag_delta',
                  '15': 'n_reboundspikes',
                 
                  '17': 'max_freq',
                  '18': 'max_inst_freq',
                  '19': 'max_inst_initial_freq'}


# params_toplot = {'1' : 'sag_delta', 
                 
#                  '5' : 'n_reboundspikes',
#                  '6' : 'reboundspike_t_toPeak',
#                  '7' : 'reboundspike_FWHM',
#                  }


plt_df['delta_vrest_to_vthres'] = plt_df['v_threshold'] - plt_df['v_rest']


# %% initialize plotting

import matplotlib as mtl
from functions.functions_plotting import save_figures, get_colors, get_figure_size, plot_half_violin


# set colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)



# set regions
regions = ['BAOT/MeA', 'MeA', 'BAOT']


# initialize figure
fig_cats, axs_cats = plt.subplots(nrows = 5,
                                  ncols = 4,
                                  figsize = get_figure_size(width = 277.25, height = 315.165),
                                  layout = 'constrained',
                                  dpi = 600,
                                  sharex = True)


# flatten numpy array of axis for easier handling
axs_cats = axs_cats.flatten()

# iterate through parameter
for key, parameter in params_toplot.items():
    
    # set axis
    ax = axs_cats[int(key)]

    swarm = sbn.swarmplot(data = plt_df,
                          x = 'Region',
                          y = parameter, 
                          ax = ax,
                          hue = 'Region', 
                          palette = region_colors,
                          size = 3,
                          color = colors_dict['primecolor'],
                          order = regions,
                          zorder = 1)
    
    # remove legend
    ax.legend().set_visible(False)
    
    # set positions of violins
    v_positions = [0, 1, 2]
    
    for r_idx, region in enumerate(regions):
        
        # get data
        violin_data = plt_df[plt_df['Region'] == region]
        violin_data = violin_data[parameter].dropna().to_list()
        
        # get mean, median, and std
        mean = np.mean(violin_data)
        median = np.median(violin_data)
        std = np.std(violin_data)
        n = len(violin_data)
    
        # set offset
        offset = 0.20
    
        # set bandwidth
        bandwidth = n**(-1./(1+2))
        # print(n**(-1./(1+4)), bandwidth)
        
        # set resolution of violin in condition
        if parameter in ['t_toPeak', 'FWHM', 'n_reboundspikes', 'n_rheobasespikes']:
            v_resolution = 0.001
        elif parameter == 'r_input':
            v_resolution = 1
        else:
            v_resolution = 0.1
        
        # set violin cutoff
        if parameter in ['v_AHP_amplitude', 'delta_vrest_to_vthres']:
            v_abs_cutoff = [np.nan, np.nan]
        else:
            v_abs_cutoff = [0, np.nan]
    
        # plot half violin
        plot_half_violin(data = violin_data, 
                         ax = ax,
                         v_resolution = v_resolution,
                         v_kde_cutoff = 0.025,
                         v_abs_cutoff = v_abs_cutoff,
                         v_bandwidth = bandwidth,
                         v_position = v_positions[r_idx],
                         v_offset = -offset,
                         v_width = 0.35,
                         v_baseline = True,
                         v_color = region_colors[region],
                         v_zorder = 0)
        
        # set x position of errorbar
        e_position = v_positions[r_idx] + offset
        
        # plot errorbar
        ax.errorbar(x = e_position,
                    y = mean,
                    yerr = std,
                    fmt='_', 
                    markersize = 7,
                    markerfacecolor = 'none',
                    capsize = 3,
                    color=colors_dict['primecolor'],
                    linewidth = 1,
                    label = '_nolegend_',
                    zorder = 2)
        
        # plot median
        ax.scatter(x = e_position,
                    y = median,
                    marker='D', 
                    s = 5,
                    color=colors_dict['primecolor'],
                    linewidth = 1,
                    label = '_nolegend_',
                    zorder = 3)
    
        
    if parameter == 'r_input':
        ax.set_ylabel('Input resistance\n' + r'[M$\Omega$]') 
        ax.set_ylim([-20, 1820])
        ax.spines['left'].set_bounds([0, 1800])
        ax.set_yticks(np.arange(0, 1800+1, 500))
        ax.set_yticks(np.arange(0, 1800+1, 100), minor = True)  
        
    elif parameter == 'tau_mem':
        ax.set_ylabel('Membrane time constant\n[ms]')
        ax.set_ylim([-2, 42])
        ax.spines['left'].set_bounds([0, 40])
        ax.set_yticks(np.arange(0, 40+1, 10))
        ax.set_yticks(np.arange(0, 40+1, 5), minor = True)
        
    elif parameter == 'c_mem':
        ax.set_ylabel('Membrane capacitance\n[pF]')
        ax.set_ylim([-2, 102])
        ax.spines['left'].set_bounds([0, 100])
        ax.set_yticks(np.arange(0, 100+1, 20))
        ax.set_yticks(np.arange(0, 100+1, 5), minor = True)

    elif parameter == 'delta_vrest_to_vthres':
        ax.set_ylabel('Delta v_rest to rheobase\nspike threshold [mV]')
        ax.set_ylim([-5-1, 45+1])
        ax.spines['left'].set_bounds([-5, 45])
        ax.set_yticks(np.arange(0, 45+1, 10))
        ax.set_yticks(np.arange(-5, 45+1, 5), minor = True)
        
    elif parameter == 'v_thres_rheobase_spike':
        ax.set_ylabel('Voltage at threshold\nof rheobase spike [mV]')
        ax.set_ylim([-87, -28])
        ax.spines['left'].set_bounds([-85, -30])
        ax.set_yticks(np.arange(-80, -30+1, 10))
        ax.set_yticks(np.arange(-85, -30+1, 5), minor = True)

    elif parameter == 'rheobase_rel':
        ax.set_ylabel('Relative rheobase\n[pA]')
        ax.set_ylim([0-2, 210+2])
        ax.spines['left'].set_bounds([0, 210])
        ax.set_yticks(np.arange(0, 210+1, 50))
        ax.set_yticks(np.arange(0, 210+1, 10), minor = True)

    elif parameter == 'n_rheobasespikes':
        ax.set_ylabel('Number of\nrheobase spikes [#]')
        ax.set_ylim([0-0.1, 13+0.1])
        ax.spines['left'].set_bounds([0, 13])
        ax.set_yticks(np.arange(0,13+1, 5))
        ax.set_yticks(np.arange(0,13+1, 1), minor = True)

    elif parameter == 't_peaks':
        ax.set_ylabel('Time to rheobase spike\n[ms]')
        ax.set_ylim([0-10, 800+10])
        ax.spines['left'].set_bounds([0, 800])
        ax.set_yticks(np.arange(0,800+1, 200))
        ax.set_yticks(np.arange(0,800+1, 50), minor = True)

    elif parameter == 'v_threshold':
        ax.set_ylabel('Voltage at threshold\nof rheobase spike [mV]')
        ax.set_ylim([-80-2, -35+2])
        ax.spines['left'].set_bounds([-80, -35])
        ax.set_yticks(np.arange(-80, -35+1, 10))
        ax.set_yticks(np.arange(-80, -35+1, 5), minor = True)
        
    elif parameter == 'v_AHP_amplitude':
        ax.set_ylabel('Rheobase spike\nthreshold to AHP [mV]')
        ax.set_ylim([-20-4, 20+4])
        ax.spines['left'].set_bounds([-20, 20])
        ax.set_yticks(np.arange(-20, 20+4, 10))
        ax.set_yticks(np.arange(-20, 20+4, 2.5), minor = True)  

    elif parameter == 'FWHM':
        ax.set_ylabel('Rheobase spike\nFWHM [ms]')
        ax.set_ylim([0.5-0.1, 2+0.1])
        ax.spines['left'].set_bounds([0.5, 2])
        ax.set_yticks(np.arange(0.5, 2+.1, 0.5))
        ax.set_yticks(np.arange(0.5, 2+.1, 0.25), minor = True)
        
    elif parameter == 't_toPeak':
        ax.set_ylabel('Rheobase spike\ntime to peak [ms]')
        ax.set_ylim([0.5-0.1, 2+0.1])
        ax.spines['left'].set_bounds([0.5, 2])
        ax.set_yticks(np.arange(0.5, 2+.1, 0.5))
        ax.set_yticks(np.arange(0.5, 2+.1, 0.25), minor = True)
 
    # # # HCN # # #
    elif parameter == 'sag_delta':
        ax.set_ylabel('Sag potential\ndelta [mV]')
        ax.set_ylim([0-1, 17.5 +1])
        ax.spines['left'].set_bounds([0, 17.5])
        ax.set_yticks(np.arange(0, 17.5+1, 5))
        ax.set_yticks(np.arange(0, 17.5+1, 2.5), minor = True)     
        
    elif parameter == 'n_reboundspikes':
        ax.set_ylabel('Number of\nrebound spikes [#]')
        ax.set_ylim([0-.5, 4+.5])
        ax.spines['left'].set_bounds([0, 4])
        ax.set_yticks(np.arange(0, 4+.5, 2))

     
    elif parameter == 'reboundspike_t_peak':
        ax.set_ylabel('Time to peak [ms]')
        ax.set_ylim([0-5, 125+5])
        ax.spines['left'].set_bounds([0, 125])
        ax.set_yticks(np.arange(0, 125+1, 25))
        ax.set_yticks(np.arange(0, 125+1, 5), minor = True)
        
    elif parameter == 'reboundspike_v_threshold':
        ax.set_ylabel('Membrane voltage\nat reboundspike threshold [mV]')
        ax.set_ylim([-85-1, -55+1])
        ax.spines['left'].set_bounds([-85, -55])
        ax.set_yticks(np.arange(-85, -55+1, 10))
        ax.set_yticks(np.arange(-85, -55+1, 5), minor = True)  
        
    elif parameter == 'reboundspike_v_amplitude':
        ax.set_ylabel('Reboundspike amplitude [mV]')
        ax.set_ylim([80-1, 130+1])
        ax.spines['left'].set_bounds([80, 130])
        ax.set_yticks(np.arange(80, 130+1, 10))
        ax.set_yticks(np.arange(80, 130+1, 5), minor = True)
        
    elif parameter == 'reboundspike_FWHM':
        ax.set_ylabel('Reboundspike FWHM [ms]')
        ax.set_ylim([0.5-0.1, 2+0.1])
        ax.spines['left'].set_bounds([0.5, 2])
        ax.set_yticks(np.arange(0.5, 2+.1, 0.5))
        ax.set_yticks(np.arange(0.5, 2+.1, 0.25), minor = True)

    # # # frequency adapatation # # #
    elif parameter == 'max_freq':
        ax.set_ylabel('Max. frequency\n(number of spikes) [Hz]')
        ax.set_ylim([0-2, 180+2])
        ax.spines['left'].set_bounds([0, 180])
        ax.set_yticks(np.arange(0, 180+1, 50))
        ax.set_yticks(np.arange(0, 180+1, 10), minor = True)

    elif parameter == 'max_inst_freq':
        ax.set_ylabel('Max. instantaneous\nspiking frequency [Hz]')
        ax.set_ylim([0-2, 180+2])
        ax.spines['left'].set_bounds([0, 180])
        ax.set_yticks(np.arange(0, 180+1, 50))
        ax.set_yticks(np.arange(0, 180+1, 10), minor = True)

    elif parameter == 'max_inst_initial_freq':
        ax.set_ylabel('Max. initial instantaneous\nspiking frequency [Hz]')
        ax.set_ylim([0-2, 180+2])
        ax.spines['left'].set_bounds([0, 180])
        ax.set_yticks(np.arange(0, 180+1, 50))
        ax.set_yticks(np.arange(0, 180+1, 10), minor = True)
        
    elif parameter == 'reboundspike_t_toPeak':
        ax.set_ylabel('Reboundspike\ntime to peak [ms]')
        ax.set_ylim([0.75-0.1, 1.75+0.1])
        ax.spines['left'].set_bounds([0.75, 1.75])
        ax.set_yticks(np.arange(0.75, 1.75+.1, 0.5))
        ax.set_yticks(np.arange(0.75, 1.75+.1, 0.25), minor = True)




# remove grid
[ax.grid(False) for ax in axs_cats]

# despine
[ax.spines[spine].set_visible(False) for ax in axs_cats for spine in ['top', 'right']]
[ax.spines['bottom'].set_bounds([0, 2]) for ax in axs_cats]

# x tick labels
[ax.set_xlabel('') for ax in axs_cats]
[ax.set_xticklabels(['BAOT/\nMeA', 'MeA', 'BAOT'], rotation = 0) for ax in axs_cats[:-4]]
[ax.set_xlim([-0.65, 2.5]) for ax in axs_cats]

# align labels
fig_cats.align_labels() 

# remove unused axes
fig_cats.delaxes(axs_cats[0])
fig_cats.delaxes(axs_cats[4])
fig_cats.delaxes(axs_cats[12])
fig_cats.delaxes(axs_cats[14])
fig_cats.delaxes(axs_cats[16])

# set font size
mtl.rcParams.update({'font.size': 12})

# display figure
plt.show()

temp_fig_dir = "C:/Users/nesseler/Desktop/Poster_iBehave"

save_figures(fig_cats, 'poster_ePhys_ cats', temp_fig_dir, darkmode_bool,
             figure_format= 'both',
             dataframe_to_save = plt_df, index_label = 'cell_ID', add_measures = True, axis_for_calcs = 0,
             groups_bool= True, groups= ['BAOT/MeA', 'MeA', 'BAOT'], groups_name= 'Region')








