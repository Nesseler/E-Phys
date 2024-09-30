# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 17:36:46 2024

@author: nesseler
"""

import pandas as pd
import scipy as sc
import numpy as np
from os.path import join

from parameters.directories_win import table_file, hierarchical_dir

from hierarchical_clustering.ePhys_hierarchical_parameters import parameters_toDrop

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# load celldescriptors
celldescriptors = pd.read_excel(join(hierarchical_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')


# %% initialize plotting

import matplotlib.pyplot as plt
import seaborn as sbn

from functions.functions_plotting import save_figures, get_colors, get_figure_size, plot_half_violin, apply_axis_settings
from parameters.axis_settings import ePhys_parameter_dict

# set colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
import matplotlib as mtl
mtl.rcParams.update({'font.size': 9})

# set region list
regions = ['BAOT/MeA', 'MeA', 'BAOT']


# %%

from functions.functions_useful import linear_func

# set parameters to plot
parameter1 = 'max_inst_freq'
parameter2 = 'max_inst_initial_freq'

# create dataframe for plotting
plt_df = pd.concat([celldescriptors.loc[:, [parameter1, parameter2]], MetaData.loc[celldescriptors.index, 'Region']], axis = 1)

# initialize figure
fig, axs = plt.subplots(nrows = 1,
                        ncols = 4,
                        figsize = get_figure_size(),
                        dpi = 600,
                        width_ratios = [1, 1, 2, 1],
                        layout = 'constrained')

# edit constrained layout
fig.get_layout_engine().set(wspace=0.10)

# set labels list
labels = ['A', 'B', 'C']


for p_idx, parameter in enumerate([parameter1, parameter2]):
    
    # set axis
    ax = axs[p_idx]
    
    # set axis title
    ax.set_title(f'{labels[p_idx]}: {parameter}',
                 fontsize=12, 
                 loc='left',
                 x = 0)
    
    # swarmplots
    swarm = sbn.swarmplot(data = plt_df,
                          x = 'Region',
                          y = parameter, 
                          ax = ax,
                          hue = 'Region', 
                          palette = region_colors,
                          size = 3,
                          order = regions,
                          zorder = 1)
    
    # remove legend
    # ax.legend().set_visible(False)
    
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
        
        # plot half violin
        plot_half_violin(data = violin_data, 
                         ax = ax,
                         v_kde_cutoff = 0.025,
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

        
        # remove spines
        [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
         
        # apply changes to y axis
        apply_axis_settings(ax, **ePhys_parameter_dict[parameter])
        
        # x tick labels
        ax.set_xticks(ticks = [0, 1, 2], labels = ['BAOT/\nMeA', 'MeA', 'BAOT'], rotation = 0)
        ax.set_xlim([-0.65, 2.5])
        ax.spines['bottom'].set_bounds([0, 2])
        

        
# # # scatter plot & correlation # # #

# set axis
ax = axs[2]

# set axis title
ax.set_title(f'{labels[2]}: Linear correlation',
             fontsize=12, 
             loc='left',
             x = 0)

# get axis settings
xdict = ePhys_parameter_dict[parameter1]
ydict = ePhys_parameter_dict[parameter2]

# initialise dict for fitting parameters
fits_dict = dict()

# iterate through regions
for r_idx, region in enumerate(regions):
    
    # get data
    scatter_data = plt_df[plt_df['Region'] == region]

    # plot scatter plot
    ax.scatter(x = scatter_data[parameter1],
               y = scatter_data[parameter2],
               s = 10,
               color = region_colors[region],
               zorder = 2,
               label = region)
    
    # fit linear curve
    popt, pcov = sc.optimize.curve_fit(linear_func, scatter_data[parameter1], scatter_data[parameter2])
    
    # write to dict
    fits_dict[region] = {'popt' : popt,
                         'pcov' : pcov}

    # plot linear fit
    x_linear_fit = np.arange(xdict['ax_min'], xdict['ax_max'] + xdict['stepminor'], xdict['stepminor'])
    
    ax.plot(x_linear_fit, 
            linear_func(x_linear_fit, *popt), 
            c = region_colors[region],
            lw = 1,
            linestyle = 'dashed',
            label = f'linear fit ({region})')
    
    
# calc and plot linear fit for all
popt, pcov = sc.optimize.curve_fit(linear_func, plt_df[parameter1], plt_df[parameter2])

# write to dict
fits_dict['all'] = {'popt' : popt,
                    'pcov' : pcov}

# plot linear fit
x_linear_fit = np.arange(xdict['ax_min'], xdict['ax_max'] + xdict['stepminor'], xdict['stepminor'])

ax.plot(x_linear_fit, 
        linear_func(x_linear_fit, *popt), 
        c = colors_dict['primecolor'],
        lw = 1,
        linestyle = 'dashed',
        label = 'linear fit (all)')


# apply axis settings
# x
apply_axis_settings(ax = ax, axis = 'x', **xdict)

# y
apply_axis_settings(ax = ax, axis = 'y', **ydict)

# plot unity line
ax.plot([xdict['ax_min'], xdict['ax_max']], [ydict['ax_min'], ydict['ax_max']],
        lw = 1,
        linestyle = 'dotted',
        color = 'gray',
        alpha = 0.5,
        zorder = 1,
        label = 'unityline')

# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# # # subplot for figure legend and fit values # # #

# get legend handles and labels
handles, labels = axs[2].get_legend_handles_labels()

# legend
legend = axs[3].legend(handles, labels,
                       loc='lower center',
                       ncol=1,
                       fontsize = 9)

# remove legend border
axs[3].get_legend().get_frame().set_linewidth(0.0)

# remove spines and ticks
[axs[3].spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom', 'left']]
axs[3].set_xticks([])
axs[3].set_yticks([])


# align labels
fig.align_labels()
        
# show plot
plt.show()

# set directory for figure
fig_dir = join(hierarchical_dir, 'temp_figs')

# save figure
save_figures(fig, f'figure-hierarchical_clustering-two_parameters-{parameter1}_v_{parameter2}', 
             save_dir = fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')