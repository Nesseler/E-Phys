#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:24:55 2024

@author: moritznesseler
"""

from functions.initialize_packages import *

# load celldescriptors
from parameters.directories_win import clustering_dir
celldescriptors = pd.read_excel(join(clustering_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')

# get cell_IDs
cell_IDs = celldescriptors.index.to_list()

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)


# %% initialize plotting

from functions.initialize_plotting import *


# %% plotting function

def plot_data_distribution(ax, distributions_data, colors_dict):
        
    # melt dataframe to plot
    distributions_data_melted = distributions_data.melt(var_name = 'parameter')
    
    
    # create list of parameters
    parameters = distributions_data.columns.to_list()
    
    # get number of parameters
    n_parameters = len(parameters)
    
    # create list of colors for parameters
    p_cmap = plt.get_cmap('viridis', n_parameters)
    
    # plot violins
        
    violins = sbn.violinplot(data = distributions_data_melted,
                              x = 'parameter',
                              y = 'value',
                              ax = ax,
                              linewidth = 1,
                              inner = None,
                              density_norm = 'width',
                              gap = 0.25,
                              width = 0.4,
                              hue = True, hue_order=[True, False], split = True)
    
    
    for p_idx, param in enumerate(parameters):
        # set line color of quarts
        # for l_idx in range(3):
        #     violins.lines[p_idx * 3 + l_idx].set_color(p_cmap(p_idx))
        
        # set edge color of violin
        violins.collections[p_idx].set_edgecolor(p_cmap(p_idx))
    
        # set facecolor of violin
        violins.collections[p_idx].set_facecolor('None')
    
    # plot swarm
    
    swarms = sbn.swarmplot(data = distributions_data_melted,
                            x = 'parameter',
                            y = 'value', 
                            ax = ax,
                            s = 1,
                            color=colors_dict['primecolor'])
    
    # plot error bar
    for p_idx, param in enumerate(parameters):
        ax.errorbar(x = p_idx+0.3,
                    y = distributions_data[param].mean(),
                    yerr = distributions_data[param].std(),
                    fmt='_', 
                    markersize = 4,
                    markerfacecolor = 'none',
                    capsize = 1,
                    color=p_cmap(p_idx),
                    linewidth = 1,
                    label = '_nolegend_')
        
    # edit seaborn legend
    ax.legend().set_visible(False)


# %% figure parameter space

# initialize figure
fig_dist, axs_dist = plt.subplots(nrows = 1,
                                  ncols = 1,
                                  layout = 'constrained',
                                  figsize = get_figure_size(),
                                  dpi = 600)

# plot
plot_data_distribution(axs_dist, distributions_data = celldescriptors, colors_dict = colors_dict)


# edit axis
# y
ydict = {'ax_min' : -200,
         'ax_max' : 1500,
         'pad' : None,
         'step' : 500,
         'stepminor' : 100,
         'label' : 'Parameter value [mV / pA / pF / Hz / ms / # / ]'}

apply_axis_settings(axs_dist, axis = 'y', **ydict)

# x
xdict = {'ax_min' : 0,
         'ax_max' : len(celldescriptors.columns) - 1,
         'pad' : None,
         'step' : 1,
         'stepminor' : 1,
         'label' : '',
         'ticklabels' : celldescriptors.columns.to_list(),
         'rotation' : 60}

apply_axis_settings(axs_dist, axis = 'x', **xdict)

# remove spines
[axs_dist.spines[spine].set_visible(False) for spine in ['top', 'right']]
   
# align labels
fig_dist.align_labels()
    
# show figure
plt.show()

# set directory for figure
fig_dist_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_dist, 'figure-violinplots-parameter_space', 
             save_dir = fig_dist_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')




# %% z-score celldescriptors

# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()


# %% figure z-scored parameter sapce

# initialize figure
fig_zdist, ax_zdist = plt.subplots(nrows = 1,
                                    ncols = 1,
                                    layout = 'constrained',
                                    figsize = get_figure_size(),
                                    dpi = 600)

# plot
plot_data_distribution(ax_zdist, distributions_data = celldescriptors_zscored, colors_dict = colors_dict)

# edit axis
# y
ydict = {'ax_min' : -6,
         'ax_max' : 6,
         'pad' : 0.2,
         'step' : 2,
         'stepminor' : 0.5,
         'label' : 'z-transformed parameter value [std]'}

apply_axis_settings(ax_zdist, axis = 'y', **ydict)

apply_axis_settings(ax_zdist, axis = 'x', **xdict)

# remove spines
[ax_zdist.spines[spine].set_visible(False) for spine in ['top', 'right']]
   
# align labels
fig_zdist.align_labels()
    
# show figure
plt.show()

# set directory for figure
fig_zdist_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_zdist, 'figure-violinplots-zscored_parameter_space', 
              save_dir = fig_zdist_dir,
              darkmode_bool = darkmode_bool,
              figure_format = 'png')


# %% figure z-scored parameter space

# initialize figure
fig_zdist, ax_zdist = plt.subplots(nrows = 1,
                                    ncols = 1,
                                    layout = 'constrained',
                                    figsize = get_figure_size(),
                                    dpi = 600)

# plot
plot_data_distribution(ax_zdist, distributions_data = celldescriptors_zscored, colors_dict = colors_dict)

# edit axis
# y
ydict = {'ax_min' : -2.5,
         'ax_max' : 2.5,
         'pad' : 0.2,
         'step' : 2,
         'stepminor' : 0.5,
         'label' : 'z-transformed parameter value [std]',
         'limits_n_0' : True}

apply_axis_settings(ax_zdist, axis = 'y', **ydict)

apply_axis_settings(ax_zdist, axis = 'x', **xdict)

# remove spines
[ax_zdist.spines[spine].set_visible(False) for spine in ['top', 'right']]
   
# align labels
fig_zdist.align_labels()
    
# show figure
plt.show()

# set directory for figure
fig_zdist_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_zdist, 'figure-violinplots-zscored_parameter_space-zoom', 
              save_dir = fig_zdist_dir,
              darkmode_bool = darkmode_bool,
              figure_format = 'both')