#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:34:59 2024

@author: moritznesseler
"""

import matplotlib.pyplot as plt
import seaborn as sbn


# rewrite function
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