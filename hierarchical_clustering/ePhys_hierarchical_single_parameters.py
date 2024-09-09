#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 14:37:27 2024

@author: moritznesseler
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

import matplotlib as mtl
import matplotlib.pyplot as plt
import seaborn as sbn

from functions.functions_plotting import save_figures, get_colors, get_figure_size


# set colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% rheobase_current 

import seaborn.objects as so


plt_df = pd.concat([celldescriptors.loc[:, ['rheobase_abs', 'rheobase_rel']], MetaData.loc[celldescriptors.index, 'Region']], axis = 1)

plt_measurement = plt_df['rheobase_abs']


fig_r, axs_r = plt.subplots(nrows = 1,
                            ncols = 3,
                            figsize = get_figure_size(),
                            dpi = 600,
                            width_ratios = [1,1,3],
                            layout = 'constrained')

axs_r[2].set_box_aspect(1)




from scipy.stats import gaussian_kde
from functions.functions_useful import linear_func



class combined_plot:
    
    def __init__(self, data, ax, x_position):
        
        self.data = data
        self.ax = ax
        self.x_position = x_position
    
    
        # get min and max of data
        data_max = np.max(data)
        data_min = np.min(data)
        data_range = data_max - data_min
        self.data_n = len(data)
        
        # get kernel density estimate
        density = gaussian_kde(data)
        
        # set range of y_axis
        padded_min = data_min - data_range
        padded_max = data_max + data_range
        # ys_interval = 1
        ys_steps = 1000
        ys = np.linspace(padded_min, padded_max, ys_steps)
        
        # test scale factor (bandwidth)
        # bandwidth is calculated by n**(-1./(d+4)), with n datapoints and d dimensions
        # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde
        # density.set_bandwidth(0.4)
        
        # get mean, median, std
        self.mean = np.mean(data)
        self.median = np.median(data)
        self.std = np.std(data)
        
        # get kde with y points
        kde = density(ys)
        
        
        self.linewidth = 1
    
        
        self.v_width = 0.25
        self.v_offset = -0.05
        self.e_offset = +0.15
        
        # normalize kde to maximum
        kde_normed = [(k / np.max(kde)) for k in kde]
        
        # remove values below 1 % on edges
        cutoff = 0.01
        xy_s = [[k, y] for k, y in zip(kde_normed, ys) if (k > cutoff or y < data_max and y > data_min)]
        
        # redefine 
        kde_clipped = [xy[0] for xy in xy_s]
        self.ys_clipped = [xy[1] for xy in xy_s]
        
        # add scale, offset ,and x position
        self.kde_withOffset = [- k * self.v_width + self.v_offset + x_position for k in kde_clipped]
        
        # plotting specifications (v - violin, s - swarm, e - errorbar)
        # colors
        self.v_color = region_colors['MeA']
        self.s_color = colors_dict['primecolor']
        # self.e_color = self.v_color
        
        # sizes
        self.s_size = 4
        self.e_msize = 7
        self.e_csize = 3
        self.e_dsize = 7
        
        

        
    def create_plot(self):
        
        # make half violin
        self.ax.plot(self.kde_withOffset, self.ys_clipped,
                     color = self.v_color,
                     linewidth = self.linewidth,
                     label = '_nolegend_',
                     zorder = 2)
        
        # axs_r[0].plot([x_position + v_offset] * len(ys), ys,
        #               color = region_colors[region])
        
        # swarm plot
        self.swarm = sbn.swarmplot(y = self.data,
                                   x = [self.x_position] * self.data_n,
                              ax = self.ax,
                              s = self.s_size,
                              color = self.s_color,
                              zorder = 1)
        
        # get offsets from swarmplot
        self.swarm_offset = self.swarm.collections[0].get_offsets()
        
        # remove swarmplot
        # self.swarm.collections[0].remove()
        
        # print(self.swarm_offset)
        
        # recreate swarmplot as scatterplot with shifted x coordinates
        # self.ax.scatter(x = self.swarm_offset[:,0] + self.x_position,
        #                 y = self.swarm_offset[:,1],
        #                 s = self.s_size,
        #                 color = self.s_color,
        #                 zorder = 1)
        
        # plot errorbar
        self.ax.errorbar(x = self.x_position + self.e_offset,
                         y = self.mean,
                         yerr = self.std,
                         fmt = '_', 
                         markersize = self.e_msize,
                         markerfacecolor = 'none',
                         capsize = self.e_csize,
                         color = self.v_color,
                         linewidth = self.linewidth,
                         label = '_nolegend_',
                         zorder = 3)
        
        self.ax.scatter(x = self.x_position + self.e_offset,
                        y = self.median,
                        marker = 'D', 
                        s = self.e_dsize,
                        color = self.v_color,
                        linewidth = self.linewidth,
                        label = '_nolegend_',
                        zorder = 3)




for r_idx, region in enumerate(['BAOT/MeA', 'MeA', 'BAOT']):
    
    
    # data = plt_df[plt_df['Region'] == region]['rheobase_abs']
    
    
    combinedplot_abs = combined_plot(data = plt_df[plt_df['Region'] == region]['rheobase_abs'], 
                                 ax = axs_r[0], 
                                 x_position = r_idx)
    
    combinedplot_abs.v_color = region_colors[region]
    
    combinedplot_abs.create_plot()
    
    
    # data = plt_df[plt_df['Region'] == region]['rheobase_abs']
    
    
    combinedplot_rel = combined_plot(data = plt_df[plt_df['Region'] == region]['rheobase_rel'], 
                                 ax = axs_r[1], 
                                 x_position = r_idx)
    
    combinedplot_rel.v_color = region_colors[region]
    
    combinedplot_rel.create_plot()



# # # scatter plot # # #

    axs_r[2].scatter(x = plt_df[plt_df['Region'] == region]['rheobase_abs'],
                     y = plt_df[plt_df['Region'] == region]['rheobase_rel'],
                     s = 5,
                     color = region_colors[region],
                     zorder = 2)
    
    axs_r[2].plot([-25, 250], [-25, 250],
                  lw = 0.5,
                  color = 'gray',
                  zorder = 1)





# fit linear curve
popt, pcov = sc.optimize.curve_fit(linear_func, plt_df['rheobase_abs'], plt_df['rheobase_rel'])

fit_dict = {'c' : 'grey', 'ls' : '--', 'lw' : 1}
axs_r[2].plot(np.arange(-25, 250), linear_func(np.arange(-25, 250), *popt), **fit_dict)
    

plt.show()