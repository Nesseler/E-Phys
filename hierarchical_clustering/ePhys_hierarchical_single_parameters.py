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
        self.ys = np.linspace(padded_min, padded_max, ys_steps)
        
        # test scale factor (bandwidth)
        # bandwidth is calculated by n**(-1./(d+4)), with n datapoints and d dimensions
        # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde
        # density.set_bandwidth(0.4)
        
        # get mean, median, std
        self.mean = np.mean(data)
        self.median = np.median(data)
        self.std = np.std(data)
        
        # get kde with y points
        kde = density(self.ys)
        
        
        self.linewidth = 1
    
        
        self.v_width = 0.25
        self.v_offset = -0.05
        self.e_offset = +0.15
        
        # normalize kde to maximum
        kde_normed = [(k / np.max(kde)) for k in kde]
        
        # remove values below 1 % on edges
        cutoff = 0.01
        xy_s = [[k, y] for k, y in zip(kde_normed, self.ys) if (k > cutoff or y < data_max and y > data_min)]
        
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
        
        self.ax.plot([self.x_position + self.v_offset] * len(self.ys_clipped), self.ys_clipped,
                      color = region_colors[region],
                      linewidth = self.linewidth,
                      label = '_nolegend_',
                      zorder = 2)
        
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



# %% 





### ToDo: function for costum violin

def plot_half_violin(data, ax,
                     data_pad_factor = 1,
                     v_resolution = 0.1,
                     v_kde_cutoff = 0.01,
                     v_abs_cutoff = [np.nan, np.nan],
                     v_position = 0,
                     v_direction = -1,
                     v_offset = -0.05,
                     v_width = 0.25,
                     v_color = 'w',
                     v_lw = 1,
                     v_baseline = False,
                     v_fill = False,
                     v_fillcolor = 'w',
                     v_filllw = 1,
                     v_zorder = 0):
    
    """
    Function plots half violin at specified position.
    Parameters:
        data : Array of data that is used for plotting.
        ax : matplotlib axis to plot violin on.
        data_pad_factor : Float, default is 1. Factor that data_range is multiplied
                          with for estimation of kernel density.
        
        v_resolution : Float (data scale), default is 0.1. Resolution of violin in 
                       data native scale.
        v_kde_cutoff : Float (kde percentage), default is 0.01. Cutoff that is used
                       to limit violin.
        v_abs_cutoff : List of floats (data scale), default is [np.nan, np.nan]. 
                       Absolute values that can be used to limit violin. Shape of 
                       list is [Min, Max]. When only one value is specified, the 
                       other must be set to np.nan and here the v_kde_cutoff is used.
                       If both values are specified the kde_cutoff must be np.nan to 
                       take effect.
        v_position : Float, default is 0. Position of violin on x axis.
        v_direction : -1 or 1, default is -1. Direction factor of violin on x
                      axis.
        v_offset : Float, default is -0.05. Offset of violin from v_position.
        v_width : Float, default is 0.25. Width of violin on x axis
        v_color : Str, default is 'w'. Color of violin.
        v_lw : Float, default is 1. Linewidth of violin.
        v_baseline : Boolean, default is False. Boolean to plot the baseline of
                     the violin.
        v_fill : Boolean, default is False. Boolean to plot the fill of the violin.
        v_fillcolor : Str, default is 'w'. Color of fill.
        v_filllw : Float, default is 1. Linewidth of fill.
        v_zorder : Integer, default is 0. Z-Order of violin.
    
    Returns:
        half_violin : List of Line and Fillobjects. Order is [violin, fill, baseline] 
    """
    
    from scipy.stats import gaussian_kde

    # function inputs
    # data = [5,60,85,100,140,-10,15,20,35,25,15,165,85,75,80,165,15,15,100,25,45,5]
    # ax = plt.gca()
    
    # # defaults
    # data_pad_factor = 1
    # v_resolution = 0.1 # in native data scale
    # v_kde_cutoff = 0.01 # in % default is 0.01
    # v_abs_cutoff = [-50, np.nan] # in native data scale, list of [min, max] for violin cutoff
    
    # v_position = 0
    # v_offset = -0.05
    # v_width = 1
    
    # v_color = 'w'
    # v_lw = 1
    # v_baseline = False
    # v_fill = True
    # v_fillcolor = v_color
    # v_filllw = v_lw
    # v_zorder = 0

    # get distribution metrics
    data_n = len(data)
    data_min = np.min(data)
    data_max = np.max(data)
    data_range = data_max - data_min
    
    # pad data range for violin that extends over range
    padded_min = data_min - (data_range * data_pad_factor)
    padded_max = data_max + (data_range * data_pad_factor)
    
    # get kernel density estimate
    density = gaussian_kde(data)
    
    # scale factor (bandwidth)
    # bandwidth is calculated by n**(-1./(d+4)), with n datapoints and d dimensions
    # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde
    # density.set_bandwidth(0.4)
  
    # get range of value for density estimation
    vs = np.arange(padded_min, padded_max + v_resolution, v_resolution)
    
    # get kde (kernel density esimate)
    kde = density(vs)
    
    # normalize kde to maximum
    kde_normed = [(k / np.max(kde)) for k in kde]
    
    # # # violin cutoff # # #
    
    # condition for set kde cutoff
    if (not np.isnan(v_kde_cutoff) and np.isnan(v_abs_cutoff).all()):
        # remove values below 1 % on edges
        kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff or v < data_max and v > data_min)]
        # print('0')
        
    # condition for set min and max absolute cutoff
    elif (np.isnan(v_kde_cutoff) and not np.isnan(v_abs_cutoff).all()):
        # remove values outside absolute cutoff
        kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (v > v_abs_cutoff[0] and v < v_abs_cutoff[1])]
        # print('1')
      
    # condition for set kde_cutoff and either abs min or max cutoff
    elif (not np.isnan(v_kde_cutoff) and np.isnan(v_abs_cutoff).any()):
        # condition for set min
        if not np.isnan(v_abs_cutoff[0]):
            # remove values below 1 % on edges
            kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff and v > v_abs_cutoff[0])]
            # print('2.1')
            
        # condition for set max 
        elif not np.isnan(v_abs_cutoff[1]):
            # remove values below 1 % on edges
            kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff and v < v_abs_cutoff[1])]
            # print('2.2')
        
    # redefine violin
    kde_clipped = [kv[0] for kv in kv_s]
    vs_clipped = [kv[1] for kv in kv_s]
    
    # add direction, scale, offset ,and x position
    kde_withOffset = [v_direction * k * v_width + v_offset + v_position for k in kde_clipped]
    
    # plot half violin
    half_violin_list = ax.plot(kde_withOffset, vs_clipped,
                               color = v_color,
                               linewidth = v_lw,
                               label = '_nolegend_',
                               zorder = 2)
    
    # violin return objects
    half_violin = half_violin_list
    
    # set coordinates for baseline
    vs_baseline = vs_clipped
    ks_baseline = [v_offset + v_position] * len(vs_clipped)
    
    if v_fill:
        v_face = ax.fill_betweenx(y = vs_clipped,
                                  x1 = kde_withOffset,
                                  x2 = ks_baseline,
                                  color = v_fillcolor,
                                  linewidth = v_filllw)
        
        # append to return object
        half_violin.append(v_face[0])
        
    
    # add endpoints of violin
    vs_baseline.insert(0, vs_baseline[0])
    vs_baseline.insert(-1, vs_baseline[-1])
    ks_baseline.insert(0, kde_withOffset[0])
    ks_baseline.insert(-1, kde_withOffset[-1])
    
    # condition to plot violin baseline
    if v_baseline:    
        baseline_list = ax.plot(ks_baseline, vs_baseline,
                                color = v_color,
                                linewidth = v_lw,
                                label = '_nolegend_',
                                zorder = v_zorder) 
        
        # append to return object
        half_violin.append(baseline_list[0])
    
    return half_violin
    
    
    







# line = plt.plot([0, 1], [0, 1])

data = [5,60,85,100,140,-10,15,20,35,25,15,165,85,75,80,165,15,15,100,25,45,5]
ax = plt.gca()

plot_half_violin(data, ax, v_lw = 0.5, v_direction=-1)

plt.show()