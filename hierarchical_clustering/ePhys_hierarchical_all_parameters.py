#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:24:55 2024

@author: moritznesseler
"""

import pandas as pd
import scipy as sc
import numpy as np
from os.path import join

from parameters.directories_win import table_file, cell_descrip_dir, hierarchical_dir

from hierarchical_clustering.ePhys_hierarchical_parameters import parameters_toDrop

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# load celldescriptors
celldescriptors = pd.read_excel(join(hierarchical_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')

# drop columns
from hierarchical_clustering.ePhys_hierarchical_parameters import parameters_toDrop
celldescriptors.drop(columns = parameters_toDrop, inplace = True)


# %% initialize plotting

import matplotlib as mtl
import matplotlib.pyplot as plt

from functions.functions_plotting import save_figures, get_colors, get_figure_size


# set colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% 

from hierarchical_clustering.ePhys_hierarchical_plotting_functions import plot_data_distribution



# initialize figure
fig_dist, axs_dist = plt.subplots(nrows = 1,
                                  ncols = 1,
                                  layout = 'constrained',
                                  figsize = get_figure_size(),
                                  dpi = 600)

# plot
plot_data_distribution(axs_dist, distributions_data = celldescriptors, colors_dict = colors_dict)

# edit axis
# x
xmin = 0
xmax = len(celldescriptors.columns) - 1
xpad = 0.6


# y
ymin = -200
ymax = 1500
ypad = 100
ystep = 500
ystepminor = 100
xticklabels = celldescriptors.columns.to_list()
ylabel = 'Parameter value [mV / pA / pF / Hz / ms / # / ]'
axistitle = 'Original parameters'
    
# set axis title
axs_dist.set_title(axistitle,
                   fontsize=12, 
                   loc='left',
                   x = 0)
    
# apply axis changes 
axs_dist.set_xlim([xmin - xpad, xmax + xpad])
axs_dist.set_xticklabels(labels = xticklabels,
                    rotation = 90)
axs_dist.spines['bottom'].set_bounds([xmin, xmax])

axs_dist.set_xlabel('')

axs_dist.set_ylim([ymin - ypad, ymax + ypad])
axs_dist.set_yticks(ticks = np.arange(0, ymax+ ystepminor, ystep))
axs_dist.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystepminor), minor = True)
axs_dist.spines['left'].set_bounds([ymin, ymax])
axs_dist.set_ylabel(ylabel)

# remove spines
[axs_dist.spines[spine].set_visible(False) for spine in ['top', 'right']]
   
# align labels
fig_dist.align_labels()
    
# show figure
plt.show()

# set directory for figure
fig_dist_dir = join(hierarchical_dir, 'temp_figs')

# save figure
save_figures(fig_dist, 'figure-hierarchical_clustering-violinplots-wo_sag', 
             save_dir = fig_dist_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')




# %% z-score celldescriptors

# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()


# %%

# initialize figure
fig_zdist, ax_zdist = plt.subplots(nrows = 1,
                                   ncols = 1,
                                   layout = 'constrained',
                                   figsize = get_figure_size(),
                                   dpi = 600)

# plot
plot_data_distribution(ax_zdist, distributions_data = celldescriptors_zscored, colors_dict = colors_dict)

# edit axis
# x
xmin = 0
xmax = len(celldescriptors.columns) - 1
xpad = 0.6


# y
ymin = -4
ymax = 6
ypad = 0.5
ystep = 2
ystepminor = 1
xticklabels = celldescriptors.columns.to_list()
ylabel = 'Z-transformed parameter value [std]'
axistitle = 'Z-transformed parameters'
    
# set axis title
ax_zdist.set_title(axistitle,
                   fontsize=12, 
                   loc='left',
                   x = 0)
    
# apply axis changes 
ax_zdist.set_xlim([xmin - xpad, xmax + xpad])
ax_zdist.set_xticklabels(labels = xticklabels,
                    rotation = 90)
ax_zdist.spines['bottom'].set_bounds([xmin, xmax])

ax_zdist.set_xlabel('')

ax_zdist.set_ylim([ymin - ypad, ymax + ypad])
ax_zdist.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystep))
ax_zdist.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystepminor), minor = True)
ax_zdist.spines['left'].set_bounds([ymin, ymax])
ax_zdist.set_ylabel(ylabel)

# remove spines
[ax_zdist.spines[spine].set_visible(False) for spine in ['top', 'right']]
   
# align labels
fig_zdist.align_labels()
    
# show figure
plt.show()

# set directory for figure
fig_zdist_dir = join(hierarchical_dir, 'temp_figs')

# save figure
save_figures(fig_zdist, 'figure-hierarchical_clustering-violinplots-z_scored-wo_sag', 
             save_dir = fig_zdist_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')