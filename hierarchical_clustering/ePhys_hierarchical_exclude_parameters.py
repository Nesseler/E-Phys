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

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# load celldescriptors
celldescriptors = pd.read_excel(join(hierarchical_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')


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
ymax = 1700
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
corr_fig_dir = join(hierarchical_dir, 'temp_figs')

# save figure
save_figures(fig_dist, 'figure-hierarchical_clustering-violinplots-all_parameters', 
             save_dir = corr_fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')


# %%



