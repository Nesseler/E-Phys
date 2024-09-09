#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 17:11:23 2024

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


# # drop columns
# celldescriptors.drop(columns = parameters_toDrop, inplace = True)


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


# %%

# # # correlation analysis # # #


# initialize figure
fig_corr_heat, axs_corr_heat = plt.subplots(nrows = 1,
                                            ncols = 2,
                                            layout = 'constrained',
                                            figsize = get_figure_size(),
                                            dpi = 600,
                                            sharey = True,
                                            sharex = True)

# set aspect ratio
[ax.set_box_aspect(1) for ax in axs_corr_heat]

# set string for color map
cmap_str = 'seismic'

# get correlation matrix of dataframe
celldescriptors_corr = celldescriptors.corr()

# plot heatmap of correlation values
sbn.heatmap(data = celldescriptors_corr,
            cmap = cmap_str,
            vmin = -1,
            vmax = 1,
            annot = False,
            cbar_kws={'label': 'Correlation coefficient', 'ticks' : [-1, 0, 1]},
            ax = axs_corr_heat[0])

# set axis title
axs_corr_heat[0].set_title('A: Pearson correlation coefficient',
                            fontsize=12, 
                            loc='left',
                            x = 0)


# plot heatmap with values above threshold only
corr_threshold = 0.8

# replace values in dataframe below threshold with nan
celldescriptors_corr_thresh = celldescriptors_corr[(celldescriptors_corr > corr_threshold) | (celldescriptors_corr < -corr_threshold)]

# plot heatmap of correlation values
sbn.heatmap(data = celldescriptors_corr_thresh,
            cmap = cmap_str,
            vmin = -1,
            vmax = 1,
            annot = False,
            cbar_kws={'label': 'Correlation coefficient',  'ticks' : [-1, -0.8, 0.8, 1]},
            ax = axs_corr_heat[1])

# set axis title
axs_corr_heat[1].set_title(r'B: Pearson correlation coefficient $\pm$' + str(corr_threshold),
                            fontsize=12, 
                            loc='left',
                            x = 0)

# show figure
plt.show()

# set directory for figure
corr_fig_dir = join(hierarchical_dir, 'temp_figs')

# save figure
save_figures(corr_fig_dir, 'figure-hierarchical_clustering-correlation_matrices-all_parameter', 
             save_dir = corr_fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')