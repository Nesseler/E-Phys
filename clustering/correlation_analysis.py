#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 17:11:23 2024

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


# %% correlation analysis

# initialize figure
fig_corr_heat, axs_corr_heat = plt.subplots(nrows = 1,
                                            ncols = 2,
                                            layout = 'constrained',
                                            figsize = get_figure_size(width = 200, height = 100),
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
                            fontsize=7, 
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
            cbar_kws={'label': 'Correlation coefficient',  'ticks' : [-1, -corr_threshold, corr_threshold, 1]},
            ax = axs_corr_heat[1])

# set axis title
axs_corr_heat[1].set_title(r'B: Pearson correlation coefficient $\pm$' + str(corr_threshold),
                            fontsize=7, 
                            loc='left',
                            x = 0)

# show figure
plt.show()

# set directory for figure
corr_fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_corr_heat, f'figure-correlation_matrices-threshold_{str(corr_threshold).replace(".", "p")}', 
             save_dir = corr_fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')