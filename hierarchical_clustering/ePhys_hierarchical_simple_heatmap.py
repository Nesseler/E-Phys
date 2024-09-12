# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:45:38 2024

@author: nesseler
"""

import pandas as pd
import scipy as sc
import numpy as np
from os.path import join

from parameters.directories_win import table_file, hierarchical_dir



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
import seaborn as sbn

from functions.functions_plotting import save_figures, get_colors, get_figure_size


# set colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% z-score celldescriptors

# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()


# %% heatmap

# initialize figure
fig_heat, ax_heat = plt.subplots(nrows = 1,
                                 ncols = 1,
                                 layout = 'constrained',
                                 figsize = get_figure_size(),
                                 dpi = 600,
                                 sharey = True,
                                 sharex = True)

heatmin = -2
heatmax = 2

# plot heatmap
sbn.heatmap(celldescriptors_zscored,
            vmin = heatmin,
            vmax = heatmax,
            center = 0,
            square = False,
            ax = ax_heat,
            xticklabels= 1, 
            cmap = 'icefire', 
            yticklabels = 3,
            linewidth = 0,
            cbar = True) 


# show figure
plt.show()

# set directory for figure
heat_fig_dir = join(hierarchical_dir, 'temp_figs')

# save figure
save_figures(fig_heat, 'figure-hierarchical_clustering-heatmap-unsorted-wo_sag', 
             save_dir = heat_fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')