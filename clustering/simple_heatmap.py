# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:45:38 2024

@author: nesseler
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


# %% z-score celldescriptors

# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()

# celldescriptors_zscored.sort_values(by = ['r_input'], inplace = True)


# %% initialize plotting

from functions.initialize_plotting import *


# %% heatmap

# initialize figure
fig_heat, ax_heat = plt.subplots(nrows = 1,
                                 ncols = 1,
                                 layout = 'constrained',
                                 figsize = get_figure_size(),
                                 dpi = 600,
                                 sharey = True,
                                 sharex = True)

heatmin = -2 #celldescriptors_zscored.min().min()
heatmax = 2  #celldescriptors_zscored.max().max()

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
heat_fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_heat, 'celldescriptors-heatmap-unsorted', 
             save_dir = heat_fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')