# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 16:55:47 2025

@author: nesseler
"""

from functions.initialize_plotting import *

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_syn_dir

# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol

# get cell_IDs
cell_IDs = get_cell_IDs_one_protocol('cc_sag', sheet_name = 'PGFs_Syn')

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# limit cell_IDs
cell_IDs = MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list()

# load passive properties
passive_properties = pd.read_excel(join(cell_descrip_syn_dir, 'cc_sag-syn-passive_properties.xlsx'), index_col = 'cell_ID')

# define regions
regions = ['BAOT/MeA', 'MeA', 'BAOT']

# %% figure settings

swarm_dict = {'size' : 4,
              'linewidth' : 0.5,
              'edgecolor' : colors_dict['seccolor']}


# figure

# init figure and axes
fig, axs = plt.subplots(nrows = 1,
                        ncols = 3,
                        layout = 'constrained',
                        figsize = get_figure_size(width = 250, height = 80),
                        dpi = 300)

for ax, prop in zip(axs, passive_properties.columns):
    
    for pos, region in enumerate(regions):
        
        # get cell_IDs for region
        region_cell_IDs = MetaData[MetaData['Region'] == region].index.to_list()
    

    
        data = passive_properties.loc[region_cell_IDs, prop]
        x = [pos] * len(passive_properties.loc[region_cell_IDs, prop])
    
        # plot data points as swarm
        sbn.swarmplot(y = data,
                      x = x,
                      ax = ax,
                      color = region_colors[region],
                      **swarm_dict)
        
        # width relative to total number of cell
        v_width = (len(region_cell_IDs) / passive_properties.shape[0]) *3
        
        # plot half violin
        plot_half_violin(data = data.to_list(), 
                         ax = ax,
                         v_direction = 1,
                         v_resolution = 0.1,
                         v_kde_cutoff = 0.1,
                         v_abs_cutoff = [0, np.nan],
                         v_position = 3.75,
                         v_offset = 0.175 * 1,
                         v_width = v_width,
                         v_baseline = False,
                         v_color = region_colors[region],
                         v_zorder = 0)  
    
        e_position = 3 + (0.25 * pos)
    
        # plot errorbar
        ax.errorbar(x = e_position,
                    y = data.mean(),
                    yerr = data.std(),
                    fmt='_', 
                    markersize = 7,
                    markerfacecolor = 'none',
                    capsize = 3,
                    color = region_colors[region],
                    linewidth = 1,
                    label = '_nolegend_',
                    zorder = 2)
        
        # plot median
        ax.scatter(x = e_position,
                   y = data.median(),
                   marker='D', 
                   s = 11,
                   color = region_colors[region],
                   linewidth = 1,
                   label = '_nolegend_',
                   zorder = 3)
        
        apply_axis_settings(ax, axis = 'y', **axis_dict[prop])
        
        # x
        ax.set_xlim([0 - 0.75, 5 + 0.5])
        ax.set_xticks(ticks = [0, 1, 2], labels = regions)
        ax.set_xticks(ticks = [], minor = True)
        ax.spines['bottom'].set_bounds([0, 2])
        ax.set_xlabel('')
    
        
        
        
        
# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# align labels
fig.align_labels()

# create saving path and save
# from parameters.directories_win import tempfigs_dir
# save_figures(fig, 'cc_sag-sagdelta', tempfigs_dir, darkmode_bool, figure_format='png')

# display figure
plt.show()