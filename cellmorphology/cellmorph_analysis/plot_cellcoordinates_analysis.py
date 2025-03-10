# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 11:26:36 2025

@author: nesseler
"""

# init plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *

# init additional packages
import numpy as np
import pandas as pd
from os.path import join, exists
from os import mkdir

# get directories
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_anaylsis

# %% set shared variables

# get field of view
from cellmorphology.cellmorph_parameters import field_of_view, max_depth

# set neurite types
neurite_types = ['neurites', 
                 'dendrites',
                 'axons']


# %% function for multiview of cell

def create_multiview_figure(title):
    
    # set colors 
    from cellmorphology.cellmorph_functions.cellmorph_init_plotting import neurite_color_dict
    neurite_color_dict = neurite_color_dict['all']

    fig, axs = plt.subplots(nrows = 2,
                            ncols = 2,
                            layout = 'constrained',
                            figsize = get_figure_size(width = 100, height = 100),
                            width_ratios = [1, max_depth/field_of_view],
                            height_ratios = [1, max_depth/field_of_view],
                            sharey = 'row',
                            sharex = 'col',
                            dpi = 300)
    
    # flatten axes array
    axs = axs.flatten()
    
    # set figure title
    fig.suptitle(title, 
                 fontsize = 9)
    
    ## legend
    
    # create points for legend
    points = [Line2D([0], [0], 
                     label = ntype, 
                     marker = '.', 
                     markersize = 1.5,  
                     markerfacecolor = neurite_color_dict[ntype], 
                     markeredgecolor = neurite_color_dict[ntype],
                     linestyle='') 
              for ntype in neurite_color_dict.keys()]
    
    # get legend handles
    h, l = axs[2].get_legend_handles_labels()
    
    axs[3].legend(points, neurite_color_dict.keys(), 
                  title = 'neurite type',
                  title_fontsize = 6,
                  frameon = False, 
                  ncol = 1, 
                  loc = 'center',
                  fontsize = 6,
                  markerscale = 5)
    
    axs[3].set_xticklabels(labels = [])
    axs[3].tick_params(axis = 'y', size = 0)
    axs[3].tick_params(axis = 'y', which = 'minor', size = 0)
    axs[3].tick_params(axis = 'x', size = 0)
    axs[3].tick_params(axis = 'x', which = 'minor', size = 0)
    
    # remove spines
    [axs[3].spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom', 'left']]
    
    ## edit axes
    # XY
    axs[0].text(x = 10, 
                y = 10, 
                s = 'XY', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[0].set_xlim([0, field_of_view])
    axs[0].set_ylim([field_of_view, 0])
    axs[0].set_ylabel('Height [µm]')
    axs[0].set_yticks(ticks = np.arange(0, field_of_view, 200))
    axs[0].set_yticks(ticks = np.arange(0, field_of_view, 25), minor = True)
    
    # ZY
    axs[1].text(x = 10, 
                y = 10, 
                s = 'ZY', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[1].set_xlim([0, max_depth])
    axs[1].set_xlabel('')
    axs[1].set_xticks(ticks = np.arange(0, max_depth, 200))
    axs[1].set_xticks(ticks = np.arange(0, max_depth, 25), minor = True)
    
    # XZ
    axs[2].text(x = 10, 
                y = 10, 
                s = 'XZ', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[2].set_xlim([0, field_of_view])
    axs[2].set_xlabel('Width [µm]')
    axs[2].set_ylim([max_depth, 0])
    axs[2].set_ylabel('Depth [µm]')
    axs[2].set_xticks(ticks = np.arange(0, field_of_view, 200))
    axs[2].set_xticks(ticks = np.arange(0, field_of_view, 25), minor = True)
    axs[2].set_yticks(ticks = np.arange(0, max_depth, 200))
    axs[2].set_yticks(ticks = np.arange(0, max_depth, 25), minor = True)
    
    # align labels
    fig.align_labels()

    return fig, axs


# %% plot cell coordinates

# cell_ID = 'E-292'
# cell_coordinates = cell_coordinates
# region = 'BAOT'


def plot_cellcoordinates(cell_ID, cell_coordinates, region = 'BAOT/MeA'):
    '''
    This function creates a figure displaying the coordinates of the reconstructed
    cell in xy, zy, and xz planes. .
    Parameters:
        cell_ID: str, like 'Exp-162', unifque cell identifier
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
        region: str, describes region of cell, default is "BAOT/MeA".
    '''

    # set colors 
    from cellmorphology.cellmorph_functions.cellmorph_init_plotting import neurite_color_dict
    neurite_color_dict = neurite_color_dict['all']
    
    # set title
    title = f'{cell_ID} cell coordinates'
    
    # init standard multiview figure
    fig, axs = create_multiview_figure(title)
    
    # get all path_IDs
    cell_path_IDs = cell_coordinates['path_IDs']
    
    # loop through all paths of cell
    for path_ID in cell_path_IDs:
        
        # get path coordinates
        path_all_coordinates = cell_coordinates['all_coor'][cell_coordinates['all_coor']['path_ID'] == path_ID]          
    
        # get label of current path
        cur_path_label = path_all_coordinates['path_label'].iloc[0]
           
        # scatter plots
        for ax, dim1, dim2, dim3 in zip(axs[:3], ['X', 'Z', 'X'], ['Y', 'Y', 'Z'], ['Z', 'X', 'Y']):
            
            # set zorder of points
            if cur_path_label == 'soma':
                scatter_plot_dict = {'s' : 2,
                                     'zorder' : 600}          
            else:
                # get average depth for zorder
                avg_depth = np.mean(path_all_coordinates[dim3])
                scatter_plot_dict = {'s' : 0.5,
                                     'zorder' : avg_depth}
                
            # all coordinates
            ax.scatter(x = path_all_coordinates[dim1],
                       y = path_all_coordinates[dim2],
                       color = neurite_color_dict[cur_path_label],
                       label = cur_path_label,
                       **scatter_plot_dict)
    
    # get figure path
    figpath = join(cellmorph_anaylsis, 'plots-cell_coordinates')
    
    # save figure
    save_figures(fig, 
                 f'{cell_ID}-cellcoordinates_xyz', 
                 figpath, 
                 darkmode_bool, 
                 figure_format='png')
    
    # display figure
    plt.show()
