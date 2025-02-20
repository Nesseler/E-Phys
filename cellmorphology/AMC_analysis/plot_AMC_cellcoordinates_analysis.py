# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:22:43 2025

@author: nesseler
"""

# init plotting
from cellmorphology.cellmorph_functions.initialize_AMC_cellmorph_plotting import *

# init additional packages
import numpy as np
import pandas as pd
from os.path import join


# %% set shared variables

# field of view dimension
max_fov_xy = 590.76
max_fov_z = 300

# set neurite types
neurite_types = ['neurites', 
                 'dendrites', 
                 'glomerular_dendrites', 
                 'nonglomerular_dendrites', 
                 'lateral_dendrites', 
                 'LOTxing_dendrites', 
                 'axons']


# %% cellcoordinates and height, width, depth

def plot_cellcoordinates(cell_ID, cell_coordinates, cell_hwd = None):
    '''
    This function creates a figure displaying the coordinates of the reconstructed
    cell in xy, zy, and xz planes. Function can be used to indicate cell height,
    width, and depth per neurite type as well.
    Parameters:
        cell_ID: str, like 'Exp-162', unifque cell identifier
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
        cell_hwd: pandas Series, describes cell height, widht, and depth per type.
                  default is None.
    '''

    fig, axs = plt.subplots(nrows = 2,
                            ncols = 2,
                            layout = 'constrained',
                            figsize = get_figure_size(width = 100, height = 100),
                            width_ratios = [1, max_fov_z/max_fov_xy],
                            height_ratios = [1, max_fov_z/max_fov_xy],
                            sharey = 'row',
                            sharex = 'col',
                            dpi = 300)

    # flatten axes array
    axs = axs.flatten()
    
    # set figure title
    fig.suptitle(f'{cell_ID} cell coordinates', 
                 fontsize = 9)
        
    
    # get all path_IDs
    cell_path_IDs = cell_coordinates['path_IDs']
    
    # loop through all paths of cell
    for path_ID in cell_path_IDs:
        
        # get path coordinates
        path_all_coordinates = cell_coordinates['all_coor'][cell_coordinates['all_coor']['path_ID'] == path_ID]
        # path_end_coordinates = cell_coordinates['end_coor'][cell_coordinates['end_coor']['path_ID'] == path_ID]           
    
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
            
    # check if height, width and depth were parsed
    if isinstance(cell_hwd, pd.Series):

        # iterate through axes
        for ax, dim, measure in zip(axs[:3], [('x', 'y'), ('z', 'y'), ('x', 'z')], [('width', 'height'), ('depth', 'height'), ('width', 'depth')]):
        
            # loop through neurite types and plot enclosing rectangles
            for ntype in neurite_types:   
                # initialize rectanlge
                rectangle = mtl.patches.Rectangle(xy = (cell_hwd[f'{ntype}-{dim[0]}_min'], cell_hwd[f'{ntype}-{dim[1]}_min']),
                                                  width = cell_hwd[f'{ntype}-{measure[0]}'],
                                                  height = cell_hwd[f'{ntype}-{measure[1]}'],
                                                  facecolor = 'None',
                                                  edgecolor = neurite_color_dict[ntype],
                                                  alpha = 0.5,
                                                  lw = 0.5)
            
                # add rectangle as patch
                ax.add_patch(rectangle)
                
        
        # update measurements label      
        m_label = 'width [µm] - height [µm] - depth [µm]'
        
        for ntype in neurite_types:
            m_label = m_label + f'\n{ntype}: {"{:.2f}".format(cell_hwd[f"{ntype}-width"])} - {"{:.2f}".format(cell_hwd[f"{ntype}-height"])} - {"{:.2f}".format(cell_hwd[f"{ntype}-depth"])}'
        
        # add measurements to plot
        axs[0].text(x = 580, y = 580, 
                    s = m_label, 
                    ha = 'right', 
                    va = 'bottom',
                    size = 4)
    
    # edit axes
    # XY
    axs[0].text(x = 10, 
                y = 10, 
                s = 'XY', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[0].set_xlim([0, max_fov_xy])
    axs[0].set_ylim([max_fov_xy, 0])
    axs[0].set_ylabel('Height [µm]')
    axs[0].set_yticks(ticks = np.arange(0, max_fov_xy, 200))
    axs[0].set_yticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)
    
    # ZY
    axs[1].text(x = 10, 
                y = 10, 
                s = 'ZY', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[1].set_xlim([0, max_fov_z])
    axs[1].set_xlabel('')
    axs[1].set_xticks(ticks = np.arange(0, max_fov_z, 200))
    axs[1].set_xticks(ticks = np.arange(0, max_fov_z, 25), minor = True)
    
    # XZ
    axs[2].text(x = 10, 
                y = 10, 
                s = 'XZ', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[2].set_xlim([0, max_fov_xy])
    axs[2].set_xlabel('Width [µm]')
    axs[2].set_ylim([max_fov_z, 0])
    axs[2].set_ylabel('Depth [µm]')
    axs[2].set_xticks(ticks = np.arange(0, max_fov_xy, 200))
    axs[2].set_xticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)
    axs[2].set_yticks(ticks = np.arange(0, max_fov_z, 200))
    axs[2].set_yticks(ticks = np.arange(0, max_fov_z, 25), minor = True)
    
    # legend
    h, l = axs[0].get_legend_handles_labels()
    
    # get unique lables
    set_label = ['soma', 'glomerular_dendrites', 'nonglomerular_dendrites', 'lateral_dendrites', 'LOTxing_dendrites', 'axons']
    
    # get and sort indices of unique lables
    label_indices = [l.index(u_l) for u_l in set_label if u_l in l]
    
    # reassign handels
    unique_handles = [h[u_li] for u_li in label_indices]
    unique_label = [l[u_li] for u_li in label_indices]
    
    axs[3].legend(unique_handles, unique_label, 
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
    
    # align labels
    fig.align_labels()
    
    # create saving path and save
    from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_analysis_dir
    
    # adjust directories for saving
    if isinstance(cell_hwd, pd.Series):
        figname = f'{cell_ID}-cell_hwd_xyz'
        figpath = join(AMCs_analysis_dir, 'height_width_depth_plots')
    else:
        figname = f'{cell_ID}-cellcoordinates_xyz'
        figpath = join(AMCs_analysis_dir, 'cell_coordinates_plots')
    
    save_figures(fig, 
                 figname, 
                 figpath, 
                 darkmode_bool, 
                 figure_format='png')
    
    # display figure
    plt.show()


