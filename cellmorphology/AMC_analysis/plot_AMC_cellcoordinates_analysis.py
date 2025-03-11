# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:22:43 2025

@author: nesseler
"""

# init plotting
from cellmorphology.AMC_analysis.initialize_AMC_cellmorph_plotting import *

# init additional packages
import numpy as np
import pandas as pd
from os.path import join, exists
from os import mkdir


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
        figpath = join(AMCs_analysis_dir, 'plots-height_width_depth')
    else:
        figname = f'{cell_ID}-cellcoordinates_xyz'
        figpath = join(AMCs_analysis_dir, 'plots-cell_coordinates')
    
    save_figures(fig, 
                 figname, 
                 figpath, 
                 darkmode_bool, 
                 figure_format='png')
    
    # display figure
    plt.show()


# %% primary, terminal and bifurcation


def plot_endpoints(cell_ID, cell_coordinates, n_primary, n_terminal, n_stems, bifurcation_ratios):
    '''
    This function creates a figure displaying the coordinates of the reconstructed
    cell in xy and the corresponding primary and terminal end points.
    Parameters:
        cell_ID: str, like 'Exp-162', unifque cell identifier
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
        n_primary: pandas Dataframe, containing the number of primary end points
        n_terminal: pandas Dataframe, containing the number of terminal end points
        n_stems: pandas Dataframe, containing the number of stem points
        bifurcation_ratios: pandas Dataframe, containing the bifurcation ratios
    '''

    fig, ax = plt.subplots(nrows = 1,
                           ncols = 1,
                           layout = 'constrained',
                           figsize = get_figure_size(width = 150, height = 100),
                           dpi = 300)
    
    # set figure title
    fig.suptitle(f'{cell_ID} primary and terminal points', 
                 fontsize = 9)
    
    # set aspect ration of plot
    ax.set_aspect(1)
    
    # plot all cell coordinates
    ax.scatter(x = cell_coordinates['all_coor'].loc[:, 'X'],
               y = cell_coordinates['all_coor'].loc[:, 'Y'],
               color = 'gray',
               s = 0.25,
               alpha = 0.1)
    
    # plot primary points (end points of primary paths)
    for path_i in cell_coordinates['pri_coor'].index.to_list():
        ax.scatter(x = cell_coordinates['pri_coor'].at[path_i, 'X'],
                   y = cell_coordinates['pri_coor'].at[path_i, 'Y'],
                   color = neurite_color_dict[cell_coordinates['pri_coor'].at[path_i, 'path_label']],
                   s = 25,
                   marker = 'x',
                   linewidths=0.5)
    
    # plot terminal points
    for path_i in cell_coordinates['end_coor'].index.to_list():
        ax.scatter(x = cell_coordinates['end_coor'].at[path_i, 'X'],
                   y = cell_coordinates['end_coor'].at[path_i, 'Y'],
                   color = neurite_color_dict[cell_coordinates['end_coor'].at[path_i, 'path_label']],
                   s = 5)
    
    # plot soma on top
    ax.scatter(x = cell_coordinates['all_coor'].loc[0, 'X'],
               y = cell_coordinates['all_coor'].loc[0, 'Y'],
               color = neurite_color_dict['soma'])
    
    
    # plane label
    ax.text(x = 10, 
            y = 10, 
            s = 'XY', 
            ha = 'left', 
            va = 'top', 
            fontsize = 9)
    
    # edit axes
    ax.set_ylim([0, max_fov_xy])
    ax.set_ylabel('Height [µm]')
    ax.set_yticks(ticks = np.arange(0, max_fov_xy, 200))
    ax.set_yticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)
    
    ax.set_xlim([0, max_fov_xy])
    ax.set_xlabel('Width [µm]')
    ax.set_xticks(ticks = np.arange(0, max_fov_xy, 200))
    ax.set_xticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)
    
    # invert y axis
    ax.invert_yaxis()
    
    # soma inset
    
    # inset marker
    box_xmin   = cell_coordinates['soma_coor'].at[0,'X']-40
    box_width  = 80 
    box_ymin   = cell_coordinates['soma_coor'].at[0,'Y']-40
    box_height = 80
       
    # add rectangle marker
    ax.add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                           width = box_width, 
                           height = box_height,
                           fill = False,
                           color = primecolor,
                           linestyle = '--',
                           lw = 0.5,
                           alpha = 0.5))
    
    ## ([left, bottom, width, height]), percentages
    ax_inset = ax.inset_axes([1.05, 0.63, 0.35, 0.35],
                              xlim=(box_xmin, box_xmin+box_width), 
                              ylim=(box_ymin, box_ymin+box_height), 
                              xticklabels=[], 
                              yticklabels=[])
    
    # edit linewidth of inset axis and its ticks
    [ax_inset.spines[spine].set_linewidth(0.5) for spine in ['left', 'bottom']]
    ax_inset.tick_params(width=0.5)
    ax_inset.tick_params(which = 'minor', width=0.25)
    
    
    # plot inset
    # plot all cell coordinates
    ax_inset.scatter(x = cell_coordinates['all_coor'].loc[:, 'X'],
                     y = cell_coordinates['all_coor'].loc[:, 'Y'],
                     color = 'gray',
                     s = 0.25,
                     alpha = 0.1)
    
    # plot primary points (end points of primary paths)
    for path_i in cell_coordinates['pri_coor'].index.to_list():
        ax_inset.scatter(x = cell_coordinates['pri_coor'].at[path_i, 'X'],
                         y = cell_coordinates['pri_coor'].at[path_i, 'Y'],
                       color = neurite_color_dict[cell_coordinates['pri_coor'].at[path_i, 'path_label']],
                       s = 25,
                       marker = 'x',
                       linewidths=0.5)
    
    # plot terminal points
    for path_i in cell_coordinates['end_coor'].index.to_list():
        ax_inset.scatter(x = cell_coordinates['end_coor'].at[path_i, 'X'],
                         y = cell_coordinates['end_coor'].at[path_i, 'Y'],
                      color = neurite_color_dict[cell_coordinates['end_coor'].at[path_i, 'path_label']],
                      s = 5)
        
    # plot soma on top
    ax_inset.scatter(x = cell_coordinates['all_coor'].loc[0, 'X'],
                     y = cell_coordinates['all_coor'].loc[0, 'Y'],
                     color = neurite_color_dict['soma'])
       
    # x
    ax_inset.set_xticks(ticks = np.arange(0, max_fov_xy, 200), labels = [])
    ax_inset.set_xticks(ticks = np.arange(0, max_fov_xy, 25), labels = [], minor = True)
    ax_inset.set_xlim([box_xmin, box_xmin + box_width])
    
    # y
    ax_inset.set_yticks(ticks = np.arange(0, max_fov_xy, 200), labels = [])
    ax_inset.set_yticks(ticks = np.arange(0, max_fov_xy, 25), labels = [], minor = True)
    ax_inset.set_ylim([box_ymin, box_ymin + box_height])
    ax_inset.invert_yaxis()
    
    # update measurements label      
    m_label = 'n_primary - n_terminal - n_stems - bifurcation_ratios'
    
    for ntype in neurite_types:
        m_label = m_label + f'\n{ntype}: {"{:.2f}".format(n_primary.at[cell_ID, f"n_primary-{ntype}"])} - {"{:.2f}".format(n_terminal.at[cell_ID, f"n_terminal-{ntype}"])} - {"{:.2f}".format(n_stems.at[cell_ID, f"n_stems-{ntype}"])} - {"{:.2f}".format(bifurcation_ratios.at[cell_ID, f"bifurcation_ratio-{ntype}"])}'
    
    # add measurements to plot
    ax.text(x = 830, y = 580, 
            s = m_label, 
            ha = 'right', 
            va = 'bottom',
            size = 4)
    
    # create saving path and save
    from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_analysis_dir
    
    # save figure
    save_figures(fig, 
                 f'{cell_ID}-primary_terminal_bifurcation', 
                 join(AMCs_analysis_dir, 'plots-primary_terminal_bifurcation'), 
                 darkmode_bool, 
                 figure_format='png')
        
    # display plot
    plt.show()


# %% terminal branches

from matplotlib.lines import Line2D

def plot_all_terminal_branches(cell_ID, terminal_branches, cell_coordinates):
    '''
    This function creates figures for each terminal branch displaying the 
    coordinates of the reconstructed terminal branch, with its measurements 
    of length and orientation.
    Parameters:
        cell_ID: str, like 'Exp-162', unifque cell identifier
        terminal_branches: pandas Dataframe, containing the measurements of all
                           terminal branches
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
    '''
    
    # group paths by label
    allcoor_paths_dict = {pathID : group for pathID, group in cell_coordinates['all_coor'].groupby('path_ID')}

    for terminal_path_ID in terminal_branches.index.to_list():
    
        # get coordinates
        pathIDs_2soma = terminal_branches.at[terminal_path_ID, 'pathIDs_2soma']
        terminal_coor = terminal_branches.loc[terminal_path_ID, ['x_end', 'y_end', 'z_end']]
        cell_somacoordinates = allcoor_paths_dict[1].loc[0, :]
    
        # initialize figure
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
        fig.suptitle(f'{cell_ID} terminal path-{terminal_path_ID}',
                     fontsize = 9)
        
        # scatter plots
        for ax, dim1, dim2 in zip(axs[:3], ['X', 'Z', 'X'], ['Y', 'Y', 'Z']):
    
            for path_ID in pathIDs_2soma:
                # get path label
                path_label = allcoor_paths_dict[path_ID]['path_label'].unique()[0]
                
                # plot coordinates
                ax.scatter(x = allcoor_paths_dict[path_ID][dim1], 
                           y = allcoor_paths_dict[path_ID][dim2],
                           s = 0.25, 
                           zorder = 1,
                           color = neurite_color_dict[path_label],
                           label = '_nolegend_')
            
            # plot terminal point
            ax.scatter(x = terminal_coor[dim1.lower() + '_end'], 
                       y = terminal_coor[dim2.lower() + '_end'],
                       s = 6,
                       marker = 'x',
                       zorder = 2,
                       color = primecolor,
                       label = '_nolegend_',
                       lw = 0.5)     
            
            # plot soma
            ax.scatter(x = cell_somacoordinates[dim1], 
                       y = cell_somacoordinates[dim2],
                       s = 30,
                       marker = '.',
                       zorder = 2,
                       color = neurite_color_dict['soma'],
                       label = '_nolegend_',
                       lw = 1)    
            
        # measurments
        m_label = [f'branch angle [deg] = {round(terminal_branches.at[terminal_path_ID ,"angle_deg"], 2)}',
                   f'branch angle [rad] = {round(terminal_branches.at[terminal_path_ID ,"angle_rad"], 2)}',
                   f'branch length [µm] = {round(terminal_branches.at[terminal_path_ID ,"length"], 2)}',
                   f'branch euclidean dist. [µm] = {round(terminal_branches.at[terminal_path_ID ,"euc_dist"], 2)}',
                   f'branch contraction = {round(terminal_branches.at[terminal_path_ID ,"contraction"], 2)}']
        
        # branch measurement
        axs[0].text(x = 580, y = 580, 
                    s = '\n'.join(m_label), 
                    ha = 'right', 
                    va = 'bottom',
                    size = 4)
        
        points = [Line2D([0], [0], 
                         label = ntype, 
                         marker = '.', 
                         markersize = 1.5,  
                         markerfacecolor = neurite_color_dict[ntype], 
                         markeredgecolor = neurite_color_dict[ntype],
                         linestyle='') 
         for ntype in neurite_color_dict.keys()]
    
        # legend
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
        
        # align labels
        fig.align_labels()
    
    
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
    
    
        # create saving path and save
        from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_analysis_dir
        
        # create figure directory
        fig_dir = join(AMCs_analysis_dir, 'plots-terminal_branches', f'{cell_ID}-terminal_branches')
        
        # test if folder exists and create new one, if necessary
        if not exists(fig_dir):
            mkdir(fig_dir)
        
        # save figure
        save_figures(fig, 
                     f'{cell_ID}-terminal_branch-{terminal_path_ID}', 
                     fig_dir, 
                     darkmode_bool, 
                     figure_format='png')
        
        # display plot
        plt.show()
        

# %% create figure function

def create_multiview_figure(title):
    
    # initialize figure
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
    fig.suptitle(title,
                 fontsize = 9)
    
    points = [Line2D([0], [0], 
                     label = ntype, 
                     marker = '.', 
                     markersize = 1.5,  
                     markerfacecolor = neurite_color_dict[ntype], 
                     markeredgecolor = neurite_color_dict[ntype],
                     linestyle='') 
     for ntype in neurite_color_dict.keys()]

    # legend
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
    
    # align labels
    fig.align_labels()


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
    
    return fig, (axs)
    

# %% neurite types

def plot_neurite_types(cell_ID, cell_coordinates):
    '''
    This function creates figures for each branch type and plots the 
    coordinates of the reconstructed branches.
    Parameters:
        cell_ID: str, like 'Exp-162', unifque cell identifier
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
    '''

    # create dictionary grouped by path label
    cell_pathcoordinates_dict = {ntype : group for ntype, group in cell_coordinates['all_coor'].groupby('path_label')}
    
    # iterate through label types
    for nlabel in cell_pathcoordinates_dict.keys():
        
        if nlabel != 'soma':
    
            # create figure
            fig, axs = create_multiview_figure(title = f'{cell_ID} neurite types')
            
        
            # scatter plots
            for ax, dim1, dim2 in zip(axs[:3], ['X', 'Z', 'X'], ['Y', 'Y', 'Z']):
                
                # plot all cell coordinates
                ax.scatter(x = cell_coordinates['all_coor'].loc[:, dim1],
                            y = cell_coordinates['all_coor'].loc[:, dim2],
                            color = 'gray',
                            s = 0.25,
                            alpha = 0.1,
                            zorder = 0)
                
                # get coordinates
                cell_pathcoordinates = cell_pathcoordinates_dict[nlabel]
                
                # plot all cell coordinates
                ax.scatter(x = cell_pathcoordinates.loc[:, dim1],
                           y = cell_pathcoordinates.loc[:, dim2],
                           color = neurite_color_dict[nlabel],
                           s = 0.25,
                           alpha = 0.1,
                           zorder = 1)
                
                # plot terminal points
                for path_i in cell_coordinates['end_coor'][cell_coordinates['end_coor']['path_label'] == nlabel].index.to_list():
                    ax.scatter(x = cell_coordinates['end_coor'].at[path_i, dim1],
                               y = cell_coordinates['end_coor'].at[path_i, dim2],
                               color = neurite_color_dict[nlabel],
                               s = 20,
                               marker = '.',
                               lw = 0.5,
                               zorder = 2)
                    
                # plot soma
                ax.scatter(x = cell_pathcoordinates_dict['soma'][dim1], 
                           y = cell_pathcoordinates_dict['soma'][dim2],
                           s = 30,
                           marker = '.',
                           zorder = 3,
                           color = neurite_color_dict['soma'],
                           label = '_nolegend_',
                           lw = 1)    
            
            # create saving path and save
            from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_analysis_dir
    
            # save figure
            save_figures(fig, 
                          f'{cell_ID}-neurite_types-{nlabel}', 
                          join(AMCs_analysis_dir,'plots-neurite_types'), 
                          darkmode_bool, 
                          figure_format='png')
            
            # display figure   
            plt.show()