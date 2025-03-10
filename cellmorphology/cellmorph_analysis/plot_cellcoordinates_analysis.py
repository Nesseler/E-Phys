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
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_anaylsis_dir

# %% set shared variables

# get field of view
from cellmorphology.cellmorph_parameters import field_of_view, max_depth

# get orientation labels
from cellmorphology.cellmorph_functions.cellmorph_polarplot_functions import orientation_labels

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


def plot_cellcoordinates(cell_ID, cell_coordinates):
    '''
    This function creates a figure displaying the coordinates of the reconstructed
    cell in xy, zy, and xz planes. .
    Parameters:
        cell_ID: str, like 'E-137', unifque cell identifier
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
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
    figpath = join(cellmorph_anaylsis_dir, 'plots-cell_coordinates')
    
    # save figure
    save_figures(fig, 
                 f'{cell_ID}-cellcoordinates_xyz', 
                 figpath, 
                 darkmode_bool, 
                 figure_format='png')
    
    # display figure
    plt.show()


# %% plot height, width, depth

def plot_cellhwd(cell_ID, cell_coordinates, cell_hwd):
    '''
    This function creates a figure displaying the coordinates of the reconstructed
    cell in xy, zy, and xz planes. Function can be used to indicate cell height,
    width, and depth per neurite type as well.
    Parameters:
        cell_ID: str, like 'E-137', unifque cell identifier
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
        cell_hwd: pandas Series, describes cell height, widht, and depth per type.
    '''

    # set colors 
    from cellmorphology.cellmorph_functions.cellmorph_init_plotting import neurite_color_dict
    neurite_color_dict = neurite_color_dict['all']
    
    # set title
    title = f'{cell_ID} cell height, width, depth'
    
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
    
    # get figure path
    figpath = join(cellmorph_anaylsis_dir, 'plots-height_width_depth')
    
    # save figure
    save_figures(fig, 
                 f'{cell_ID}-cell_hwd_xyz', 
                 figpath, 
                 darkmode_bool, 
                 figure_format='png')
    
    # display figure
    plt.show()


# %% plot terminal branches

def plot_all_terminal_branches(cell_ID, cell_coordinates, terminal_branches):
    '''
    This function creates figures for each terminal branch displaying the 
    coordinates of the reconstructed terminal branch, with its measurements 
    of length and orientation.
    Parameters:
        cell_ID: str, like 'E-137', unifque cell identifier
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
        terminal_branches: pandas Dataframe, containing the measurements of all
                           terminal branches
    '''

    # set colors 
    from cellmorphology.cellmorph_functions.cellmorph_init_plotting import neurite_color_dict
    neurite_color_dict = neurite_color_dict['all']
    
    # group paths by label
    allcoor_paths_dict = {pathID : group for pathID, group in cell_coordinates['all_coor'].groupby('path_ID')}
    
    for terminal_path_ID in terminal_branches.index.to_list():
    
        # get coordinates
        pathIDs_2soma = terminal_branches.at[terminal_path_ID, 'pathIDs_2soma']
        terminal_coor = terminal_branches.loc[terminal_path_ID, ['x_end', 'y_end', 'z_end']]
        cell_somacoordinates = allcoor_paths_dict[1].loc[0, :]
    
        # set title
        title = f'{cell_ID} terminal path-{terminal_path_ID}'
    
        # init standard multiview figure
        fig, axs = create_multiview_figure(title)
    
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
                    f'branch contraction = {round(terminal_branches.at[terminal_path_ID ,"contraction"], 2)}',
                    f'bin id = {terminal_branches.at[terminal_path_ID ,"bin_id"]}',
                    f'orientation = {orientation_labels[terminal_branches.at[terminal_path_ID ,"bin_id"]]}',
                    ]
        
        # branch measurement
        axs[0].text(x = 580, y = 580, 
                    s = '\n'.join(m_label), 
                    ha = 'right', 
                    va = 'bottom',
                    size = 4)
            
        # create figure directory
        fig_dir = join(cellmorph_anaylsis_dir, 'plots-terminal_branches', f'{cell_ID}-terminal_branches')
        
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
        
        
# %% plot primary and terminal points


def plot_endpoints(cell_ID, cell_coordinates, n_primary, n_terminal, bifurcation_ratios):
    '''
    This function creates a figure displaying the coordinates of the reconstructed
    cell in xy and the corresponding primary and terminal end points.
    Parameters:
        cell_ID: str, like 'E-137', unifque cell identifier
        cell_coordinates: dict, includes all-, end-, soma coordinates dataframes
                          as well as a list of path_IDs
        n_primary: pandas Dataframe, containing the number of primary end points
        n_terminal: pandas Dataframe, containing the number of terminal end points
        bifurcation_ratios: pandas Dataframe, containing the bifurcation ratios
    '''

    # set colors 
    from cellmorphology.cellmorph_functions.cellmorph_init_plotting import neurite_color_dict
    neurite_color_dict = neurite_color_dict['all']

    # init figure
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
    ax.set_ylim([0, field_of_view])
    ax.set_ylabel('Height [µm]')
    ax.set_yticks(ticks = np.arange(0, field_of_view, 200))
    ax.set_yticks(ticks = np.arange(0, field_of_view, 25), minor = True)
    
    ax.set_xlim([0, field_of_view])
    ax.set_xlabel('Width [µm]')
    ax.set_xticks(ticks = np.arange(0, field_of_view, 200))
    ax.set_xticks(ticks = np.arange(0, field_of_view, 25), minor = True)
    
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
    ax_inset.set_xticks(ticks = np.arange(0, field_of_view, 200), labels = [])
    ax_inset.set_xticks(ticks = np.arange(0, field_of_view, 25), labels = [], minor = True)
    ax_inset.set_xlim([box_xmin, box_xmin + box_width])
    
    # y
    ax_inset.set_yticks(ticks = np.arange(0, field_of_view, 200), labels = [])
    ax_inset.set_yticks(ticks = np.arange(0, field_of_view, 25), labels = [], minor = True)
    ax_inset.set_ylim([box_ymin, box_ymin + box_height])
    ax_inset.invert_yaxis()
    
    # update measurements label      
    m_label = 'n_primary - n_terminal - bifurcation_ratios'
    
    for ntype in neurite_types:
        m_label = m_label + f'\n{ntype}: {"{:.2f}".format(n_primary.at[cell_ID, f"n_primary-{ntype}"])} - {"{:.2f}".format(n_terminal.at[cell_ID, f"n_terminal-{ntype}"])} - {"{:.2f}".format(bifurcation_ratios.at[cell_ID, f"bifurcation_ratio-{ntype}"])}'
    
    # add measurements to plot
    ax.text(x = 830, y = 580, 
            s = m_label, 
            ha = 'right', 
            va = 'bottom',
            size = 4)
    
    # get figure path
    figpath = join(cellmorph_anaylsis_dir, 'plots-primary_terminal_bifurcation')
    
    # save figure
    save_figures(fig, 
                 f'{cell_ID}-primary_terminal_bifurcation', 
                 figpath, 
                 darkmode_bool, 
                 figure_format='png')
    
    # display figure
    plt.show()


# %% cell polar plot (absolute)

def plot_polar_plot_abs(cell_ID, terminal_branches):
    '''
    This function creates an polar histogram of the orientation of all terminal
    branches for one cell.
    Parameters:
        cell_ID: str, like 'E-137', unifque cell identifier
        terminal_branches: pandas Dataframe, containing the measurements of all
                           terminal branches
    '''

    from cellmorphology.cellmorph_functions.cellmorph_polarplot_functions import create_polar_histogram
    
    # sort dataframe of terminal branches measurements to plot histogram
    terminal_branches.sort_values('length', inplace = True)
    
    # init figure
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},
                           layout = 'constrained',
                           height_ratios= [1],
                           width_ratios=[1],
                           figsize = get_figure_size(width = 100, height = 100),
                           dpi = 300)
    
    # set title
    fig.suptitle(cell_ID)
    
    # create polar plot on axis
    create_polar_histogram(ax, terminal_branches)
    
    # set grid
    ax.grid(True, alpha = 0.5)
    
    # set grid behind plot
    ax.set_axisbelow(True)
    
    # get legend handles
    h, l = ax.get_legend_handles_labels()
    
    # get unique labels and handles
    unique_l = np.unique(l)
    unique_h = [h[l.index(ntype)] for ntype in unique_l]
    
    # legend
    ax.legend(unique_h, unique_l, 
              title = 'neurite type',
              title_fontsize = 6,
              frameon = False, 
              ncol = 2, 
              loc = 'upper center',
              fontsize = 6,
              bbox_to_anchor=(0.5, -0.08))
    
    
    # get figure path
    figpath = join(cellmorph_anaylsis_dir, 'plots-polar_plot_abs')
    
    # save figure
    save_figures(fig, 
                 f'{cell_ID}-polar_plot_abs', 
                 figpath, 
                 darkmode_bool, 
                 figure_format='png')
    
    # display figure
    plt.show()
