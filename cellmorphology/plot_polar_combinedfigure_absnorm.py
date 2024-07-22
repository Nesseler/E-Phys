# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 17:43:11 2024

@author: nesseler
"""

import pandas as pd
from os import mkdir
from os.path import join, exists
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np



from parameters.directories_win import table_file, cell_morph_descrip_dir, cell_morph_traces_coordinates_dir, cell_morph_plots_dir

from functions.functions_plotting import get_figure_size, get_colors, set_font_sizes, save_figures
from cellmorphology.functions_cellmorph import clean_OnPath_column_to_path_ID_n_label, calc_polar_histo_binangles

# settings
# get colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)


# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')


# %%

cell_IDs = ['E-137']

# cell_IDs = MetaData[MetaData['reconstructed'] == 1].index.to_list()

for cell_ID in cell_IDs:
    
    
    # load dataframes
    
    ### terminal branches measurements
    # diretory
    terminal_branch_measurements_path = join(cell_morph_descrip_dir, 'terminal_branches_measurements', f'{cell_ID}-terminal_branches.xlsx')
    
    # dataframe
    terminal_branches_df = pd.read_excel(terminal_branch_measurements_path, index_col = 'path_ID')
    
    
    ### coordinates
    # all coordinates
    all_coordinates_path = join(cell_morph_traces_coordinates_dir, f'{cell_ID}-all_coordinates.csv')
    cell_allcoordinates = pd.read_csv(all_coordinates_path)
    clean_OnPath_column_to_path_ID_n_label(cell_allcoordinates)
    
    
    # end / last / terminal coordinates
    last_coordinates_path = join(cell_morph_traces_coordinates_dir, f'{cell_ID}-last_coordinates.csv')
    cell_endcoordinates = pd.read_csv(last_coordinates_path)
    clean_OnPath_column_to_path_ID_n_label(cell_endcoordinates)
    
    
    # check if x coordinates need to be flipped
    if MetaData.at[cell_ID, 'to_be_x_flipped']:
        print('flipping x coordinates')
        cell_allcoordinates['X'] = 590.76 - cell_allcoordinates['X']
        cell_endcoordinates['X'] = 590.76 - cell_endcoordinates['X']
    
    
    # get soma coordinates
    soma_coordinates = cell_endcoordinates[cell_endcoordinates['path_ID'] == 1]
    
    
    # get all path_IDs
    path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()
    
    
    # %% combined figure for one cell
    
    fig_all, ax_all = plt.subplots(nrows = 2,
                                   ncols = 2,
                                   layout = 'constrained',
                                   figsize = get_figure_size(width = 165.5))
    
    # set figure title
    fig_all.suptitle(cell_ID)
    
    
    
    
    ### subplot 1: XY coordinates ###
    
    # loop through path
    for path_ID in path_IDs:
        
        # get all coordinates
        path_all_coordinates = cell_allcoordinates[cell_allcoordinates['path_ID'] == path_ID]
        path_end_coordinates = cell_endcoordinates[cell_endcoordinates['path_ID'] == path_ID]           
    
        # get label of current path
        cur_path_label = path_all_coordinates['path_label'].iloc[0]
        
        # set colors for paths
        path_color_dict = {'dendrite' : {'path_color' : colors_dict['primecolor'], 'end_color' : 'r'},
                            'axon' :     {'path_color' : 'b', 'end_color' : 'lightcoral'},
                            'soma' :     {'path_color' : colors_dict['color2'], 'end_color' : colors_dict['color2']}}
    
        path_color = path_color_dict[cur_path_label]['path_color']
        end_color = path_color_dict[cur_path_label]['end_color']
    
    
        
        ax_all.flat[0].scatter(path_all_coordinates['X'], path_all_coordinates['Y'],
                                s = 0.5, c = path_color, label = 'cell')
            
        ax_all.flat[0].scatter(path_end_coordinates['X'], path_end_coordinates['Y'], 
                                s = 0.75, c = end_color, label = 'end points')
        
    ax_all.flat[0].set_box_aspect(1)
    
    ax_all.flat[0].scatter(soma_coordinates['X'], soma_coordinates['Y'], s = 10, color = colors_dict['color2'], label = 'soma')
    ax_all.flat[0].text(x = 10, y = 10, s = 'XY', ha = 'left', va = 'top')
    ax_all.flat[0].set_xlim([0, 590.76])
    ax_all.flat[0].set_ylim([590.76, 0])
    
    ax_all.flat[0].set_xticks(np.arange(0, 590, 200))
    ax_all.flat[0].set_yticks(np.arange(0, 590, 200))
    
    
    
    # define function to change projection type of subplot specific subplot
    def change_projection(fig, axs, ax_tochange, projection = 'polar'):
    
        rows, cols, start, stop = ax_tochange.get_subplotspec().get_geometry()
    
        axs.flat[start].remove()
        axs.flat[start] = fig.add_subplot(rows, cols, start+1, projection=projection)
    
    # change projection of subplots
    [change_projection(fig_all, ax_all, ax, 'polar') for ax in ax_all.flat[1:]]
    
    
    
    ### all polar plot ### 
    
    
    # set colors
    neurite_colordict = {'axon' : 'b', 'dendrite' : colors_dict['primecolor']}
    
    # initilize polar histogram
    n_bins = 8
    binsize = (2 * np.pi) / n_bins
    bins_angles = calc_polar_histo_binangles(n_bins)
    
    # get total number of neurites in bin
    total_n_neurites = terminal_branches_df.drop(index = 1).shape[0]
    total_n_dendrites = terminal_branches_df[terminal_branches_df['path_label'] == 'dendrite'].shape[0]
    total_n_axons = terminal_branches_df[terminal_branches_df['path_label'] == 'axon'].shape[0]
    
    # get max number of neurites in bin
    max_n_neurites = terminal_branches_df.groupby('bin_id').size().max()
    max_n_dendrites = terminal_branches_df[terminal_branches_df['path_label'] == 'dendrite'].groupby('bin_id').size().max()
    max_n_axons = terminal_branches_df[terminal_branches_df['path_label'] == 'axon'].groupby('bin_id').size().max()
    
    
    # axons & dendrites
    
    # define array with number of previouse numbers of branches in bin
    bottom = [0] * n_bins
    
    # skip (drop) path 1, i.e. soma
    branch_idc = terminal_branches_df.drop(index = 1).index.to_list()
    
    # loop through all branches to assign specific color for length            
    for branch_idx in branch_idc:
    
        # get angles of branches
        branch_length = terminal_branches_df.at[branch_idx, "length"]
        branch_bin = terminal_branches_df.at[branch_idx, "bin_id"].astype(int)
        branch_label = terminal_branches_df.at[branch_idx, "path_label"]
        
        # create empty bins and assign branch to bin
        hist_angles_occu = [0] * n_bins
        hist_angles_occu[branch_bin] = 1 / total_n_neurites
        
        # plot histogram as barplot
        ax_all.flat[1].bar(bins_angles, hist_angles_occu, bottom = bottom,
                           width = binsize, 
                           align = 'edge',
                           edgecolor = 'none',
                           color = neurite_colordict[branch_label])
            
        # add to bottom list for next step
        bottom = np.add(bottom, hist_angles_occu)
        
    # redefine bottom as neurite occurrences
    neurite_occu = bottom
    
    # x axis
    ax_all.flat[1].set_xticks(np.arange(0, np.pi*2, np.pi / 4))
    ax_all.flat[1].set_xticklabels(['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv'])
    
    # grid
    ax_all.flat[1].grid(True, alpha = 0.5)
    
    # yaxis
    ax_all.flat[1].set_ylim([0, 1])
    ax_all.flat[1].set_yticks(ticks = [1.0])
    
    
    # colorcoded dendrites or axons
    
    
    ### color code for length of branches ###
    # initialise color code
    norm_min = 0
    norm_max = 1000
    cmap_str = 'inferno'
    norm = mtl.colors.Normalize(norm_min, norm_max)
    cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
    
    # colorbar
    # fig_all.colorbar(cmap, ax = ax_all.flat[3], label = 'Terminal branch length [µm]', orientation='horizontal')
    fig_all.colorbar(cmap, ax = ax_all.flat[3], label = 'Terminal branch length [µm]', orientation='vertical')
    
    
    def plot_colorcoded_polar_normed(polar_occurances_df, max_n_neurites, ax):
    
        # define array with number of previouse numbers of branches in bin
        bottom = [0] * n_bins
        
        # skip (drop) path 1, i.e. soma
        if 1 in polar_occurances_df.index.to_list():
            branch_idc = polar_occurances_df.drop(index = 1).index.to_list()
        else:
            branch_idc = polar_occurances_df.index.to_list()
        
        # loop through all branches to assign specific color for length            
        for branch_idx in branch_idc:
        
            # get angles of branches
            branch_length = polar_occurances_df.at[branch_idx, "length"]
            branch_bin = polar_occurances_df.at[branch_idx, "bin_id"].astype(int)
            
            # create empty bins and assign branch to bin
            hist_angles_occu = [0] * n_bins
            hist_angles_occu[branch_bin] = 1 / max_n_neurites
            
            # plot histogram as barplot
            ax.bar(bins_angles, hist_angles_occu, bottom = bottom,
                    width = binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = cmap.to_rgba(branch_length))
                
            # add to bottom list for next step
            bottom = np.add(bottom, hist_angles_occu)
    
        return bottom
    
    # plot dendrites
    dendrite_occu = plot_colorcoded_polar_normed(terminal_branches_df[terminal_branches_df['path_label'] == 'dendrite'], 
                                                 max_n_neurites= total_n_dendrites, ax = ax_all.flat[2])        
    
    # plot dendrites
    axon_occu = plot_colorcoded_polar_normed(terminal_branches_df[terminal_branches_df['path_label'] == 'axon'], 
                                             max_n_neurites= total_n_axons, ax = ax_all.flat[3])           
    
    # axis for polar plots
    
    for ax in ax_all.flat[1:]:
        
        # x axis
        ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
        ax.set_xticklabels(['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv'])
        
        ax.grid(True, alpha = 0.5)
    
        
    for ax, occu, total_n, max_n in zip(ax_all.flat[1:], [neurite_occu, dendrite_occu, axon_occu], [total_n_neurites, total_n_dendrites, total_n_axons], [max_n_neurites, max_n_dendrites, max_n_axons]):
        # get max occurrences
        max_occurrence = np.max(occu)
        
        # no dendrites / axons
        if total_n == 0:
            total_n = 1
            max_occurrence = 1
        
        # calculate label steps
        step = 1 / total_n
     
        if ((max_n % 2) == 0) and (max_n < 8):
            start = step * 2
            step = step * 2
        elif ((max_n % 2) == 0) and (max_n > 8):
            start = step * 2
            step = step * 4
        elif ((max_n % 2) != 0) and (max_n > 8):
            start = step * 3
            step = step * 4
        else:
            start = step
            step = step *2
            
        yticks = np.arange(start, max_occurrence + step/5, step)
            
        ylabels = [str(round(label_float * 100, 1)) + ' %' for label_float in yticks]
        
        # y axis          
        ax.set_ylim([0, max_occurrence])
        ax.set_yticks(ticks = yticks, labels = ylabels, fontsize = 9)
        ax.set_rlabel_position(80)
    
    
    fig_all.align_labels()
    
    set_font_sizes(small_font_size= 9, large_font_size= 12)
    
    plt.show()
    
    fig_dir = join(cell_morph_plots_dir, 'polar_plots',  'coordinates_neurites_dendrites_axons')
    save_figures(fig_all, f'{cell_ID}-polar_plot-coordinates_neurites_dendrites_axons-normed_to_type', save_dir= fig_dir,
                 figure_format='both')
    
    
    
    
