# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 13:53:45 2024

@author: nesseler
"""

import pandas as pd
from os import mkdir
from os.path import join, exists
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import math

from parameters.directories_win import cell_morph_traces_coordinates_dir, cell_morph_plots_dir, table_file, cell_morph_descrip_dir
from getter.get_onlyfiles_list import get_onlyfiles_list

from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size
from cellmorphology.functions_cellmorph import clean_OnPath_column_to_path_ID_n_label

# get onlyfiles list
onlyfiles = get_onlyfiles_list(cell_morph_traces_coordinates_dir)

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# extract only filenames with all coordinates
onlyfiles_allcoor = [f for f in onlyfiles if 'all_coordinates' in f]
onlyfiles_endcoor = [f for f in onlyfiles if 'last_coordinates' in f]

vplots_bool = False

cell_IDs = []


# dataframes for exporting
# write dataframe that contains readout from polar plot
polar_plot_occurrances = pd.DataFrame(columns = cell_IDs, index = ['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv'])



# test 
# onlyfiles = onlyfiles[:2]

for cell_in_files in range(int(len(onlyfiles)/2)):

    all_coor_filename = onlyfiles_allcoor[cell_in_files]
    
    # get cell ID from filename
    cell_ID = all_coor_filename[:5]
    cell_IDs.append(cell_ID)
    print(f'Started: {cell_ID}')
    
    # read all coordinates file
    all_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, all_coor_filename)) 
    end_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, onlyfiles_endcoor[cell_in_files])) 
    
    # clean up dataframes
    clean_OnPath_column_to_path_ID_n_label(all_coordinates)
    clean_OnPath_column_to_path_ID_n_label(end_coordinates)
    
    # check if x coordinates need to be flipped
    if MetaData.at[cell_ID, 'to_be_x_flipped']:
        print('flipping x coordinates')
        all_coordinates['X'] = 590.76 - all_coordinates['X']
        end_coordinates['X'] = 590.76 - end_coordinates['X']
    
    # get total number of paths
    n_paths = all_coordinates['path_ID'].max()
    
    # get colors
    darkmode_bool = True
    colors_dict, region_colors = get_colors(darkmode_bool)
    
    
    # %% plot cell with all branches and end points marked
    
    if vplots_bool:
    
        ratio = 300 / 590.76
        
        fig_cell, axs_cell = plt.subplots(nrows = 2,
                                              ncols = 2,
                                              height_ratios= [1, ratio],
                                              width_ratios= [1, ratio],
                                              sharey = 'row',
                                              sharex = 'col',
                                              figsize = get_figure_size(width = 165.5))
        set_font_sizes()
        plt.subplots_adjust(wspace=0, hspace=0)
        axs_cell = axs_cell.flatten()
        fig_cell.delaxes(axs_cell[-1])
        fig_cell.suptitle(cell_ID)
        
        # XY
        axs_cell[0].scatter(all_coordinates['X'], all_coordinates['Y'], s = 0.5, color = colors_dict['primecolor'], label = 'cell')
        axs_cell[0].scatter(end_coordinates['X'], end_coordinates['Y'], s = 0.75, color = 'r', label = 'end points')
        axs_cell[0].scatter(end_coordinates['X'][0], end_coordinates['Y'][0], s = 10, color = colors_dict['color2'], label = 'soma')
        axs_cell[0].text(x = 10, y = 10, s = 'XY', ha = 'left', va = 'top')
        axs_cell[0].set_ylim([590.76, 0])
        axs_cell[0].set_ylabel('Height [µm]')
        axs_cell[0].legend(fontsize = 9)
        
        # ZY
        axs_cell[1].scatter(all_coordinates['Z'], all_coordinates['Y'], s = 0.5, color = colors_dict['primecolor'])
        axs_cell[1].scatter(end_coordinates['Z'], end_coordinates['Y'], s = 0.75, color = 'r')
        axs_cell[1].scatter(end_coordinates['Z'][0], end_coordinates['Y'][0], s = 10, color = colors_dict['color2'])
        axs_cell[1].text(x = 10, y = 10, s = 'ZY', ha = 'left', va = 'top')
        axs_cell[1].set_xlim([0, 300])
        axs_cell[1].tick_params(axis = 'y', size = 0)
        axs_cell[1].set_xlabel('Depth [µm]')
        
        # XZ
        axs_cell[2].scatter(all_coordinates['X'], all_coordinates['Z'], s = 0.5, color = colors_dict['primecolor'])
        axs_cell[2].scatter(end_coordinates['X'], end_coordinates['Z'], s = 0.75, color = 'r')
        axs_cell[2].scatter(end_coordinates['X'][0], end_coordinates['Z'][0], s = 10, color = colors_dict['color2'])
        axs_cell[2].text(x = 10, y = 10, s = 'XZ', ha = 'left', va = 'top')
        axs_cell[2].set_xlim([0, 590.76])
        axs_cell[2].set_ylim([300, 0])
        axs_cell[2].set_ylabel('Depth [µm]')
        axs_cell[2].set_xlabel('Width [µm]')
        
        [ax.grid(False) for ax in axs_cell]
        
        cell_fig_dir = join(cell_morph_plots_dir, 'cell_coordinates')
        save_figures(fig_cell, f'{cell_ID}-cell_coordinates_xyz', cell_fig_dir, darkmode_bool)
        
        plt.show()    
    
    ### calculate xy angle from end point ###
    
    # get soma coordinates
    soma_coordinates = end_coordinates[end_coordinates['path_ID'] == 1]
    
    # define vector to compare to
    reference_vector = pd.Series({'X' : 1., 'Y' : 0.})
    
    # define dataframe for all vectors to endpoints
    terminal_branches_df = pd.DataFrame(soma_coordinates.rename(columns = {'X': 'x_end', 'Y': 'y_end', 'Z': 'z_end'}).set_index('path_ID'))
    
    
    def dotproduct(v1, v2):
      return sum((a*b) for a, b in zip(v1, v2))
    
    def length(v):
      return math.sqrt(dotproduct(v, v))
    
    def angle(v1, v2):
      return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
    
    
    
    # skip first to avoid the soma
    for idx_end_point in end_coordinates.index[1:]:
        
        path_ID_end_point = end_coordinates.at[idx_end_point, 'path_ID']
        
        end_point_coordinates = end_coordinates[end_coordinates['path_ID'] == path_ID_end_point]
    
        vector_coordinates = pd.concat([soma_coordinates, end_point_coordinates], axis = 0)
        vector_coordinates.drop(columns = ['path_ID', 'path_label'], inplace = True)
        
        xy_vector = vector_coordinates.drop(columns = ['Z'])
        
        
        terminal_branch_df = pd.DataFrame({'x_end': end_point_coordinates.at[idx_end_point,'X'],
                                           'y_end': end_point_coordinates.at[idx_end_point,'Y'],
                                           'z_end': end_point_coordinates.at[idx_end_point,'Z'],
                                           'x_diff': end_point_coordinates.at[idx_end_point, 'X'] - soma_coordinates.at[0, 'X'],
                                           'y_diff': end_point_coordinates.at[idx_end_point, 'Y'] - soma_coordinates.at[0, 'Y'],
                                           'z_diff': end_point_coordinates.at[idx_end_point, 'Z'] - soma_coordinates.at[0, 'Z'],
                                           'path_label': end_point_coordinates.at[idx_end_point,'path_label']},
                                          index = [path_ID_end_point])
        
        terminal_branch_df.index.name = 'path_ID'
        
        xy_diff_vector = terminal_branch_df.loc[path_ID_end_point, 'x_diff':'y_diff']
        
        
        if xy_diff_vector['y_diff'] > 0:
            deg = 360 - math.degrees(angle(reference_vector, xy_diff_vector))
        else:
            deg = math.degrees(angle(reference_vector, xy_diff_vector))
    
        terminal_branch_df.at[path_ID_end_point, 'angle_deg'] = deg
        terminal_branch_df.at[path_ID_end_point, 'angle_rad'] = math.radians(deg)
    
        # write to dataframe
        terminal_branches_df = pd.concat([terminal_branches_df, terminal_branch_df], axis = 0)
        


    
    # %% functions for branch reconstruction
    
    
    def node_in_path(first_node_coor, pot_parent_path_coor):
        # create mask where elements are True that correspond to the first node corrdinates
        coor_mask = pot_parent_path_coor == first_node_coor
        
        # get intersection point
        # where all coordinates are the same
        intersect_mask = coor_mask.query('X == True & Y == True & Z == True')
        
        # test if parent, i.e.: mask does or doesn't contain True values
        if intersect_mask.empty:
            in_path_bool = False
            
        else:
            in_path_bool = True
            
        return in_path_bool
    
    
    def find_parent_path(path_ID, path_IDs_to_search, all_coordinates):
        # get path coordinates
        path_coor = all_coordinates[all_coordinates['path_ID'] == path_ID]
    
        # get first node
        firstnode = path_coor.iloc[0, :]
            
        # loop through path_IDs_to_search
        for pot_parent_path_ID in path_IDs_to_search:
            
            # skip branch itself 
            if pot_parent_path_ID != path_ID:
        
                # coordinates of potential parent path
                pot_parent_path_coor = all_coordinates[all_coordinates['path_ID'] == pot_parent_path_ID]
    
                # test if node is in potential parent path 
                if node_in_path(firstnode, pot_parent_path_coor):
                    
                    # create mask where elements are True that correspond to the first node corrdinates
                    coor_mask = pot_parent_path_coor == firstnode
       
                    # get intersection point
                    # where all coordinates are the same
                    intersect_mask = coor_mask.query('X == True & Y == True & Z == True')
       
                    # get intersection index
                    intersect_index = intersect_mask.index[0]
       
                    # get intersection coordinates
                    intersect_coor = pot_parent_path_coor.loc[intersect_index]
       
                    # get all coordinates until intersection point
                    ## get all indices of parent path until intersection index
                    parent_indices = [i for i in pot_parent_path_coor.index.to_list() if i <= intersect_index] 
                    
                    # avoid finding the root of another bifurcatio
                    if len(parent_indices) > 1 or intersect_coor['path_ID'].astype(int) == 1:
                        parent_path_ID = intersect_coor['path_ID'].astype(int)
        
        return parent_path_ID
    
    
    # start with terminal branch coordinates
    # path_ID_terminal_branch = 3
    # terminal_path_coor = all_coordinates[all_coordinates['path_ID'] == path_ID_terminal_branch]    
    # terminal_path_firstnode = terminal_path_coor.iloc[0, :]
        
    
    
    
    
    # print(find_parent_path(13, path_IDs, all_coordinates))
    
    # %%
    
    def calc_length_of_branch(branch_coor):
        """ function calculates the length of an branch by the sum of all 3d
        euclidean distances
        https://en.wikipedia.org/wiki/Euclidean_distance
        """
    
        # calculate the differences between points on each axis
        diff = branch_coor.diff(axis = 0)
        
        # calculate the square of given distance
        sq_diff = diff**2
        
        # calculate sum for each point to get distance
        sum_sq_diff = sq_diff.sum(axis = 1)
        
        # calculate the square root of the sum of each difference squared
        sqrt_sum_sq_diff = np.sqrt(sum_sq_diff)
        
        # calculate length as sum of all euclidean distances
        length = sqrt_sum_sq_diff.sum()
        
        return length
    
    
    # %% reconstruct branch terminal to soma
    
    ## list of all paths except current terminal path
    path_IDs = np.arange(1, n_paths + 1, 1)
    
    terminal_path_ID = 24
    
    terminal_path_IDs = end_coordinates['path_ID'][end_coordinates['path_ID'] != 1]
    
    for terminal_path_ID in terminal_path_IDs:
            
        # coordinates of terminal branch to reconstruct until soma
        terminal_path_coor = all_coordinates[all_coordinates['path_ID'] == terminal_path_ID]
        
        # order of coordinates starts with terminal branch endpoint 
        # to be reversed later when calculating the entire branch
        branch_coor_terminal_to_soma = pd.DataFrame()
        
        # reversed and reindex coordinates of terminal path
        rr_terminal_path_coor = terminal_path_coor[::-1].reset_index(drop=True)
        
        # concatenate coordinates of terminal path
        branch_coor_terminal_to_soma = pd.concat([branch_coor_terminal_to_soma, rr_terminal_path_coor])
        
        # start while loop from terminal path
        path_ID = terminal_path_ID
        parent_path_ID = terminal_path_ID
        
        while parent_path_ID != 1:

            parent_path_ID = find_parent_path(path_ID, path_IDs, all_coordinates)
             
            # coordinates of potential parent path
            parent_path_coor = all_coordinates[all_coordinates['path_ID'] == parent_path_ID]
            
            # get coordinates of first node
            first_node = all_coordinates[all_coordinates['path_ID'] == path_ID].iloc[0, :]
            
            # create mask where elements are True that correspond to the first node corrdinates
            coor_mask = parent_path_coor == first_node
            
            # get intersection point
            # where all coordinates are the same
            intersect_mask = coor_mask.query('X == True & Y == True & Z == True')
            
            # get intersection index
            intersect_index = intersect_mask.index[0].item()
            
            # get intersection coordinates
            intersect_coor = parent_path_coor.loc[intersect_index]
            
            # get all coordinates until intersection point
            ## get all indices of parent path until intersection index
            parent_indices = [i for i in parent_path_coor.index.to_list() if i <= intersect_index]
            
            # get coordinates
            parent_path_coor_until_intersect = all_coordinates.loc[parent_indices]
            
            # reversed and reindex
            rr_parent_path_coor_until_intersect = parent_path_coor_until_intersect[::-1].reset_index(drop=True)
            
            # concatenate to all coordinates from terminal to soma
            branch_coor_terminal_to_soma = pd.concat([branch_coor_terminal_to_soma, rr_parent_path_coor_until_intersect])
            
            # reset index following concatenation
            branch_coor_terminal_to_soma.reset_index(drop=True, inplace = True)
            
            # set path_ID to parent for next step in while loop
            path_ID = parent_path_ID
        
        
        ### calculate length
        terminal_branch_length = calc_length_of_branch(branch_coor_terminal_to_soma.drop(columns = ['path_ID', 'path_label']))
        
        terminal_branches_df.at[terminal_path_ID ,'length'] = terminal_branch_length
           
        ### calculate euclidean distance between origin and terminal of branch
        ## can be calculated with same function but just two coordinates
        line_coordinates = branch_coor_terminal_to_soma.iloc[[0, -1]]
        
        terminal_branch_euc = calc_length_of_branch(line_coordinates.drop(columns = ['path_ID', 'path_label']))
        
        # write to dataframe
        terminal_branches_df.at[terminal_path_ID ,'euc_dist'] = terminal_branch_euc
        
        ### calculate branch contraction
        terminal_branch_contraction = terminal_branch_euc / terminal_branch_length
        
        # write to dataframe
        terminal_branches_df.at[terminal_path_ID ,'contraction'] = terminal_branch_contraction
        
        if vplots_bool:
        
            ### branch verification plot ###
            ratio = 300 / 590.76
            
            fig_branch, axs_branch = plt.subplots(nrows = 2,
                                                  ncols = 2,
                                                  height_ratios= [1, ratio],
                                                  width_ratios= [1, ratio],
                                                  sharey = 'row',
                                                  sharex = 'col',
                                                  figsize = get_figure_size(width = 165.5))
            
            
            
            plt.subplots_adjust(wspace=0, hspace=0)
            
            
            axs_branch = axs_branch.flatten()
            fig_branch.delaxes(axs_branch[-1])
            
            fig_branch.suptitle(f'{cell_ID} terminal_path_ID: {terminal_path_ID}')
            
            ## XY
            # branch
            axs_branch[0].scatter(branch_coor_terminal_to_soma['X'],
                                  branch_coor_terminal_to_soma['Y'],
                                  s = 0.5,
                                  color = colors_dict['primecolor'])
            
            # soma
            axs_branch[0].scatter(end_coordinates[end_coordinates['path_ID'] == 1]['X'], 
                                  end_coordinates[end_coordinates['path_ID'] == 1]['Y'],
                                  s = 2,
                                  color = colors_dict['color2'])
            
            # end point
            axs_branch[0].scatter(end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['X'], 
                                  end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Y'],
                                  s = 2,
                                  color = 'r')
            
            # connection line
            axs_branch[0].plot(line_coordinates['X'],
                                line_coordinates['Y'],
                                color = colors_dict['primecolor'],
                                alpha = 0.5)
             
            axs_branch[0].set_xlim([0, 590.76])
            axs_branch[0].set_ylim([590.76, 0])
            axs_branch[0].set_ylabel('Height [µm]')
            axs_branch[0].text(x = 10, y = 10, s = 'XY', ha = 'left', va = 'top')
            
            # branch measurement
            axs_branch[0].text(x = 10, y = 590, 
                               s = f'branch angle [deg] = {terminal_branches_df.at[terminal_path_ID ,"angle_deg"]}\nbranch angle [rad] = {terminal_branches_df.at[terminal_path_ID ,"angle_rad"]}\nbranch length [µm] = {terminal_branches_df.at[terminal_path_ID ,"length"]}\nbranch euclidean dist. [µm] = {terminal_branches_df.at[terminal_path_ID ,"euc_dist"]}\nbranch contraction = {terminal_branches_df.at[terminal_path_ID ,"contraction"]}\n', 
                                ha = 'left', va = 'bottom',
                                size = 6)
               
                
            ## YZ
            # branch
            axs_branch[1].scatter(branch_coor_terminal_to_soma['Z'],
                                  branch_coor_terminal_to_soma['Y'],
                                  s = 0.5,
                                  color = colors_dict['primecolor'])
            
            # soma
            axs_branch[1].scatter(end_coordinates[end_coordinates['path_ID'] == 1]['Z'], 
                                  end_coordinates[end_coordinates['path_ID'] == 1]['Y'],
                                  s = 2,
                                  color = colors_dict['color2'])
            
            # end point
            axs_branch[1].scatter(end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Z'], 
                                  end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Y'],
                                  s = 2,
                                  color = 'r')
             
            axs_branch[1].set_xlim([0, 300])
            axs_branch[1].tick_params(axis = 'y', size = 0)
            axs_branch[1].set_xlabel('Depth [µm]')
            axs_branch[1].text(x = 10, y = 10, s = 'ZY', ha = 'left', va = 'top')
            
            ## XZ
            # branch
            axs_branch[2].scatter(branch_coor_terminal_to_soma['X'],
                                  branch_coor_terminal_to_soma['Z'],
                                  s = 0.5,
                                  color = colors_dict['primecolor'])
            
            # soma
            axs_branch[2].scatter(end_coordinates[end_coordinates['path_ID'] == 1]['X'], 
                                  end_coordinates[end_coordinates['path_ID'] == 1]['Z'],
                                  s = 2,
                                  color = colors_dict['color2'])
            
            # end point
            axs_branch[2].scatter(end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['X'], 
                                  end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Z'],
                                  s = 2,
                                  color = 'r')
             
            axs_branch[2].set_xlim([0, 590.76])
            axs_branch[2].set_ylim([300, 0])
            axs_branch[2].set_ylabel('Depth [µm]')
            axs_branch[2].set_xlabel('Width [µm]')
        
            
            axs_branch[2].text(x = 10, y = 10, s = 'XZ', ha = 'left', va = 'top')
            
            [ax.grid(False) for ax in axs_branch]
            
            plt.show()
            
            fig_branch_dir = join(cell_morph_plots_dir, 'terminal_branches', cell_ID)
            if not exists(fig_branch_dir):
                mkdir(fig_branch_dir)
            
            save_figures(fig_branch, f'{cell_ID}_{terminal_path_ID}-terminal_branch_measurements', fig_branch_dir, darkmode_bool)
    
    
    # # %% 2D histogram of angle and length of terminal branches
    
    # # fig_hist, ax_hist = plt.subplots(subplot_kw={'projection': 'polar'})
    # fig_hist, ax_hist = plt.subplots()
    
    # fig_hist.suptitle(cell_ID)
    
    # # get angles of branches
    # branch_angles_rad = terminal_branches_df["angle_rad"].to_numpy()
    
    # ## start with double the number of desired binsizes to end up with histogram
    # ## that doesnt split at 0
    # # initialize values for histogram
    # n_bins = 8
    # binsize = (2 * np.pi) / n_bins
    # hist_bins = np.arange(0, 2*np.pi + binsize, binsize)
    
    # # get histogram
    # hist_angles_occu, bins_angles = np.histogram(branch_angles_rad, hist_bins)
    
    
    
    # # plot histogram as barplot
    # ax_hist.bar(bins_angles[:-1], hist_angles_occu, width = binsize, align = 'edge')
    
    # ax_hist.set_xticks(bins_angles)
    # ax_hist.set_xticklabels(np.arange(0, 360 + 1, 360 / n_bins, dtype = int), rotation = 45)
    
    # %% absolute polar plot
    # 2D histogram of angle and length of terminal branches
    
    # branch categorization 
    ### histogram ###
    ## start with double the number of desired binsizes to end up with histogram
    ## that doesnt split at 0
    # initialize values for histogram
    resul_n_bins = 8
    resul_binsize = (2 * np.pi) / resul_n_bins
    n_bins = resul_n_bins * 2
    binsize = (2 * np.pi) / n_bins
    hist_bins = np.arange(0, 2*np.pi + binsize, binsize)
    
    # loop through all branches to assign specific color for length
    # skip (drop) path 1, i.e. soma
    for branch_idx in terminal_branches_df.drop(index = 1).index.to_list():
    
        # get angles of branches
        branch_angles_rad = terminal_branches_df.at[branch_idx, "angle_rad"]    
    
        # get histogram
        hist_angles_occu_pre, bins_angles_pre = np.histogram(branch_angles_rad, hist_bins)
    
        # define empty list to populate after combining count in bins
        hist_angles_occu = [None] * resul_n_bins
        bins_angles = [None] * resul_n_bins
        
        # define bin idc to be combined
        idc_bins = np.arange(0, n_bins)
        idc_bins = np.roll(idc_bins, 1)
        
        # loop through resulting bins to combine previouse bins
        for bin_idx in range(resul_n_bins):
            idc_bins_resul = idc_bins[(bin_idx*2):(bin_idx*2)+1+1]
            
            hist_angles_occu[bin_idx] = np.sum(hist_angles_occu_pre[idc_bins_resul])
            bins_angles[bin_idx] = bins_angles_pre[idc_bins[(bin_idx*2)]]
    
        # get bin_id of current branch    
        terminal_branches_df.at[branch_idx, 'bin_id'] = np.argmax(hist_angles_occu)
    
    
    
    # %%
    
    
    # sort dataframe of terminal branches measurements to plot histogram
    terminal_branches_df.sort_values('length', inplace = True)
    
    fig_hist, ax_hist = plt.subplots(subplot_kw={'projection': 'polar'},
                                      layout = 'constrained',
                                      height_ratios= [1],
                                      width_ratios=[1],
                                      figsize = get_figure_size(width = 185.5))
    
    # ax_hist.set_theta_offset(-np.pi / 8)
    
    # fig_hist, ax_hist = plt.subplots(layout = 'constrained')
    
    fig_hist.suptitle(cell_ID)
    

    ### color code for length of branches ###
    # initialise color code
    norm_min = 0
    norm_max = 1000
    cmap_str = 'plasma'
    norm = mtl.colors.Normalize(norm_min, norm_max)
    cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
    
    # colorbar
    fig_hist.colorbar(cmap, ax = ax_hist, label = 'Terminal branch length [µm]')
    
    # define array with number of previouse numbers of branches in bin
    bottom = [0] * resul_n_bins
    
    # loop through all branches to assign specific color for length
    # skip (drop) path 1, i.e. soma
    for branch_idx in terminal_branches_df.drop(index = 1).index.to_list():
    
        # get angles of branches
        branch_angles_rad = terminal_branches_df.at[branch_idx, "angle_rad"]
        branch_length = terminal_branches_df.at[branch_idx, "length"]
        branch_bin = terminal_branches_df.at[branch_idx, "bin_id"].astype(int)
        
        # create empty bins and assign branch to bin
        hist_angles_occu = [0] * resul_n_bins
        hist_angles_occu[branch_bin] = 1
        
        # plot histogram as barplot
        ax_hist.bar(bins_angles, hist_angles_occu, bottom = bottom,
                    width = resul_binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = cmap.to_rgba(branch_length))
            
        # add to bottom list for next step
        bottom = np.add(bottom, hist_angles_occu)
        
    
    # plot terminal branch lines
    for i in terminal_branches_df.index:
        ax_hist.plot([terminal_branches_df.at[i, 'angle_rad'], terminal_branches_df.at[i, 'angle_rad']], [0.5, 2.5], 'w', alpha = 0.25)
    
    # x axis
    ax_hist.set_xticks(np.arange(0, np.pi*2, np.pi / 2))
    ax_hist.set_xticklabels(['p', 'd', 'a', 'v'])
    
    # y axis
    ymax = np.ceil(max(bottom)/5)*5
    
    if ymax <= 15:
        ymax = 15
    
    ax_hist.set_yticks(ticks = np.arange(0, ymax + 1, 5))
    # ax_hist.set_yticks(np.arange(0, round_to_base(max(bottom)+1, 5)+1, 1), minor = True)
    
    
    ax_hist.grid(True, alpha = 0.5)
    
    plt.show()    
    
    # set polar plots directory
    polar_plots_dir = join(cell_morph_plots_dir, 'polar_plots',  'absolute_polar_plots')
    save_figures(fig_hist, f'{cell_ID}-absolute_polar_plot-terminal_branch_orientation-colorcoded_length', polar_plots_dir, darkmode_bool)


# %% normalised polar plot

    max_hist_angles_occu = np.max(bottom)

    fig_hist_norm, ax_hist_norm = plt.subplots(subplot_kw={'projection': 'polar'},
                                                layout = 'constrained',
                                                height_ratios= [1],
                                                width_ratios=[1],
                                                figsize = get_figure_size(width = 185.5))
    # set title
    fig_hist_norm.suptitle(cell_ID)
    
    # colorbar
    fig_hist_norm.colorbar(cmap, ax = ax_hist_norm, label = 'Terminal branch length [µm]')
    
    # define bottom (empty bins for histogram)
    bottom_norm = [0] * resul_n_bins
    norm_bin_height = 1 / max_hist_angles_occu
    
    # loop through all branches to assign specific color for length
    # skip (drop) path 1, i.e. soma
    for branch_idx in terminal_branches_df.drop(index = 1).index.to_list():
    
        # get angles of branches
        branch_angles_rad = terminal_branches_df.at[branch_idx, "angle_rad"]
        branch_length = terminal_branches_df.at[branch_idx, "length"]
        branch_bin = terminal_branches_df.at[branch_idx, "bin_id"].astype(int)
        
        # create empty bins and assign branch to bin
        hist_angles_occu = [0] * resul_n_bins
        hist_angles_occu[branch_bin] = norm_bin_height
        
        # plot histogram as barplot
        ax_hist_norm.bar(bins_angles, hist_angles_occu, bottom = bottom_norm,
                         width = resul_binsize, 
                         align = 'edge',
                         edgecolor = 'none',
                         color = cmap.to_rgba(branch_length))
            
        # add to bottom list for next step
        bottom_norm = np.add(bottom_norm, hist_angles_occu)
        
    # plot terminal branch lines
    for i in terminal_branches_df.index:
        ax_hist_norm.plot([terminal_branches_df.at[i, 'angle_rad'], terminal_branches_df.at[i, 'angle_rad']], [norm_bin_height*0.2, norm_bin_height*0.8], 'w', alpha = 0.25)

    # x axis
    ax_hist_norm.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
    ax_hist_norm.set_xticklabels(['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv'])

    # grid
    ax_hist_norm.grid(True, alpha = 0.5)
    
    # yaxis
    ax_hist_norm.set_ylim([0, 1])
    yticks = np.arange(0, 1 + norm_bin_height, norm_bin_height)
    ax_hist_norm.set_yticks(ticks = [1.0])
    
    plt.show()    
    
    # set polar plots directory
    norm_polar_plots_dir = join(cell_morph_plots_dir, 'polar_plots', 'normalized_polar_plots')
    save_figures(fig_hist_norm, f'{cell_ID}-normalized_polar_plot-terminal_branch_orientation-colorcoded_length', norm_polar_plots_dir, darkmode_bool)


# %% 

    # write histogram to dataframe
    polar_plot_occurrances[cell_ID] = bottom

    # save terminal_branches_df that contains all measurements of terminal branches for each cell
    terminal_branches_df.to_excel(join(cell_morph_descrip_dir, 'terminal_branches_measurements', f'{cell_ID}-terminal_branches.xlsx'), index_label = 'path_ID')   

    
# %% save dataframe

# reset index with orientation angle of bins in rad
polar_plot_occurrances.index = bins_angles
polar_plot_occurrances.to_excel(join(cell_morph_descrip_dir, 'polar_plot_occurrances.xlsx'), index_label = 'orientation_rad')

print('Finished!')

# %% todo list



# necessary functions:
    # coordinate in path
        # output: in_path_bool, in_path_ID, idx_intersect


# TODO:
    # polar plot with length-colorcoded lines
    # polar plot with histogram
        # bin size (quaters or )
    # histogram
    # calc length of terminal branches
        # get intersect points
    
    # get sizes of corresponding stack
        ## append functionality to coordinates exporting script

    




