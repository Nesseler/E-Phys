# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 13:53:45 2024

@author: nesseler
"""

import pandas as pd
from os.path import join
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import math

from parameters.directories_win import cell_morph_traces_coordinates_dir, cell_morph_plots_dir, table_file
from getter.get_onlyfiles_list import get_onlyfiles_list

from functions.functions_plotting import get_colors, save_figures
from functions.functions_useful import round_to_base

# get onlyfiles list
onlyfiles = get_onlyfiles_list(cell_morph_traces_coordinates_dir)

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# extract only filenames with all coordinates
onlyfiles_allcoor = [f for f in onlyfiles if 'all_coordinates' in f]
onlyfiles_endcoor = [f for f in onlyfiles if 'last_coordinates' in f]

# test 
cell_in_files = 1
all_coor_filename = onlyfiles_allcoor[cell_in_files]

# get cell ID from filename
cell_ID = all_coor_filename[:5]

# read all coordinates file
all_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, all_coor_filename)) 
end_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, onlyfiles_endcoor[cell_in_files])) 

# clean up dataframes
def clean_coordinates_df(coordinates_df):
    coordinates_df['path_ID'] = [int(s.replace("Path ", "").replace("(", "").replace(")", "").replace("-", "").replace('1 [Single Point]', '1')) for s in coordinates_df['On Path']]
    coordinates_df.drop(columns = ['On Path'], inplace = True)

clean_coordinates_df(all_coordinates)

# end_coordinates['On Path'] = end_coordinates['On Path'].replace('Path (1) [Single Point]', '1')
clean_coordinates_df(end_coordinates)

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

fig_cell, ax_cell = plt.subplots()

fig_cell.suptitle(cell_ID)
ax_cell.scatter(all_coordinates['X'], all_coordinates['Y'], s = 0.5, color = colors_dict['primecolor'])
ax_cell.scatter(end_coordinates['X'], end_coordinates['Y'], s = 0.75, color = 'r')
ax_cell.scatter(end_coordinates['X'][0], end_coordinates['Y'][0], s = 2, color = colors_dict['color2'])
ax_cell.invert_yaxis()
ax_cell.set_box_aspect(1)
ax_cell.set_xlim([0, 590.76])
ax_cell.set_xlabel('Position [µm]')
ax_cell.set_ylim([590.76, 0])
ax_cell.set_ylabel('Position [µm]')


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
    vector_coordinates.drop(columns = ['path_ID'], inplace = True)
    
    ax_cell.plot(vector_coordinates['X'], vector_coordinates['Y'], colors_dict['primecolor'], alpha = 0.5)
    
    xy_vector = vector_coordinates.drop(columns = ['Z'])
    
    
    terminal_branch_df = pd.DataFrame({'x_end': end_point_coordinates.at[idx_end_point,'X'],
                                       'y_end': end_point_coordinates.at[idx_end_point,'Y'],
                                       'z_end': end_point_coordinates.at[idx_end_point,'Z'],
                                       'x_diff': end_point_coordinates.at[idx_end_point, 'X'] - soma_coordinates.at[0, 'X'],
                                       'y_diff': end_point_coordinates.at[idx_end_point, 'Y'] - soma_coordinates.at[0, 'Y'],
                                       'z_diff': end_point_coordinates.at[idx_end_point, 'Z'] - soma_coordinates.at[0, 'Z']},
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
    
    

plt.show()


# %% rough polar plot

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

fig.suptitle(cell_ID)

for i in terminal_branches_df.index:
    
    ax.plot([terminal_branches_df.at[i, 'angle_rad'], terminal_branches_df.at[i, 'angle_rad']], [0.1, 1], 'w', alpha = 0.5)
    ax.set_rticks([1])
    ax.set_xticks(np.arange(0, np.pi*2, np.pi/4))
    
ax.grid(False)

plt.show()


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
        intersect_index = intersect_mask.index.item()
        
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
    terminal_branch_length = calc_length_of_branch(branch_coor_terminal_to_soma.drop(columns = ['path_ID']))
    
    terminal_branches_df.at[terminal_path_ID ,'length'] = terminal_branch_length
       
    ### calculate euclidean distance between origin and terminal of branch
    ## can be calculated with same function but just two coordinates
    line_coordinates = branch_coor_terminal_to_soma.iloc[[0, -1]]
    
    terminal_branch_euc = calc_length_of_branch(line_coordinates.drop(columns = ['path_ID']))
    
    # write to dataframe
    terminal_branches_df.at[terminal_path_ID ,'euc_dist'] = terminal_branch_euc
    
    ### calculate branch contraction
    terminal_branch_contraction = terminal_branch_euc / terminal_branch_length
    
    # write to dataframe
    terminal_branches_df.at[terminal_path_ID ,'contraction'] = terminal_branch_contraction
    
    
#     ### branch verification plot ###
#     ratio = 300 / 590.76
    
#     fig_branch, axs_branch = plt.subplots(nrows = 2,
#                                           ncols = 2,
#                                           height_ratios= [1, ratio],
#                                           width_ratios= [1, ratio],
#                                           sharey = 'row',
#                                           sharex = 'col',
#                                           figsize = [5, 5])
    
    
    
#     plt.subplots_adjust(wspace=0, hspace=0)
    
    
#     axs_branch = axs_branch.flatten()
#     fig_branch.delaxes(axs_branch[-1])
    
#     fig_branch.suptitle(f'{cell_ID} terminal_path_ID: {terminal_path_ID}')
    
#     ## XY
#     # branch
#     axs_branch[0].scatter(branch_coor_terminal_to_soma['X'],
#                           branch_coor_terminal_to_soma['Y'],
#                           s = 0.5,
#                           color = colors_dict['primecolor'])
    
#     # soma
#     axs_branch[0].scatter(end_coordinates[end_coordinates['path_ID'] == 1]['X'], 
#                           end_coordinates[end_coordinates['path_ID'] == 1]['Y'],
#                           s = 2,
#                           color = colors_dict['color2'])
    
#     # end point
#     axs_branch[0].scatter(end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['X'], 
#                           end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Y'],
#                           s = 2,
#                           color = 'r')
    
#     # connection line
#     axs_branch[0].plot(line_coordinates['X'],
#                        line_coordinates['Y'],
#                        color = colors_dict['primecolor'],
#                        alpha = 0.5)
     
#     axs_branch[0].set_xlim([0, 590.76])
#     axs_branch[0].set_ylim([590.76, 0])
#     axs_branch[0].set_ylabel('Position [µm]')
#     axs_branch[0].text(x = 10, y = 10, s = 'XY', ha = 'left', va = 'top')
    
#     # branch measurement
#     axs_branch[0].text(x = 10, y = 590, 
#                        s = f'branch angle [deg] = {terminal_branches_df.at[terminal_path_ID ,"angle_deg"]}\n\
# branch angle [rad] = {terminal_branches_df.at[terminal_path_ID ,"angle_rad"]}\n\
# branch length [µm] = {terminal_branches_df.at[terminal_path_ID ,"length"]}\n\
# branch euclidean dist. [µm] = {terminal_branches_df.at[terminal_path_ID ,"euc_dist"]}\n\
# branch contraction = {terminal_branches_df.at[terminal_path_ID ,"contraction"]}\n\
# ', 
#                        ha = 'left', va = 'bottom',
#                        size = 6)
       
        
#     ## YZ
#     # branch
#     axs_branch[1].scatter(branch_coor_terminal_to_soma['Z'],
#                           branch_coor_terminal_to_soma['Y'],
#                           s = 0.5,
#                           color = colors_dict['primecolor'])
    
#     # soma
#     axs_branch[1].scatter(end_coordinates[end_coordinates['path_ID'] == 1]['Z'], 
#                           end_coordinates[end_coordinates['path_ID'] == 1]['Y'],
#                           s = 2,
#                           color = colors_dict['color2'])
    
#     # end point
#     axs_branch[1].scatter(end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Z'], 
#                           end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Y'],
#                           s = 2,
#                           color = 'r')
     
#     axs_branch[1].set_xlim([0, 300])
#     axs_branch[1].tick_params(axis = 'y', size = 0)
    
#     axs_branch[1].text(x = 10, y = 10, s = 'ZY', ha = 'left', va = 'top')
    
#     ## XZ
#     # branch
#     axs_branch[2].scatter(branch_coor_terminal_to_soma['X'],
#                           branch_coor_terminal_to_soma['Z'],
#                           s = 0.5,
#                           color = colors_dict['primecolor'])
    
#     # soma
#     axs_branch[2].scatter(end_coordinates[end_coordinates['path_ID'] == 1]['X'], 
#                           end_coordinates[end_coordinates['path_ID'] == 1]['Z'],
#                           s = 2,
#                           color = colors_dict['color2'])
    
#     # end point
#     axs_branch[2].scatter(end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['X'], 
#                           end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Z'],
#                           s = 2,
#                           color = 'r')
     
#     axs_branch[2].set_ylim([300, 0])
#     axs_branch[2].set_xlabel('Position [µm]')
    
#     axs_branch[2].text(x = 10, y = 10, s = 'XZ', ha = 'left', va = 'top')
    
#     [ax.grid(False) for ax in axs_branch]
    
#     plt.show()


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

# %% 2D histogram of angle and length of terminal branches

# sort dataframe of terminal branches measurements to plot histogram
terminal_branches_df.sort_values('length', inplace = True)

fig_hist, ax_hist = plt.subplots(subplot_kw={'projection': 'polar'},
                                  layout = 'constrained',
                                  height_ratios= [1],
                                  width_ratios=[1])

# ax_hist.set_theta_offset(-np.pi / 8)

# fig_hist, ax_hist = plt.subplots(layout = 'constrained')

fig_hist.suptitle(cell_ID)

### histogram ###
## start with double the number of desired binsizes to end up with histogram
## that doesnt split at 0
# initialize values for histogram
resul_n_bins = 8
resul_binsize = (2 * np.pi) / resul_n_bins
n_bins = resul_n_bins * 2
binsize = (2 * np.pi) / n_bins
hist_bins = np.arange(0, 2*np.pi + binsize, binsize)

### color code for length of branches ###
# initialise color code
norm_min = 0
norm_max = 1000
cmap_str = 'plasma'
norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# colorbar
fig_hist.colorbar(cmap, ax = ax_hist, label = 'Length [µm]')

# define array with number of previouse numbers of branches in bin
bottom = [0] * n_bins

# loop through all branches to assign specific color for length
for branch_idx in terminal_branches_df.index.to_list():

    # get angles of branches
    branch_angles_rad = terminal_branches_df.at[branch_idx, "angle_rad"]
    branch_length = terminal_branches_df.at[branch_idx, "length"]
    
    # get histogram
    hist_angles_occu_pre, bins_angles_pre = np.histogram(branch_angles_rad, hist_bins)
    
    hist_angles_occu = [None] * resul_n_bins
    bins_angles = [None] * resul_n_bins
    idc_bins = np.arange(0, n_bins)
    idc_bins = np.roll(idc_bins, 1)
    
    for bin_idx in range(resul_n_bins):
        idc_bins_resul = idc_bins[(bin_idx*2):(bin_idx*2)+1+1]
        
        hist_angles_occu[bin_idx] = np.sum(hist_angles_occu_pre[idc_bins_resul])
        bins_angles[bin_idx] = bins_angles_pre[idc_bins[(bin_idx*2)+1]]
        
        print(idc_bins_resul, hist_angles_occu[bin_idx],  hist_angles_occu_pre[idc_bins_resul])
    
    bins_angles = list(bins_angles) + [bins_angles[-1] + resul_binsize]
    
    break

    # plot histogram as barplot
    ax_hist.bar(bins_angles_pre[:-1], hist_angles_occu_pre, bottom = bottom,
                width = binsize, 
                align = 'edge', 
                color = cmap.to_rgba(branch_length), 
                edgecolor = 'none')
    
    # add to bottom list for next step
    # bottom = np.add(bottom, hist_angles_occu_pre)

# ax_hist.set_xticks(np.arange(0, np.pi * 2, np.pi / 2))

# ax_hist.set_xticklabels(['p', 'd', 'a', 'v'], rotation = 45)


for i in terminal_branches_df.index:
    
    ax_hist.plot([terminal_branches_df.at[i, 'angle_rad'], terminal_branches_df.at[i, 'angle_rad']], [0.1, 5], 'w', alpha = 0.5)


# ax_hist.set_xticklabels(np.arange(360 / (8*2), 360 + (360 / (8*2)), 360 / resul_n_bins, dtype = int), rotation = 45)

# create overly complicated list of labels
# ls = []
# for l in np.arange(0, round_to_base(max(bottom)+1, 5)+1, 1):
#     if l % 5 == 0:
#         ls.append(l)
#     else:
#         ls.append(None)

# y axis
ax_hist.set_yticks(ticks = np.arange(0, round_to_base(max(bottom)+1, 5)+1, 5))
# ax_hist.set_yticks(np.arange(0, round_to_base(max(bottom)+1, 5)+1, 1), minor = True)
# ax_hist.set_ylabel('N terminal branches')

ax_hist.grid(True, alpha = 0.5)

# ax_hist.set_axisbelow(True)

# set polar plots directory
polar_plots_dir = join(cell_morph_plots_dir, 'polar_plots')

save_figures(fig_hist, f'{cell_ID}-polar_plot-terminal_branch_orientation-colorcoded_length', polar_plots_dir, darkmode_bool)


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

    




