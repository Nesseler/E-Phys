# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:15:07 2024

@author: nesseler
"""

# initialize packages
import numpy as np
import math


# %% data import

def clean_OnPath_column_to_path_ID_n_label(coordinates_dataframe):
    
    onPath_column = coordinates_dataframe['On Path']

    for i, txt in enumerate(onPath_column):
        path_ID = [int(s.replace("-", "").strip("()")) for s in txt.split() if s.replace("-", "").strip("()").isdigit()]
        
        if len(path_ID) > 1:
            raise ValueError('More than one integer in path_ID')
        elif len(path_ID) < 1:
            raise ValueError('Less than one integer in path_ID')
        elif len(path_ID) == 1:
            path_ID = path_ID[-1]
        
        
        # txt_ls = txt.split()
            
        if 'axon' in txt:
            path_label = 'axon'
        elif 'soma' in txt:
            path_label = 'soma'
        else:
            path_label = 'dendrite'
            
        if path_ID == 1:
            path_label = 'soma'
    
        coordinates_dataframe.at[i, 'path_ID'] = path_ID
        coordinates_dataframe.at[i, 'path_label'] = path_label
        
    coordinates_dataframe.drop(columns = ['On Path'], inplace = True)

    return coordinates_dataframe


# %% caculations

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
      

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


# %% branch reconstruction


def find_parent_path(path_ID, path_IDs_to_search, allcoor_paths_dict):
    
    # get path coordinates
    path_coor = allcoor_paths_dict[path_ID]

    # get first node
    firstnode = tuple(path_coor.iloc[0][['X', 'Y', 'Z']])
    
    # set parent path Id and initial index for the potential parent
    parent_path_ID = None
    pot_parent_idx = 0
    
    # loop through path_IDs_to_search
    while not parent_path_ID:
    
        pot_parent_path_ID = path_IDs_to_search[pot_parent_idx]
             
        # skip branch itself 
        if pot_parent_path_ID != path_ID:
    
            # coordinates of potential parent path
            pot_parent_path_coor = allcoor_paths_dict[pot_parent_path_ID]

            # test if node is in potential parent path 
            # create mask where elements are True that correspond to the first node coordinates
            coor_mask = (pot_parent_path_coor[['X', 'Y', 'Z']] == firstnode).all(axis = 1)
    
            if coor_mask.any():
                
                # get the first True index
                intersect_index = coor_mask.idxmax()
                
                # get intersection coordinates
                intersect_coor = pot_parent_path_coor.loc[intersect_index]
   
                # get all coordinates until intersection point
                ## get all indices of parent path until intersection index
                parent_indices = [i for i in pot_parent_path_coor.index.to_list() if i <= intersect_index] 
                
                # avoid finding the root of another bifurcation
                if len(parent_indices) > 1 or intersect_coor['path_ID'].astype(int) == 1:
                    parent_path_ID = intersect_coor['path_ID'].astype(int)
                    
        pot_parent_idx += 1           
                                   
    return parent_path_ID, intersect_index



# %% polar plot

def calc_polar_histo_binangles(n_bins = 8):
    '''
    Function calculates border angles of bins in rad for the polar plot histogram.
    Parameters: 
        n_bins (int): Number of bins in polar histogram. Default is 8.
    Returns:
        bin_angles (list): List of border angles in rad.
    '''
    
    # step size is set bin number of resulting bins and 2*pi
    bin_stepsize = (2 * np.pi) / n_bins
    
    # start point: half of step size
    # because of rotated polar bins
    bin_start = bin_stepsize / 2
    
    # calc bin borders
    bin_angles = np.arange(bin_start, np.pi * 2, bin_stepsize)

    # roll bin borders
    bin_angles = np.roll(bin_angles, 1)    

    return bin_angles, bin_stepsize



def plot_colorcoded_polar_normed(polar_occurances_df, max_n_neurites, ax, n_bins, cmap):
    """
    Function uses terminal_brances dataframe to plot the occurrances of branches color-coded in 
    polar plots.
    """
    
    # get bins angles and binsize
    bins_angles, binsize = calc_polar_histo_binangles(n_bins)
    
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
    




# %% legacy functions

# def node_in_path(first_node_coor, pot_parent_path_coor):
#     # create mask where elements are True that correspond to the first node corrdinates
#     coor_mask = pot_parent_path_coor == first_node_coor
    
#     # get intersection point
#     # where all coordinates are the same
#     intersect_mask = coor_mask.query('X == True & Y == True & Z == True')
    
#     # test if parent, i.e.: mask does or doesn't contain True values
#     if intersect_mask.empty:
#         in_path_bool = False
        
#     else:
#         in_path_bool = True
        
#     return in_path_bool


    
    