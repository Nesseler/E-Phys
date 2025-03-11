# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 16:53:53 2025

@author: nesseler
"""

import numpy as np
import pandas as pd

# %% polar plots parameters

# set orientation labels
orientation_labels = ['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv']

# set number of bins
n_bins = 8
binsize = (2 * np.pi) / n_bins
bins = np.arange(binsize / 2, 2 * np.pi, binsize)

# roll bins to start with posterior bin
hist_bins = np.roll(bins, 1)

# create temporary bins for assignment
temp_nbins = n_bins * 2
temp_binsize = (2 * np.pi) / temp_nbins
temp_bins = np.arange(0, (2 * np.pi) + temp_binsize, temp_binsize)


# %% polar plots functions

def assign_bin2branches(branch_angles_rad):
    '''
    This function assigns a given terminal branch angle to the corresponding
    bin in the polar plot.
    Parameters:
        branch_angles_rad: pandas Series, describes to angle of each terminal
                           branch in rad.
    Return:
        bin_Ids: pandas Series (like input), describes to assiged bins
    '''
    
    # check for type
    if type(branch_angles_rad) != pd.Series:
        raise ValueError('Did not recieve pandas Series!')
    
    # create Series like input
    bin_Ids = pd.Series(index = branch_angles_rad.index, name = 'bin_id')
    
    # iterate through branches
    for branch_id, angle in branch_angles_rad.items():
    
        # get assignment to temp bins
        temp_assign, _ = np.histogram([angle], temp_bins)
        
        # roll and combine
        assign = np.roll(temp_assign, 1).reshape((n_bins, 2)).sum(axis = 1)
        
        # get index
        assigned_bin = assign.argmax()
        
        # write to Series
        bin_Ids[branch_id] = assigned_bin

    return bin_Ids


def assign_bin2branch(angle_rad):
    '''
    This function assigns a given terminal branch angle to the corresponding
    bin in the polar plot.
    Parameters:
        angle_rad: float, describes to angle of on terminal branch in rad.
    Return:
        assigned_bin: int, describes to assiged bin
    '''
    
    # check for type
    if type(angle_rad) != float:
        raise ValueError('Did not recieve float!')
    
    # get assignment to temp bins
    temp_assign, _ = np.histogram([angle_rad], temp_bins)
    
    # roll and combine
    assign = np.roll(temp_assign, 1).reshape((n_bins, 2)).sum(axis = 1)
    
    # get index
    assigned_bin = int(assign.argmax())
        
    return assigned_bin


def get_orientation_occurances(angles_rad):
    '''
    This function calculates the orientation occurances per bin for a list of
    orientations in rad.
    Parameters:
        angles_rad: list of floats, describes the orientation angles for a 
                    series of terminal branches.
    Returns:
        occu: list, occurances of angles per bins
    '''
    # get occurances per bin for temporary bins
    temp_occu, _ = np.histogram(angles_rad, temp_bins)
    
    # get occurances per bin for original bins
    occu = np.roll(temp_occu, 1).reshape((n_bins, 2)).sum(axis = 1)

    return occu


# %% polar plots plotting


def create_polar_histogram(ax, terminal_branches):

    # set colors 
    from cellmorphology.cellmorph_functions.cellmorph_init_plotting import neurite_color_dict
    neurite_color_dict = neurite_color_dict['all']
             
    # define array with number of previouse numbers of branches in bin
    bottom = [0] * n_bins
    
    # get terminal branch Ids
    if 1 in terminal_branches.index:
        terminalbranch_ids = terminal_branches.index.drop(index = 1).to_list()
    else:
        terminalbranch_ids = terminal_branches.index.to_list()
    
    # loop through all branches to assign specific color for length
    for branch_idx in terminalbranch_ids:
    
        # get angles of branches
        branch_label = terminal_branches.at[branch_idx, "path_label"]
        branch_angles_rad = terminal_branches.at[branch_idx, "angle_rad"]
        branch_bin = terminal_branches.at[branch_idx, "bin_id"]
        
        # create empty bins and assign branch to bin
        hist_angles_occu = [0] * n_bins
        hist_angles_occu[branch_bin] = 1
        
        # plot histogram as barplot
        ax.bar(hist_bins, hist_angles_occu, bottom = bottom,
               width = binsize, 
               align = 'edge',
               edgecolor = 'none',
               color = neurite_color_dict[branch_label],
               label = branch_label)
            
        # add to bottom list for next step
        bottom = np.add(bottom, hist_angles_occu)
        
    # x axis
    ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
    ax.set_xticklabels(['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv'])
    
    # y axis
    ymax = np.ceil(max(bottom)/5)*5

    # if ymax <= 15:
    #     ymax = 15

    ax.set_yticks(ticks = np.arange(0, ymax + 1, 5))
    
    return bottom


