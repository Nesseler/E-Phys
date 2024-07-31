# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 16:48:58 2024

@author: nesseler
"""

import pandas as pd
from os import mkdir
from os.path import join, exists
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np


from parameters.directories_win import table_file, cell_morph_descrip_dir, cell_morph_traces_coordinates_dir, cell_morph_plots_dir

from functions.functions_import import get_onlyfiles_list

# create directory path
terminal_branches_measurements_dir = join(cell_morph_descrip_dir, 'terminal_branches_measurements')

# get all file names from directory
terminal_branches_measurements_onlyfiles = get_onlyfiles_list(terminal_branches_measurements_dir)

# get cell_IDs from list of files
cell_IDs = [file_str[:5] for file_str in terminal_branches_measurements_onlyfiles]

# create dataframe from files list and cell_IDs
files_df = pd.DataFrame(terminal_branches_measurements_onlyfiles, columns= ['filename'],  index = cell_IDs)
files_df.index.name = 'cell_ID'

# tests for loop
# cell_ID = cell_IDs[0]

# dataframe to export to
n_terminal_points_df = pd.DataFrame(columns = ['neuritic_terminals', 'dendritic_terminals', 'axonic_terminals'],
                                    index = cell_IDs)
n_terminal_points_df.index.name = 'cell_ID'

# loop through all cell_IDs
for cell_ID in cell_IDs:
    
    # get filename
    filename = files_df.at[cell_ID, 'filename']

    # load terminal branches dataframe
    terminal_branches_df = pd.read_excel(join(terminal_branches_measurements_dir, filename))
    
    # get number of neuritic, dendritic, and axonic terminal points
    n_neuritic = terminal_branches_df.shape[0]-1
    n_dendritic = terminal_branches_df[terminal_branches_df['path_label'] == 'dendrite'].shape[0]
    n_axonic = terminal_branches_df[terminal_branches_df['path_label'] == 'axon'].shape[0]
    
    # quality check
    if n_dendritic + n_axonic != n_neuritic:
        raise ValueError(f'Number of dendritic and axonic terminal points does not match! {cell_ID}')

    # write to dataframe
    n_terminal_points_df.loc[cell_ID] = {'neuritic_terminals' : n_neuritic,
                                         'dendritic_terminals' : n_dendritic,
                                         'axonic_terminals' : n_axonic}

# save dataframe to directory
n_terminal_points_df.to_excel(join(cell_morph_descrip_dir, 'n_terminal_points.xlsx'), 
                              index_label= 'cell_ID')


# %% primary points

from cellmorphology.functions_cellmorph import clean_OnPath_column_to_path_ID_n_label

# import end points coordinates of primary paths
# get filenames of primary_last_coordinates
primary_last_coordinates_onlyfiles = [filename for filename in get_onlyfiles_list(cell_morph_traces_coordinates_dir) if 'primary_last_coordinates' in filename]

# create dataframe with filenames
primary_filenames_df = pd.DataFrame({'filename' : primary_last_coordinates_onlyfiles},
                                    index = [filename[:5] for filename in primary_last_coordinates_onlyfiles])
primary_filenames_df.index.name = 'cell_ID'

# dataframe to export to
n_primary_points_df = pd.DataFrame(columns = ['neuritic_primaries', 'dendritic_primaries', 'axonic_primaries'],
                                    index = cell_IDs)
n_primary_points_df.index.name = 'cell_ID'

# loop through cell_IDs
for cell_ID in cell_IDs:

    # get path to coordinates file
    primary_coordinates_path = join(cell_morph_traces_coordinates_dir, primary_filenames_df.at[cell_ID, 'filename'])
    
    # load coordinates file
    primary_last_coordinates = pd.read_csv(primary_coordinates_path)
    
    # clean coordinates dataframe
    clean_OnPath_column_to_path_ID_n_label(primary_last_coordinates)
    
    # get length of coordinates file as number of primary paths
    n_primary_neuritic = primary_last_coordinates.shape[0]
    n_primary_dendritic = primary_last_coordinates[primary_last_coordinates['path_label'] == 'dendrite'].shape[0]
    n_primary_axonic = primary_last_coordinates[primary_last_coordinates['path_label'] == 'axon'].shape[0]
    
    # quality check
    if n_primary_dendritic + n_primary_axonic != n_primary_neuritic:
        raise ValueError(f'Number of dendritic and axonic primary points does not match! {cell_ID}')
        
    # write to dataframe
    n_primary_points_df.loc[cell_ID] = {'neuritic_primaries' : n_primary_neuritic,
                                        'dendritic_primaries' : n_primary_dendritic,
                                        'axonic_primaries' : n_primary_axonic}
    
# save dataframe to directory
n_primary_points_df.to_excel(join(cell_morph_descrip_dir, 'n_primary_points.xlsx'), 
                             index_label= 'cell_ID')


# %% bifurcation points

# calculate number of bifurcations as difference in number of terminal and primary points 
n_bifurcation_points_df = pd.DataFrame({'n_bifurations_neuritic' : n_terminal_points_df['neuritic_terminals'] - n_primary_points_df['neuritic_primaries'],
                                        'n_bifurations_dendritic' : n_terminal_points_df['dendritic_terminals'] - n_primary_points_df['dendritic_primaries'],
                                        'n_bifurations_axonic' : n_terminal_points_df['axonic_terminals'] - n_primary_points_df['axonic_primaries']},
                                       index= cell_IDs)

# save dataframe to directory
n_bifurcation_points_df.to_excel(join(cell_morph_descrip_dir, 'n_bifurcation_points.xlsx'), 
                                 index_label= 'cell_ID')


# %% bifurcation ratio
# calculate bifurcation ratio as ratio between number of terminal and primary points

# create dataframe to export to
bifurcation_ratio_df = pd.DataFrame(columns = ['bifurcation_ratio_neuritic', 'bifurcation_ratio_dendritic', 'bifurcation_ratio_axonic'],
                                    index = cell_IDs)
bifurcation_ratio_df.index.name = 'cell_ID'

# loop through cell_IDs to avoid division by zero
for cell_ID in cell_IDs:
    
    # get terminal points
    n_neuritic_terminals = n_terminal_points_df.at[cell_ID, 'neuritic_terminals']
    n_dendritic_terminals = n_terminal_points_df.at[cell_ID, 'dendritic_terminals']
    n_axonic_terminals = n_terminal_points_df.at[cell_ID, 'axonic_terminals']
    
    # get primary points
    n_neuritic_primaries = n_primary_points_df.at[cell_ID, 'neuritic_primaries']
    n_dendritic_primaries = n_primary_points_df.at[cell_ID, 'dendritic_primaries']
    n_axonic_primaries = n_primary_points_df.at[cell_ID, 'axonic_primaries']
    
    # neuritic bif ratio
    bif_ratio_neuritic = n_neuritic_terminals / n_neuritic_primaries
    
    # dendritic bif ratio
    bif_ratio_dendritic = n_dendritic_terminals / n_dendritic_primaries
    
    # axonic bif ratio
    if n_axonic_primaries == 0 and n_axonic_terminals == 0:
        bif_ratio_axonic = np.nan
    elif n_axonic_primaries == 0 and n_axonic_terminals > 0:
        bif_ratio_axonic = n_axonic_terminals / 1
    else:
        bif_ratio_axonic = n_axonic_terminals / n_axonic_primaries

    # write to dataframe
    bifurcation_ratio_df.loc[cell_ID] = {'bifurcation_ratio_neuritic' : bif_ratio_neuritic,
                                         'bifurcation_ratio_dendritic' : bif_ratio_dendritic,
                                         'bifurcation_ratio_axonic' : bif_ratio_axonic}

# save dataframe to directory
bifurcation_ratio_df.to_excel(join(cell_morph_descrip_dir, 'bifurcation_ratios.xlsx'),
                              index_label= 'cell_ID')















































