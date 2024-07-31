# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 09:47:58 2024

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
from cellmorphology.functions_cellmorph import clean_OnPath_column_to_path_ID_n_label, calc_length_of_branch

# path measurements directory
path_measurements_dir = join(cell_morph_descrip_dir, 'path_measurements')


# %%

# get onlyfiles list of coordinates folder
onlyfiles = get_onlyfiles_list(cell_morph_traces_coordinates_dir)

# get cell_IDs from list of all coordinates filenames
cell_IDs = [filename[:5] for filename in [f for f in onlyfiles if 'all_coordinates' in f]]

# get all filenames
filenames_df = pd.DataFrame({'all_coordinates_filename' : [f for f in onlyfiles if 'all_coordinates' in f],
                             'terminal_points_coordinates_filename' : [f for f in onlyfiles if 'terminal_last_coordinates' in f],
                             'primary_points_coordinates_filename' : [f for f in onlyfiles if 'primary_last_coordinates' in f],
                             'soma_coordinates_filename' : [f for f in onlyfiles if 'soma_coordinates' in f]},
                            index = cell_IDs)



# %%

# create dataframe for all cell_IDs
total_cable_length_df = pd.DataFrame(columns=['total_cable_length-neurites', 'total_cable_length-dendrites', 'total_cable_length-axons'], index = cell_IDs)
total_cable_length_df.index.name = 'cell_IDs'


# cell_ID = 'E-056'
for cell_ID in cell_IDs:
    # get soma coordinates
    # soma_coordinates = end_coordinates[end_coordinates['path_ID'] == 1]
    soma_coor_filename = filenames_df.at[cell_ID, 'soma_coordinates_filename']
    soma_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, soma_coor_filename)) 
    clean_OnPath_column_to_path_ID_n_label(soma_coordinates)
    
    ## all coordinates
    # get filename and load csv
    all_coor_filename = filenames_df.at[cell_ID, 'all_coordinates_filename']
    all_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, all_coor_filename))
    
    # clean up dataframe
    clean_OnPath_column_to_path_ID_n_label(all_coordinates)
    
    # remove soma coordinates from all coordinates dataframe
    soma_index_range = soma_coordinates.isin(all_coordinates).index
    all_coordinates.drop(soma_coordinates.isin(all_coordinates).index, inplace = True)
    
    # plot cell with all branches and end points marked
    path_IDs = all_coordinates['path_ID'].drop_duplicates().astype(int).to_list()
    
    # create dataframe for all paths
    path_measurements_df = pd.DataFrame(columns=['path_length', 'path_label'], index = path_IDs)
    path_measurements_df.index.name = 'path_IDs'
    
    # loop through path_IDs
    for path_ID in path_IDs:
        
        # coordinates and label of current path
        path_coordinates = all_coordinates[all_coordinates['path_ID'] == path_ID]
        path_label = all_coordinates[all_coordinates['path_ID'] == path_ID]['path_label'].drop_duplicates().values[0]
    
        # calculate the euclidean distance of current path
        path_length = calc_length_of_branch(path_coordinates.drop(columns = ['path_ID', 'path_label']))
    
        # write to dataframe
        path_measurements_df.at[path_ID, 'path_length'] = path_length
        path_measurements_df.at[path_ID, 'path_label'] = path_label
        
    
    # save dataframe
    path_measurements_df.to_excel(join(path_measurements_dir, f'{cell_ID}-path_measurements.xlsx'), index_label= 'path_ID')
    
    # calc total cable length for cell
    total_cable_length_df.at[cell_ID, 'total_cable_length-neurites'] = path_measurements_df['path_length'].sum()
    total_cable_length_df.at[cell_ID, 'total_cable_length-dendrites'] = path_measurements_df[path_measurements_df['path_label'] == 'dendrite']['path_length'].sum()
    total_cable_length_df.at[cell_ID, 'total_cable_length-axons'] = path_measurements_df[path_measurements_df['path_label'] == 'axon']['path_length'].sum()
    
    # feedback
    print(f'{cell_ID} done')
    
# save dataframe
total_cable_length_df.to_excel(join(cell_morph_descrip_dir, 'total_cable_length.xlsx'), index_label= 'path_ID')

# feedback
print('Finished!')





