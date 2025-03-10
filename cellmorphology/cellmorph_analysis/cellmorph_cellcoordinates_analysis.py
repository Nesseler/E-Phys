# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 10:14:45 2025

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# get metrics saving directory
from cellmorphology.cellmorph_functions.cellmorph_dir import cell_morph_traces_coordinates_dir as cellmorph_coordinates_dir

# get files in traces folder
from functions.get_onlyfiles_list import get_onlyfiles_list
onlyfiles = get_onlyfiles_list(cellmorph_coordinates_dir)

# get cell_IDs from list of all coordinates filenames
cell_IDs = [filename[:5] for filename in [f for f in onlyfiles if 'all_coordinates' in f]]

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# set neurite types
neurite_types = ['neurites', 
                 'dendrites',
                 'axons']

vplots = True


# %% check wether cell has been anaylzed

# cellIDs_toAnalyze = ['E-292']
cellIDs_toAnalyze = ['E-137']

# %% define output

height_width_depth = pd.DataFrame(columns = [ntype + '-' + col for ntype in neurite_types for col in ['x_min', 'y_min', 'z_min', 'width', 'height', 'depth']], 
                                  index = cellIDs_toAnalyze)

total_cable_length = pd.DataFrame(columns = [ntype + '-total_cable_length' for ntype in neurite_types], 
                                  index = cellIDs_toAnalyze)

n_stems = pd.DataFrame(columns = [f'n_stems-{ntype}' for ntype in neurite_types],
                       index = cellIDs_toAnalyze)

bifurcation_ratios = pd.DataFrame(columns = [f'bifurcation_ratio-{ntype}' for ntype in neurite_types],
                                  index = cellIDs_toAnalyze)

# rename index
for df in [height_width_depth, total_cable_length, n_stems, bifurcation_ratios]:
    df.index.name = 'cell_ID'


# %% start analysis (per cell)

for cell_ID in cellIDs_toAnalyze:
    
    # feedback
    print(f'Started with {cell_ID}')
    
    # get region of cell
    region = MetaData.at[cell_ID, 'Region']


    # %% load coordinates
    
    # all coordinates dict for one cell
    cell_coordinates = dict.fromkeys(['all_coor', 'end_coor', 'soma_coor', 'path_IDs'])

    # import function that cleans up imported coordinates table
    from cellmorphology.cellmorph_functions.cellmorph_functions import clean_OnPath_column_to_path_ID_n_label
     
    ### coordinates
    # all coordinates
    all_coordinates_path = join(cellmorph_coordinates_dir, f'{cell_ID}-all_coordinates.csv')
    cell_allcoordinates = pd.read_csv(all_coordinates_path)
    clean_OnPath_column_to_path_ID_n_label(cell_allcoordinates)
    
    # primary last coordinates
    primary_coordinates_path = join(cellmorph_coordinates_dir, f'{cell_ID}-primary_last_coordinates.csv')
    cell_primarycoordinates = pd.read_csv(primary_coordinates_path)
    clean_OnPath_column_to_path_ID_n_label(cell_primarycoordinates) 
    
    # end / last / terminal coordinates
    last_coordinates_path = join(cellmorph_coordinates_dir, f'{cell_ID}-terminal_last_coordinates.csv')
    cell_endcoordinates = pd.read_csv(last_coordinates_path)
    clean_OnPath_column_to_path_ID_n_label(cell_endcoordinates) 
    
    # get soma coordinates
    soma_coordinates = cell_allcoordinates[cell_allcoordinates['path_ID'] == 1]
     
    # get all path_IDs
    path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()
    
    # set dict for all coordinates
    cell_coordinates = {'all_coor' : cell_allcoordinates,
                        'pri_coor' : cell_primarycoordinates,
                        'end_coor' : cell_endcoordinates,
                        'soma_coor': soma_coordinates,
                        'path_IDs' : path_IDs}
   
    # %% plot coordinates

    if vplots:   
        # load plotting function    
        from cellmorphology.cellmorph_analysis.plot_cellcoordinates_analysis import plot_cellcoordinates
        
        # plot 
        plot_cellcoordinates(cell_ID = cell_ID, 
                             cell_coordinates = cell_coordinates,
                             region = region)
        
    # %% height, width, depth
    
    
    # %% plot height, width, depth


    # %% total_cable_length
    
    
    # %% reconstruct terminal branches
        # from end point (terminal point) to soma
        
        
    # %% number of primary and terminal points
    
    
    # %% number of stems and bifuracation ratio
    
    
    # %% plot terminal branches orientation
    
    
    # %% get axon carrying dendrite
    
    