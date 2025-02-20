# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 14:57:41 2025

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# get metrics saving directory
from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_metrics_dir

# get MetaData
from cellmorphology.AMC_analysis.AMC_analysis_import import get_cells_list
MetaData = get_cells_list()

# get cell_IDs to be analyzed
cell_IDs = MetaData.query('coordinates == "Yes" & paths_checked == "Yes"').index.to_list()
cell_IDs = ['Exp-162']

# set neurite types
neurite_types = ['neurites', 
                 'dendrites', 
                 'glomerular_dendrites', 
                 'nonglomerular_dendrites', 
                 'lateral_dendrites', 
                 'LOTxing_dendrites', 
                 'axons']

vplots = True

# %% load coordinates

print('cell coordinates ...')

# coordinates dict for all cells
coordinates_dict = dict.fromkeys(cell_IDs)

# all coordinates dict for one cell
cell_coordintes_dict = dict.fromkeys(['all_coor', 'end_coor', 'soma_coor', 'path_IDs'])

# get directory of cell coordinates
from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_coordinates_dir

# import function that cleans up imported coordinates table
from cellmorphology.AMC_analysis.AMC_functions import clean_OnPath_column_to_path_ID_n_label_AMCs


for cell_ID in tqdm(cell_IDs):       
    ### coordinates
    # all coordinates
    all_coordinates_path = join(AMCs_coordinates_dir, f'{cell_ID}-all_coordinates.csv')
    cell_allcoordinates = pd.read_csv(all_coordinates_path)
    clean_OnPath_column_to_path_ID_n_label_AMCs(cell_allcoordinates)
    
    # primary last coordinates
    primary_coordinates_path = join(AMCs_coordinates_dir, f'{cell_ID}-primary_last_coordinates.csv')
    cell_primarycoordinates = pd.read_csv(primary_coordinates_path)
    clean_OnPath_column_to_path_ID_n_label_AMCs(cell_primarycoordinates) 
    
    # end / last / terminal coordinates
    last_coordinates_path = join(AMCs_coordinates_dir, f'{cell_ID}-terminal_last_coordinates.csv')
    cell_endcoordinates = pd.read_csv(last_coordinates_path)
    clean_OnPath_column_to_path_ID_n_label_AMCs(cell_endcoordinates) 
    
    # get soma coordinates
    soma_coordinates = cell_allcoordinates[cell_allcoordinates['path_ID'] == 1]
     
    # get all path_IDs
    path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()
    
    # set dict for all coordinates
    cell_coordintes_dict = {'all_coor' : cell_allcoordinates,
                            'pri_coor' : cell_primarycoordinates,
                            'end_coor' : cell_endcoordinates,
                            'soma_coor': soma_coordinates,
                            'path_IDs' : path_IDs}
   
    # # add to dict 
    coordinates_dict[cell_ID] = cell_coordintes_dict


# %% plot coordinates

if vplots:   
    # load plotting function    
    from cellmorphology.AMC_analysis.plot_AMC_cellcoordinates_analysis import plot_cellcoordinates
    
    # plot
    for cell_ID in tqdm(cell_IDs):  
        plot_cellcoordinates(cell_ID = cell_ID, cell_coordinates = coordinates_dict[cell_ID])


# %% get height, width and depth per type

print('height, width, and depth ...')

# create list of columns for dataframe
columns = [ntype + '-' + col for ntype in neurite_types for col in ['x_min', 'y_min', 'z_min', 'width', 'height', 'depth']]

# define output
height_width_depth = pd.DataFrame(columns=columns, 
                                  index = cell_IDs)
height_width_depth.index.name = 'cell_ID'


for cell_ID in tqdm(cell_IDs):  
    
    # get cell coordinates
    cell_allcoordinates = coordinates_dict[cell_ID]['all_coor']

    # get coordinates min and max
    for ntype in neurite_types:
        
        # limit coordinates to specific neurite type
        if ntype == 'neurites':
            cell_allcoordinates_pertype = cell_allcoordinates
        elif ntype == 'dendrites':
            cell_allcoordinates_pertype = cell_allcoordinates[cell_allcoordinates['path_label'] != 'axons']
        else:
            cell_allcoordinates_pertype = cell_allcoordinates[cell_allcoordinates['path_label'] == ntype]
        
        # get x, y, z min and max
        xmin = cell_allcoordinates_pertype['X'].min()
        ymin = cell_allcoordinates_pertype['Y'].min()
        zmin = cell_allcoordinates_pertype['Z'].min()
        xmax = cell_allcoordinates_pertype['X'].max()
        ymax = cell_allcoordinates_pertype['Y'].max()
        zmax = cell_allcoordinates_pertype['Z'].max()
        
        # write to dataframe
        # mins
        height_width_depth.at[cell_ID, ntype + '-' + 'x_min'] = xmin
        height_width_depth.at[cell_ID, ntype + '-' + 'y_min'] = ymin
        height_width_depth.at[cell_ID, ntype + '-' + 'z_min'] = zmin
        
        # calc width and height
        height_width_depth.at[cell_ID, ntype + '-' + 'width'] = xmax - xmin
        height_width_depth.at[cell_ID, ntype + '-' + 'height'] = ymax - ymin
        height_width_depth.at[cell_ID, ntype + '-' + 'depth'] = zmax - zmin

# save dataframe
height_width_depth.to_excel(join(AMCs_metrics_dir, 'height_width_depth.xlsx'),
                            index_label = 'cell_ID')


# %% plot height, width, and depth per type

if vplots:      
    
    # plot with height, width, and depth
    for cell_ID in tqdm(cell_IDs):  
        plot_cellcoordinates(cell_ID = cell_ID, 
                             cell_coordinates = coordinates_dict[cell_ID], 
                             cell_hwd = height_width_depth.loc[cell_ID, :])


# %% number of primary & terminal points & bifurcation ratio

n_primary_points = 

# %%

# init plotting
from cellmorphology.cellmorph_functions.initialize_AMC_cellmorph_plotting import *


plt.figure(dpi = 300)


plt.scatter(x = cell_allcoordinates.loc[:, 'X'],
            y = cell_allcoordinates.loc[:, 'Y'],
            color = 'gray',
            s = 0.25)


for path_i in cell_primarycoordinates.index.to_list():
    plt.scatter(x = cell_primarycoordinates.at[path_i, 'X'],
                y = cell_primarycoordinates.at[path_i, 'Y'],
                color = neurite_color_dict[cell_primarycoordinates.at[path_i, 'path_label']],
                s = 15,
                marker = 'x',
                linewidths=0.5)

for path_i in cell_endcoordinates.index.to_list():
    plt.scatter(x = cell_endcoordinates.at[path_i, 'X'],
                y = cell_endcoordinates.at[path_i, 'Y'],
                color = neurite_color_dict[cell_endcoordinates.at[path_i, 'path_label']],
                s = 3)

plt.scatter(x = soma_coordinates.at[0, 'X'],
            y = soma_coordinates.at[0, 'Y'],
            color = neurite_color_dict[soma_coordinates.at[0, 'path_label']])

plt.gca().invert_yaxis()

plt.show()