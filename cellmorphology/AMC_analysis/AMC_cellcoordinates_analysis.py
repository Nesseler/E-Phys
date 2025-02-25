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
cell_IDs = ['Exp-161', 'Exp-162']

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


for cell_ID in cell_IDs:       
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


# %% plot types of neurites

if vplots:   
    # load plotting function    
    from cellmorphology.AMC_analysis.plot_AMC_cellcoordinates_analysis import plot_neurite_types
    
    # plot
    for cell_ID in cell_IDs:  
        plot_neurite_types(cell_ID = cell_ID, 
                           cell_coordinates = coordinates_dict[cell_ID])


# %% get height, width and depth

print('\nheight, width, and depth ...')

# define output
height_width_depth = pd.DataFrame(columns = [ntype + '-' + col for ntype in neurite_types for col in ['x_min', 'y_min', 'z_min', 'width', 'height', 'depth']], 
                                  index = cell_IDs)
height_width_depth.index.name = 'cell_ID'


for cell_ID in cell_IDs:  
    
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


# %% plot height, width, and depth

if vplots:
    # load plotting function    
    from cellmorphology.AMC_analysis.plot_AMC_cellcoordinates_analysis import plot_cellcoordinates
    
    # plot with height, width, and depth
    for cell_ID in tqdm(cell_IDs):  
        plot_cellcoordinates(cell_ID = cell_ID, 
                             cell_coordinates = coordinates_dict[cell_ID], 
                             cell_hwd = height_width_depth.loc[cell_ID, :])


# %% total_cable_length

print('\ntotal cable length ...')

# import necessary function
from cellmorphology.cellmorph_functions.cellmorph_functions import calc_length_of_branch

# define output
total_cable_length = pd.DataFrame(columns = [ntype + '-total_cable_length' for ntype in neurite_types], 
                                  index = cell_IDs)
total_cable_length.index.name = 'cell_ID'


for cell_ID in cell_IDs:
    
    # get coordinates of cell
    cell_allcoordinates = coordinates_dict[cell_ID]['all_coor']
       
    # get all path_IDs of cell
    path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()
    
    # create dataframe to save path measurements
    path_measurements = pd.DataFrame(columns = ['path_length', 'path_label'], 
                                     index = path_IDs)
    path_measurements.index.name = 'path_IDs'
    
    # loop through path_IDs
    for path_ID in path_IDs:
        
        # coordinates and label of current path
        path_coordinates = cell_allcoordinates[cell_allcoordinates['path_ID'] == path_ID]
        path_label = cell_allcoordinates[cell_allcoordinates['path_ID'] == path_ID]['path_label'].drop_duplicates().values[0]
    
        # calculate the euclidean distance of current path
        path_length = calc_length_of_branch(path_coordinates.drop(columns = ['path_ID', 'path_label']))
    
        # write to dataframe
        path_measurements.at[path_ID, 'path_length'] = path_length
        path_measurements.at[path_ID, 'path_label'] = path_label

    # calc total cable length for cell
    for ntype in neurite_types:
        
        # limit coordinates to specific neurite type
        if ntype == 'neurites':
            tcl = path_measurements['path_length'].sum()
        elif ntype == 'dendrites':
            tcl = path_measurements.query(f'path_label != "axons"')['path_length'].sum()
        else:
            tcl = path_measurements[path_measurements['path_label'] == ntype]['path_length'].sum()
        
        # write to dataframe
        total_cable_length.at[cell_ID, f'{ntype}-total_cable_length'] = tcl

# save dataframe
total_cable_length.to_excel(join(AMCs_metrics_dir, 'total_cable_length.xlsx'),
                            index_label = 'cell_ID')


# %% reconstruct terminal branches
    # from end point (terminal point) to soma
    
print('\nreconstructing terminal branches ...')
    
# define dictionary to keep terminal branches measurements for all cells
all_terminal_branches = dict.fromkeys(cell_IDs)

# import necessary functions
from cellmorphology.cellmorph_functions.cellmorph_functions import find_parent_path
from cellmorphology.cellmorph_functions.cellmorph_functions import angle
from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_analysis_dir


for cell_ID in cell_IDs:

    # get coordinates of cell
    cell_allcoordinates = coordinates_dict[cell_ID]['all_coor']
    cell_terminalcoordinates = coordinates_dict[cell_ID]['end_coor']
    cell_somacoordinates = coordinates_dict[cell_ID]['soma_coor'].loc[0, :]
    
    # reset index
    cell_terminalcoordinates = cell_terminalcoordinates.set_index('path_ID')

    # get path IDs
    path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()
    terminal_pathIDs = cell_terminalcoordinates.index.astype(int).to_list()
    
    # define output
    terminal_branches = pd.DataFrame(columns = ['x_end', 'y_end', 'z_end', 'path_label', 'pathIDs_2soma', 
                                                'x_diff', 'y_diff', 'z_diff', 
                                                'angle_deg', 'angle_rad', 
                                                'length', 'euc_dist', 'contraction'],
                                     index = terminal_pathIDs)
    terminal_branches.index.name = 'terminal_pathIDs'
    
    # raise warning if soma is included
    if 1 in terminal_pathIDs:
        warnings.warn('Soma coordinates in terminal points coordinates!')
        
    # set reference vector
    reference_vector = pd.Series({'X' : 1., 'Y' : 0.})
    
    # precompute dictionary for path access
    allcoor_paths_dict = {pathID : group for pathID, group in cell_allcoordinates.groupby('path_ID')}

    # iterate through terminal paths
    for terminal_pathID in tqdm(terminal_pathIDs):
        
        ### orientation of terminal points ###
        
        # get terminal point coordinates
        terminal_point_coor = cell_terminalcoordinates.loc[terminal_pathID, :]
        
        # get vector from coordinates
        terminal_xy_vector = terminal_point_coor.drop(index = ['Z', 'path_label'])
        
        # write to dataframe
        terminal_point = {'x_end': terminal_point_coor['X'],
                          'y_end': terminal_point_coor['Y'],
                          'z_end': terminal_point_coor['Z'],
                          'x_diff': terminal_point_coor['X'] - cell_somacoordinates['X'],
                          'y_diff': terminal_point_coor['Y'] - cell_somacoordinates['Y'],
                          'z_diff': terminal_point_coor['Z'] - cell_somacoordinates['Z'],
                          'path_label': terminal_point_coor['path_label']}
        
        terminal_branches.loc[terminal_pathID, list(terminal_point.keys())] = list(terminal_point.values())
        
        # calc difference vector
        xy_diff_vector = terminal_branches.loc[terminal_pathID, ['x_diff','y_diff']]
        
        # check quadrant of difference vector
        if xy_diff_vector['y_diff'] > 0:
            deg = 360 - math.degrees(angle(reference_vector, xy_diff_vector))
        else:
            deg = math.degrees(angle(reference_vector, xy_diff_vector))
    
        # write angles to dataframe
        terminal_branches.at[terminal_pathID, 'angle_deg'] = deg
        terminal_branches.at[terminal_pathID, 'angle_rad'] = math.radians(deg)
    
        ### terminal branch reconstruction ###
        
        # coordinates of terminal branch to reconstruct until soma
        terminal_path_coor = allcoor_paths_dict[terminal_pathID]
        
        # order of coordinates starts with terminal branch endpoint 
        # to be reversed later when calculating the entire branch
        branch_terminal2soma = pd.DataFrame()
        
        # reversed and reindex coordinates of terminal path
        branch_terminal2soma = terminal_path_coor[::-1].reset_index(drop=True)
            
        # start while loop from terminal path
        # while loop finishes when soma path id is reached
        path_ID = terminal_pathID
        parent_pathID = terminal_pathID
        
        # set list of path ID that describe terminal branch
        terminal_branch_pathIDs_2soma = list()
        
        while parent_pathID != 1:
            
            # find parent of current path
            parent_pathID, intersect_index = find_parent_path(path_ID, path_IDs, allcoor_paths_dict)
                        
            # get parent path coordinates
            parent_path_coor = allcoor_paths_dict[parent_pathID]
        
            # get all indices of parent path until intersection index
            parent_indices = [i for i in parent_path_coor.index.to_list() if i <= intersect_index]
            
            # get all coordinates until intersection point
            parent_path_coor_until_intersect = cell_allcoordinates.loc[parent_indices]
            
            # reversed and reindex
            rr_parent_path_coor_until_intersect = parent_path_coor_until_intersect[::-1].reset_index(drop=True)
            
            # concatenate to all coordinates from terminal to soma
            branch_terminal2soma = pd.concat([branch_terminal2soma, rr_parent_path_coor_until_intersect])
            
            # reset index following concatenation
            branch_terminal2soma.reset_index(drop=True, inplace = True)
            
            # write parent path to list
            terminal_branch_pathIDs_2soma.append(path_ID)
            
            # set path_ID to parent for next step in while loop
            path_ID = parent_pathID
          
        # write list of path IDs to dataframe
        terminal_branches.at[terminal_pathID ,'pathIDs_2soma'] = terminal_branch_pathIDs_2soma
        
        # calculate length
        terminal_branch_length = calc_length_of_branch(branch_terminal2soma.drop(columns = ['path_ID', 'path_label']))
        
        # calculate euclidean distance between origin and terminal of branch
        line_coordinates = branch_terminal2soma.iloc[[0, -1]]
        
        # can be calculated with same function but just two coordinates
        terminal_branch_euc = calc_length_of_branch(line_coordinates.drop(columns = ['path_ID', 'path_label']))
        
        # write to dataframe
        terminal_branches.at[terminal_pathID ,'length'] = terminal_branch_length
        terminal_branches.at[terminal_pathID ,'euc_dist'] = terminal_branch_euc
        terminal_branches.at[terminal_pathID ,'contraction'] = terminal_branch_euc / terminal_branch_length        
        
    # write dataframe to dictionary
    all_terminal_branches[cell_ID] = terminal_branches
    
    # save dataframe
    terminal_branches.to_excel(join(AMCs_analysis_dir, 'metrics-terminal_branches' , f'{cell_ID}-terminal_branches.xlsx'),
                                index_label = 'terminal_pathIDs')
    

# %% plot terminal branches

if vplots:
    # load plotting function    
    from cellmorphology.AMC_analysis.plot_AMC_cellcoordinates_analysis import plot_all_terminal_branches
    
    # plot with height, width, and depth
    for cell_ID in tqdm(cell_IDs):  
        plot_all_terminal_branches(cell_ID = cell_ID,
                                   terminal_branches = all_terminal_branches[cell_ID],
                                   cell_coordinates = coordinates_dict[cell_ID])


# %% number of primary and terminal points

print('\nnumber of primary and terminal points ...')

# initialize dataframes
n_primary = pd.DataFrame(columns = [f'n_primary-{ntype}' for ntype in neurite_types],
                          index = cell_IDs)
n_terminal = pd.DataFrame(columns = [f'n_terminal-{ntype}' for ntype in neurite_types],
                          index = cell_IDs)

# rename index columns
for df in [n_primary, n_terminal]:
    df.index.name = 'cell_ID'

# iterate through cells
for cell_ID in cell_IDs:
    
    # get cell coordinates
    cell_pricoordinates = coordinates_dict[cell_ID]['pri_coor']
    cell_endcoordinates = coordinates_dict[cell_ID]['end_coor']
    
    for ntype in neurite_types:
        
        # limit coordinates to specific neurite type
        if ntype == 'neurites':
            cell_pricoordinates_pertype = cell_pricoordinates
            cell_endcoordinates_pertype = cell_endcoordinates
            
        elif ntype == 'dendrites':
            cell_pricoordinates_pertype = cell_pricoordinates[cell_pricoordinates['path_label'] != 'axons']
            cell_endcoordinates_pertype = cell_endcoordinates[cell_endcoordinates['path_label'] != 'axons']
            
        else:
            cell_pricoordinates_pertype = cell_pricoordinates[cell_pricoordinates['path_label'] == ntype]
            cell_endcoordinates_pertype = cell_endcoordinates[cell_endcoordinates['path_label'] == ntype]       
        
        # get number of points
        n_primary_pertype = cell_pricoordinates_pertype.shape[0]
        n_terminal_pertype = cell_endcoordinates_pertype.shape[0]
        
        # write to dataframe
        n_primary.at[cell_ID, f'n_primary-{ntype}'] = n_primary_pertype
        n_terminal.at[cell_ID, f'n_terminal-{ntype}'] = n_terminal_pertype
    
# save dataframe
n_primary.to_excel(join(AMCs_analysis_dir, 'metrics' , f'n_primary.xlsx'),
                   index_label = 'cell_ID')
n_terminal.to_excel(join(AMCs_analysis_dir, 'metrics' , f'n_terminal.xlsx'),
                   index_label = 'cell_ID')
  
            
# %% number of stems and bifuracation ratio

print('\nnumber of stems and bifuracation ratio ...')

# define output
n_stems = pd.DataFrame(columns = [f'n_stems-{ntype}' for ntype in neurite_types],
                       index = cell_IDs)
bifurcation_ratios = pd.DataFrame(columns = [f'bifurcation_ratio-{ntype}' for ntype in neurite_types],
                                  index = cell_IDs)

# rename index
for df in [n_stems, bifurcation_ratios]:
    df.index.name = 'cell_ID'

# iterate through cells
for cell_ID in cell_IDs:

    # get terminal branches dataframe
    terminal_branches = all_terminal_branches[cell_ID]
    
    # get all coordinates for cell
    cell_allcoordinates = coordinates_dict[cell_ID]['all_coor']
    
    # get all coordinates grouped by pathID
    allcoor_paths_dict = {pathID : group for pathID, group in cell_allcoordinates.groupby('path_ID')}
    
    # get assigned path labels for each path
    pathlabels_dict = {pathID : group['path_label'].iat[0] for pathID, group in cell_allcoordinates.groupby('path_ID')}

    # iterate through neurite type
    for ntype in neurite_types:
        
        # number of stems
        
        # get only terminal_branches of neurite type
        if ntype == 'neurites':
            ntype_terminalbranches = terminal_branches
        elif ntype == 'dendrites':
            ntype_terminalbranches = terminal_branches[terminal_branches['path_label'] != 'axons']
        else:
            ntype_terminalbranches = terminal_branches[terminal_branches['path_label'] == ntype]
      
        # create a list of start points
        stems = pd.DataFrame(index = ['X', 'Y', 'Z'])

        # iterate through terminal points
        for ti, terminal_pathID in enumerate(ntype_terminalbranches.index.to_list()):
            
            # get the paths to soma
            pathIDs_2soma = ntype_terminalbranches.at[terminal_pathID, 'pathIDs_2soma']
        
            # get terminal label
            terminal_label = ntype_terminalbranches.at[terminal_pathID, 'path_label']
    
            # iterate through path to soma
            for pathID_idx, path_ID in enumerate(pathIDs_2soma):
                
                # get path label
                cur_pathlabel = pathlabels_dict[path_ID]
                
                # check if paths labels match
                if cur_pathlabel != terminal_label:
                    
                    # get previouse path ID
                    prev_pathID = pathIDs_2soma[pathID_idx-1]
                    
                    # get stem coordinates
                    stems_coor = allcoor_paths_dict[prev_pathID].loc[:, ['X', 'Y', 'Z']].iloc[0, :]

                    break
                
                # check if loop reaches all paths to soma 
                elif pathID_idx == len(pathIDs_2soma)-1:
                    
                    # define stem coordinates
                    stems_coor = allcoor_paths_dict[path_ID].loc[:, ['X', 'Y', 'Z']].iloc[0, :]
                
            # write coordinates of first node to stems dataframe
            stems.loc[:, stems_coor.name] = stems_coor
                            
        # get number of stems from dataframe index and write to dataframe
        n_stems.at[cell_ID, f'n_stems-{ntype}'] = stems.T.index.shape[0]

        # bifurcation ratio
        
        # check if cell has neurite type
        if ntype_terminalbranches.shape[0] > 0:
            
            # calc bifurcation ratio
            bifurcation_ratio_pertype = n_terminal.at[cell_ID, f'n_terminal-{ntype}'] / n_stems.at[cell_ID, f'n_stems-{ntype}']
            
        else:
            bifurcation_ratio_pertype = np.nan
        
        # write to dataframe
        bifurcation_ratios.at[cell_ID, f'bifurcation_ratio-{ntype}'] = bifurcation_ratio_pertype

# save dataframe
n_stems.to_excel(join(AMCs_analysis_dir, 'metrics' , f'n_stems.xlsx'),
                 index_label = 'cell_ID')
bifurcation_ratios.to_excel(join(AMCs_analysis_dir, 'metrics' , f'bifurcation_ratios.xlsx'),
                            index_label = 'cell_ID')

# %% plot primary & terminal points & bifurcation ratio

if vplots:   
    # load plotting function    
    from cellmorphology.AMC_analysis.plot_AMC_cellcoordinates_analysis import plot_endpoints
    
    # plot
    for cell_ID in tqdm(cell_IDs):  
        plot_endpoints(cell_ID = cell_ID, 
                        cell_coordinates = coordinates_dict[cell_ID],
                        n_primary = n_primary,
                        n_terminal = n_terminal,
                        n_stems = n_stems,
                        bifurcation_ratios = bifurcation_ratios)


# %% terminal branches orientation



