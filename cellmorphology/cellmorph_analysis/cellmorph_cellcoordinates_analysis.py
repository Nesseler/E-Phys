# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 10:14:45 2025

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# get directories
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_coordinates_dir, cellmorph_analysis_dir

# get parameters & functions
from cellmorphology.cellmorph_functions.cellmorph_parameters import field_of_view
from cellmorphology.cellmorph_functions.cellmorph_polarplot_functions import assign_bin2branch, get_orientation_occurances

# get files in traces folder
from functions.get_onlyfiles_list import get_onlyfiles_list
onlyfiles = get_onlyfiles_list(cellmorph_coordinates_dir)

# get cell_IDs from list of all coordinates filenames
cell_IDs = [filename[:5] for filename in [f for f in onlyfiles if 'all_coordinates' in f]]

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# set neurite types
neurite_types = ['neurites', 'dendrites', 'axons']

vplots = True


# %% check for new cells to be analyzed

# load anaylsis worksheet
from parameters.directories_win import table_file
analyzed = pd.read_excel(table_file,
                         sheet_name = 'analyzed',
                         index_col = 'cell_ID')

# get list of cell_IDs already analyzed
analyzed_cell_IDs = analyzed.loc[analyzed['cellmorpho_cellcoordinates'].notna()]['cellmorpho_cellcoordinates'].index.to_list()

# redefine cell_IDs list
cellIDs_toAnalyze = [cell_ID for cell_ID in cell_IDs if cell_ID not in analyzed_cell_IDs]

# raise error
if len(cellIDs_toAnalyze) == 0:
    raise ValueError('Nothing new to analyze!')


# cellIDs_toAnalyze = ['E-217', 'E-218', 'E-219', 'E-222', 'E-238', 'E-239', 'E-240', 'E-242', 'E-243', 'E-244', 'E-245', 'E-236', 'E-277', 'E-280', 'E-284', 'E-290', 'E-292', 'E-296', 'E-297']
cellIDs_toAnalyze = ['E-111', 'E-137', 'E-147', 'E-162']
# cellIDs_toAnalyze = cell_IDs

# %% define output

height_width_depth = pd.DataFrame(columns = [ntype + '-' + col for ntype in neurite_types for col in ['x_min', 'y_min', 'z_min', 'width', 'height', 'depth']], 
                                  index = cellIDs_toAnalyze)

total_cable_length = pd.DataFrame(columns = [ntype + '-total_cable_length' for ntype in neurite_types], 
                                  index = cellIDs_toAnalyze)

n_primary = pd.DataFrame(columns = [f'n_primary-{ntype}' for ntype in neurite_types],
                          index = cellIDs_toAnalyze)

n_terminal = pd.DataFrame(columns = [f'n_terminal-{ntype}' for ntype in neurite_types],
                          index = cellIDs_toAnalyze)

bifurcation_ratios = pd.DataFrame(columns = [f'bifurcation_ratio-{ntype}' for ntype in neurite_types],
                                  index = cellIDs_toAnalyze)

axon_origins = pd.DataFrame(columns = ['axon_origin'],
                            index = cellIDs_toAnalyze)

circ_stats = pd.DataFrame(columns = [f'terminal_angles_rad-{ntype}' for ntype in neurite_types] 
                                    + [f'{stat}_{unit}-{ntype}'  for unit in ['rad', 'deg'] for ntype in neurite_types for stat in ['circmean', 'circstd', 'circvar']],
                          index = cellIDs_toAnalyze)

from cellmorphology.cellmorph_functions.cellmorph_polarplot_functions import orientation_labels
orientations = pd.DataFrame(columns = [f'{ntype}-{orientation}' for ntype in neurite_types for orientation in orientation_labels],
                            index = cellIDs_toAnalyze)

AcD = pd.DataFrame(columns = ['dendrites_pathIDs', 'axons_pathIDs', 'terminal_pathIDs', 'dendrites-angles_rad' , 'circmean_rad', 'circstd_rad', 'circmean_deg'] + orientation_labels,
                   index = cellIDs_toAnalyze)

# rename index
for df in [height_width_depth, total_cable_length, n_primary, n_terminal, bifurcation_ratios, axon_origins, circ_stats, AcD, orientations]:
    df.index.name = 'cell_ID'

# %% start analysis (per cell)

for cell_ID in tqdm(cellIDs_toAnalyze):
        
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
     
    # check if x coordinates need to be flipped
    if MetaData.at[cell_ID, 'to_be_x_flipped']:
        cell_allcoordinates.loc[:, 'X']     = field_of_view - cell_allcoordinates['X']
        cell_primarycoordinates.loc[:, 'X'] = field_of_view - cell_primarycoordinates['X']
        cell_endcoordinates.loc[:, 'X']     = field_of_view - cell_endcoordinates['X']
        soma_coordinates.loc[:, 'X']        = field_of_view - soma_coordinates['X']
    
    # get all path_IDs
    path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()
    
    # set dict for all coordinates
    cell_coordinates = {'all_coor' : cell_allcoordinates,
                        'pri_coor' : cell_primarycoordinates,
                        'end_coor' : cell_endcoordinates,
                        'soma_coor': soma_coordinates,
                        'path_IDs' : path_IDs}
   
    # plot coordinates
    if vplots:   
        # load plotting function    
        from cellmorphology.cellmorph_analysis.plot_cellcoordinates_analysis import plot_cellcoordinates
        
        # plot 
        plot_cellcoordinates(cell_ID = cell_ID, 
                             cell_coordinates = cell_coordinates)
        
        
    # %% height, width, depth
    
    # get cell coordinates
    cell_allcoordinates = cell_coordinates['all_coor']

    # get coordinates min and max
    for ntype in neurite_types:
        
        # limit coordinates to specific neurite type
        if ntype == 'neurites':
            cell_allcoordinates_pertype = cell_allcoordinates
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
    
    # plot height, width, depth
    if vplots:   
        # load plotting function    
        from cellmorphology.cellmorph_analysis.plot_cellcoordinates_analysis import plot_cellhwd
        
        # plot 
        plot_cellhwd(cell_ID = cell_ID, 
                     cell_coordinates = cell_coordinates,
                     cell_hwd = height_width_depth.loc[cell_ID, :])
    

    # %% total_cable_length
    
    # import necessary function
    from cellmorphology.cellmorph_functions.cellmorph_functions import calc_length_of_branch
    
    # get coordinates of cell
    cell_allcoordinates = cell_coordinates['all_coor']
       
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
        else:
            tcl = path_measurements[path_measurements['path_label'] == ntype]['path_length'].sum()
        
        # write to dataframe
        total_cable_length.at[cell_ID, f'{ntype}-total_cable_length'] = tcl
    
    
    # %% reconstruct terminal branches
        # from end point (terminal point) to soma
    
    # import necessary functions
    from cellmorphology.cellmorph_functions.cellmorph_functions import find_parent_path
    from cellmorphology.cellmorph_functions.cellmorph_functions import angle
    
    # get coordinates of cell
    cell_allcoordinates = cell_coordinates['all_coor']
    cell_terminalcoordinates = cell_coordinates['end_coor']
    cell_somacoordinates = cell_coordinates['soma_coor'].loc[0, :]
    
    # reset index
    cell_terminalcoordinates = cell_terminalcoordinates.set_index('path_ID')

    # get path IDs
    path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()
    terminal_pathIDs = cell_terminalcoordinates.index.astype(int).to_list()
    
    # get dict for all path and their labels
    unique_pathidc = cell_allcoordinates['path_ID'].drop_duplicates().index
    pathlabel_dict = dict(zip(path_IDs, cell_allcoordinates.loc[unique_pathidc, 'path_label']))
    
    # define output
    terminal_branches = pd.DataFrame(columns = ['x_end', 'y_end', 'z_end', 'path_label', 
                                                'pathIDs_2soma', 'parent_idx_allcoor',
                                                'x_diff', 'y_diff', 'z_diff', 
                                                'angle_deg', 'angle_rad', 
                                                'length', 'euc_dist', 'contraction',
                                                'bin_id'],
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
    for terminal_pathID in terminal_pathIDs: #[39]: #
        
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
            
        # convert to radians
        rad = math.radians(deg)
    
        # write angles to dataframe
        terminal_branches.at[terminal_pathID, 'angle_deg'] = deg
        terminal_branches.at[terminal_pathID, 'angle_rad'] = rad
    
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
        
        # set list of path ID that describe terminal branch (and list for start and stop index in all coordinates)
        terminal_branch_pathIDs_2soma = list()
        terminal_branch_parent_idx_allcoor = [(terminal_path_coor.index[0], terminal_path_coor.index[-1])]
        
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
            terminal_branch_parent_idx_allcoor.append((parent_indices[0], parent_indices[-1]))
            
            # set path_ID to parent for next step in while loop
            path_ID = parent_pathID
          
        # write list of path IDs to dataframe
        terminal_branches.at[terminal_pathID ,'pathIDs_2soma'] = terminal_branch_pathIDs_2soma
        terminal_branches.at[terminal_pathID ,'parent_idx_allcoor'] = terminal_branch_parent_idx_allcoor
        
        # calculate length
        terminal_branch_length = calc_length_of_branch(branch_terminal2soma.drop(columns = ['path_ID', 'path_label']))
        
        # calculate euclidean distance between origin and terminal of branch
        line_coordinates = branch_terminal2soma.iloc[[0, -1]]
        
        # can be calculated with same function but just two coordinates
        terminal_branch_euc = calc_length_of_branch(line_coordinates.drop(columns = ['path_ID', 'path_label']))
        
        # get orientation histogram bin assignment
        assigned_bin = assign_bin2branch(rad)
        
        # write to dataframe
        terminal_branches.at[terminal_pathID ,'length'] = terminal_branch_length
        terminal_branches.at[terminal_pathID ,'euc_dist'] = terminal_branch_euc
        terminal_branches.at[terminal_pathID ,'contraction'] = terminal_branch_euc / terminal_branch_length        
        terminal_branches.at[terminal_pathID ,'bin_id'] = assigned_bin
    
    # save dataframe
    terminal_branches.to_excel(join(cellmorph_analysis_dir, 'metrics-terminal_branches' , f'{cell_ID}-terminal_branches.xlsx'),
                                index_label = 'terminal_pathIDs')
        
    # plot terminal branches
    if False:
        # load plotting function    
        from cellmorphology.cellmorph_analysis.plot_cellcoordinates_analysis import plot_all_terminal_branches
        
        # plot terminal branches
        plot_all_terminal_branches(cell_ID = cell_ID,
                                   cell_coordinates = cell_coordinates,
                                   terminal_branches = terminal_branches)
        
        
    # %% orientation histogram and circular mean
    
    for ntype in neurite_types:
        
        # get orientation per type
        if ntype == 'neurites':
            angles_rad = terminal_branches['angle_rad'].to_list()
        else:
            angles_rad = terminal_branches[terminal_branches['path_label'] == ntype]['angle_rad'].to_list()

        # calc occurances
        orientation_occu_pertype = get_orientation_occurances(angles_rad)

        # write to dataframe
        orientations.loc[cell_ID, [f'{ntype}-{orientation}' for orientation in orientation_labels]] = orientation_occu_pertype
        
        # calc circular stats
        circmean_pertype = sc.stats.circmean(angles_rad)
        circstd_pertype = sc.stats.circstd(angles_rad)
        circvar_pertype = sc.stats.circvar(angles_rad)
        
        # write to dataframe
        circ_stats.at[cell_ID, f'terminal_angles_rad-{ntype}'] = angles_rad
        circ_stats.at[cell_ID, f'circmean_rad-{ntype}'] = circmean_pertype
        circ_stats.at[cell_ID, f'circstd_rad-{ntype}'] = circstd_pertype
        circ_stats.at[cell_ID, f'circvar_rad-{ntype}'] = circvar_pertype
        circ_stats.at[cell_ID, f'circmean_deg-{ntype}'] = np.rad2deg(circmean_pertype)
        circ_stats.at[cell_ID, f'circstd_deg-{ntype}'] = np.rad2deg(circstd_pertype)
        circ_stats.at[cell_ID, f'circvar_deg-{ntype}'] = np.rad2deg(circvar_pertype)
        
    # plot cell polar plots
    if vplots:
        # load plotting function    
        from cellmorphology.cellmorph_analysis.plot_cellcoordinates_analysis import plot_polar_plot_abs
        
        # plot terminal branches
        plot_polar_plot_abs(cell_ID = cell_ID,
                            cell_coordinates = cell_coordinates,
                            terminal_branches = terminal_branches,
                            circ_stats = circ_stats)
    
    
    # %% number of primary, terminal points, bifurcation ratio
    
    # get cell coordinates
    cell_pricoordinates = cell_coordinates['pri_coor']
    cell_endcoordinates = cell_coordinates['end_coor']
    
    for ntype in neurite_types:
        
        # limit coordinates to specific neurite type
        if ntype == 'neurites':
            cell_pricoordinates_pertype = cell_pricoordinates
            cell_endcoordinates_pertype = cell_endcoordinates
                       
        else:
            cell_pricoordinates_pertype = cell_pricoordinates[cell_pricoordinates['path_label'] == ntype]
            cell_endcoordinates_pertype = cell_endcoordinates[cell_endcoordinates['path_label'] == ntype]       
        
        # get number of points
        n_primary_pertype = cell_pricoordinates_pertype.shape[0]
        n_terminal_pertype = cell_endcoordinates_pertype.shape[0]
        
        # write to dataframe
        n_primary.at[cell_ID, f'n_primary-{ntype}'] = n_primary_pertype
        n_terminal.at[cell_ID, f'n_terminal-{ntype}'] = n_terminal_pertype
        
        ### bifurcation ratio
        if ntype == 'axons':
            
            # condition for dendritic axon
            if n_primary_pertype == 0 and n_terminal_pertype > 0:
                bif_ratio_pertype = n_terminal_pertype / 1
                cell_axonorigin = 'dendritic'
                
            # condition for somatic axon
            elif n_primary_pertype == 1 and n_terminal_pertype > 0:
                bif_ratio_pertype = n_terminal_pertype / n_primary_pertype
                cell_axonorigin = 'somatic'
            
            # condition for no axons
            elif n_primary_pertype == 0 and n_terminal_pertype == 0:
                bif_ratio_pertype = np.nan
                cell_axonorigin = np.nan
                
            else:
                raise ValueError('Bifurcation ratio not calculated')
            
        else:
            bif_ratio_pertype = n_terminal_pertype / n_primary_pertype
            
        # write to dataframe
        bifurcation_ratios.at[cell_ID, f'bifurcation_ratio-{ntype}'] = bif_ratio_pertype
    
    # write to dataframe
    axon_origins.at[cell_ID, 'axon_origin'] = cell_axonorigin
        
    # plot number of primary and terminal points and bifurcation ratios
    if vplots:
        # load plotting function    
        from cellmorphology.cellmorph_analysis.plot_cellcoordinates_analysis import plot_endpoints
        
        # plot terminal branches
        plot_endpoints(cell_ID = cell_ID,
                       cell_coordinates = cell_coordinates,
                       n_primary = n_primary, 
                       n_terminal = n_terminal,
                       bifurcation_ratios = bifurcation_ratios)
    
    
    # %% get axon carrying dendrite
    
    # precompute dictionary for path label access
    terminalbranches_pertype = {pathID : group for pathID, group in terminal_branches.groupby('path_label')}
    
    # check axon origin type
    if cell_axonorigin == 'dendritic':
    
        # get axonic terminal branches
        axon_terminalbranches = terminalbranches_pertype['axons']
        
        # create empty list
        axon_stems = set()

        # iterate through axonic terminal branches    
        for terminal_pathID, terminal_branch in axon_terminalbranches.iterrows():
            
            # get last path_ID
            l_pathID = terminal_branch['pathIDs_2soma'][-1]
            
            # set as axon_stem
            if l_pathID not in axon_stems:
                axon_stems.add(l_pathID)
    
        # all terminal axons should be carried by same dendritic path
        if len(axon_stems) > 1:
            raise ValueError('More than one axon stem!')
        else:
            axonstem_pathID = list(axon_stems)[0]
    
        # create list of axon carrying dendrite terminal path_IDs
        AcD_terminalpathIDs = list()
    
        # get all terminal paths with same stem
        for terminal_pathID, terminal_branch in terminal_branches.iterrows():

            # add to list if in paths to soma
            if axonstem_pathID in terminal_branch['pathIDs_2soma']:
                AcD_terminalpathIDs.append(terminal_pathID)
                
        # create dataframe of terminal branches
        AcD_terminalbranches = terminal_branches.loc[AcD_terminalpathIDs, :]
    
        # create list of axon carrying dendrite path_IDs
        AcD_pathIDs = list()
        Axon_pathIDs = list()
    
        # get all path_IDs on axon carrying dendrites
        for terminal_pathID, terminal_branch in AcD_terminalbranches.iterrows():
            
            # iterate through paths 2 soma
            for path_ID in terminal_branch['pathIDs_2soma'][::-1]:
                
                # get path label
                path_label = pathlabel_dict[path_ID]
            
                # add to list if in paths to soma
                if path_label == 'dendrites':
                    
                    if path_ID not in AcD_pathIDs:
                        AcD_pathIDs.append(path_ID)
                
                elif path_label == 'axons':

                    if path_ID not in Axon_pathIDs:
                        Axon_pathIDs.append(path_ID)
        
        # write to dataframe
        AcD.at[cell_ID, 'dendrites_pathIDs'] = AcD_pathIDs
        AcD.at[cell_ID, 'axons_pathIDs'] = Axon_pathIDs
        AcD.at[cell_ID, 'terminal_pathIDs'] = AcD_terminalbranches.index.to_list()
        
        # calc circular mean for AcD terminal points
        AcDterminals_angles_rad = AcD_terminalbranches[AcD_terminalbranches['path_label'] == 'dendrites']['angle_rad'].to_list()
        circmean_rad = sc.stats.circmean(AcDterminals_angles_rad)
        circstd_rad = sc.stats.circstd(AcDterminals_angles_rad)
        circmean_deg = np.rad2deg(circmean_rad)
        
        # get orientation of dendrites
        AcD_occu = get_orientation_occurances(AcDterminals_angles_rad)
        
        # write to dataframe
        AcD.at[cell_ID, 'dendrites-angles_rad'] = AcDterminals_angles_rad
        AcD.at[cell_ID, 'circmean_rad'] = circmean_rad
        AcD.at[cell_ID, 'circstd_rad'] = circstd_rad
        AcD.at[cell_ID, 'circmean_deg'] = circmean_deg
        AcD.loc[cell_ID, orientation_labels] = AcD_occu
    
        # plot axon carrying dendrite plot
        if vplots:
            # load plotting function    
            from cellmorphology.cellmorph_analysis.plot_cellcoordinates_analysis import plot_AcD
            
            # plot terminal branches
            plot_AcD(cell_ID = cell_ID,
                     cell_coordinates = cell_coordinates,
                     AcD = AcD)
 

# %% save dataframes

# get directory
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_metrics_dir

# tobe saved ('filename' : variable)
export_vars = {'height_width_depth' : height_width_depth,
                'total_cable_length' : total_cable_length,
                'n_primary' : n_primary,
                'n_terminal' : n_terminal,
                'bifurcation_ratios' : bifurcation_ratios,
                'axon_origins' : axon_origins,
                'circ_stats' : circ_stats,
                'AcD' : AcD,
                'orientations' : orientations}

export_extension = '.xlsx'

# iterate through export variables
for export_name, export_var in export_vars.items():
    
    # get rows and cols
    rows = export_var.index.to_list()
    cols = export_var.columns.to_list()

    # try loading and writing or create new file
    try:
        loaded_export_var = pd.read_excel(join(cellmorph_metrics_dir, export_name + export_extension),
                                          index_col = 'cell_ID')
      
        try:
            # combine both dataframes
            loaded_export_var.loc[rows, cols] = export_var.loc[rows, cols].values
    
        # if values not yet in file, append/concat new values
        except KeyError:
            loaded_export_var = pd.concat([loaded_export_var, export_var.loc[rows, cols]], axis = 0)
                
        # sort by index
        loaded_export_var.sort_index(inplace = True)
        
        # save activity dataframe
        loaded_export_var.to_excel(join(cellmorph_metrics_dir, export_name + export_extension), 
                                    index_label = 'cell_ID')
    
    except FileNotFoundError:
        # save activity dataframe
        export_var.to_excel(join(cellmorph_metrics_dir, export_name + export_extension), 
                            index_label = 'cell_ID')

# %% update analyzed cells

# from functions.update_database import update_analyzed_sheet
    
# update_analyzed_sheet(cellIDs_toAnalyze, PGF = 'cellmorpho_cellcoordinates')


# %% calc population orientation

# re-load orientations dataframe to include all cells
orientations = pd.read_excel(join(cellmorph_metrics_dir, 'orientations.xlsx'),
                             index_col = 'cell_ID')

# initialize dataframe for populations
population_orientation = pd.DataFrame(columns = orientation_labels,
                                      index = [f'{ntype}-{region}-{measure}' for measure in ['abs', 'norm_toall', 'norm_totype'] for region in ['MeA', 'BAOT'] for ntype in neurite_types])

# for region in ['BAOT', 'MeA']:
    
#     # get region cell_IDs
#     region_cellIDs = MetaData[MetaData['Region'] == region].index.to_list()
    
#     # limit to reconstructed cells in region
#     region_cellIDs = [cell_ID for cell_ID in cell_IDs if cell_ID in region_cellIDs]
    
#     # get number of neurites
#     occu_neurites = orientations.loc[region_cellIDs, [f'neurites-{o_label}' for o_label in orientation_labels]]
    
#     # get total number
#     n_neurites = occu_neurites.sum(axis = 1)
    
#     for ntype in neurite_types:
        
#         # get orientation labels per type
#         orientation_labels_pertype = [f'{ntype}-{o_label}' for o_label in orientation_labels]
        
#         # get orientation occurances per region
#         occu_pertype_percell = orientations.loc[region_cellIDs, orientation_labels_pertype]
#         occu_pertype = occu_pertype_percell.sum(axis = 0)

#         # calc total number 
#         n_pertype = occu_pertype_percell.sum(axis = 1)
        
#         # normalize to total number per type
#         occu_norm_totype_percell = occu_pertype_percell.div(n_pertype, axis = 'rows')
        
#         # drop rows (cells) with zeros (condition for cells with no axons)
#         occu_norm_totype_percell = occu_norm_totype_percell.loc[(occu_norm_totype_percell!=0).any(axis=1)]
        
#         # calc average over all cells per bin
#         occu_norm_totype = occu_norm_totype_percell.mean(axis = 0)

#         # normalize to number of neurites
#         occu_norm_toall_percell = occu_pertype_percell.div(n_neurites, axis = 'rows')
#         occu_norm_toall = occu_norm_toall_percell.mean(axis = 0)

#         # write to dataframe
#         population_orientation.loc[f'{ntype}-{region}-abs', :] = occu_pertype.to_list()
#         population_orientation.loc[f'{ntype}-{region}-norm_totype', :] = occu_norm_totype.to_list()
#         population_orientation.loc[f'{ntype}-{region}-norm_toall', :] = occu_norm_toall.to_list()
        
# # save dataframe
# population_orientation.to_excel(join(cellmorph_metrics_dir, 'population_orientation.xlsx'), 
#                                 index_label = 'ntype-region-measure')


# %% calc population circular stats


