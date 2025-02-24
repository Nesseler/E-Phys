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

vplots = False

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
    
from cellmorphology.cellmorph_functions.initialize_AMC_cellmorph_plotting import *


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
        
    # save dataframe
    terminal_branches.to_excel(join(AMCs_analysis_dir, 'metrics-terminal_branches' , f'{cell_ID}-terminal_branches.xlsx'),
                                index_label = 'terminal_pathIDs')
    

# %%


# TODO: vplot endpoint orientation


# init plotting
from cellmorphology.cellmorph_functions.initialize_AMC_cellmorph_plotting import *


# field of view dimension
max_fov_xy = 590.76
max_fov_z = 300


# 
allcoor_paths_dict = {pathID : group for pathID, group in coordinates_dict[cell_ID]['all_coor'].groupby('path_ID')}


for terminal_branch_ID in [51]:#terminal_branches.index.to_list():

    # path_label = terminal_branches.at[terminal_branch_ID, 'path_label']
    pathIDs_2soma = terminal_branches.at[terminal_branch_ID, 'pathIDs_2soma'] 

    fig, axs = plt.subplots(nrows = 2,
                            ncols = 2,
                            layout = 'constrained',
                            figsize = get_figure_size(width = 100, height = 100),
                            width_ratios = [1, max_fov_z/max_fov_xy],
                            height_ratios = [1, max_fov_z/max_fov_xy],
                            sharey = 'row',
                            sharex = 'col',
                            dpi = 300)
    
    # flatten axes array
    axs = axs.flatten()
    
    # set figure title
    fig.suptitle(f'{cell_ID} terminal path-{terminal_branch_ID}',
                 fontsize = 9)

    for path_ID in pathIDs_2soma:
        
        path_label = allcoor_paths_dict[path_ID]['path_label'].unique()[0]
        
        scatter_paths = axs[0].scatter(x = allcoor_paths_dict[path_ID]['X'], 
                                       y = allcoor_paths_dict[path_ID]['Y'],
                                       s = 0.25, 
                                       zorder = 1,
                                       color = neurite_color_dict[path_label],
                                       alpha = 0.5)








    # edit axes
    # XY
    axs[0].text(x = 10, 
                y = 10, 
                s = 'XY', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[0].set_xlim([0, max_fov_xy])
    axs[0].set_ylim([max_fov_xy, 0])
    axs[0].set_ylabel('Height [µm]')
    axs[0].set_yticks(ticks = np.arange(0, max_fov_xy, 200))
    axs[0].set_yticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)
    
    # ZY
    axs[1].text(x = 10, 
                y = 10, 
                s = 'ZY', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[1].set_xlim([0, max_fov_z])
    axs[1].set_xlabel('')
    axs[1].set_xticks(ticks = np.arange(0, max_fov_z, 200))
    axs[1].set_xticks(ticks = np.arange(0, max_fov_z, 25), minor = True)
    
    # XZ
    axs[2].text(x = 10, 
                y = 10, 
                s = 'XZ', 
                ha = 'left', 
                va = 'top', 
                fontsize = 9)
    axs[2].set_xlim([0, max_fov_xy])
    axs[2].set_xlabel('Width [µm]')
    axs[2].set_ylim([max_fov_z, 0])
    axs[2].set_ylabel('Depth [µm]')
    axs[2].set_xticks(ticks = np.arange(0, max_fov_xy, 200))
    axs[2].set_xticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)
    axs[2].set_yticks(ticks = np.arange(0, max_fov_z, 200))
    axs[2].set_yticks(ticks = np.arange(0, max_fov_z, 25), minor = True)



    plt.show()



















# # %% reconstruction plot

# %matplotlib qt5

# # Initialise figure
# fig_ani, ax_ani = plt.subplots(1, 1, 
#                                figsize = get_figure_size(width=100, height = 100),
#                                layout = 'constrained',
#                                sharex=True)

# ## initialise scatter plot
# scatter_paths = ax_ani.scatter(x = [], 
#                                   y = [],
#                                   s = 0.25, 
#                                   zorder = 1,
#                                   color = 'r',
#                                   alpha = 0.5)

# scatter_finpaths = ax_ani.scatter(x = [], 
#                                   y = [],
#                                   s = 0.25, 
#                                   zorder = 0,
#                                   color = 'grey',
#                                   alpha = 0.5)

# finpaths = list()

# scatter_terminal = ax_ani.scatter(x = [], 
#                                      y = [],
#                                      s = 20, 
#                                      c = 'r',
#                                      zorder = 1,
#                                      marker = 'x')


# # plane label
# ax_ani.text(x = 10, 
#         y = 10, 
#         s = 'XY', 
#         ha = 'left', 
#         va = 'top', 
#         fontsize = 9)

# max_fov_xy = 590.76

# # edit axes
# ax_ani.set_ylim([0, max_fov_xy])
# ax_ani.set_ylabel('Height [µm]')
# ax_ani.set_yticks(ticks = np.arange(0, max_fov_xy, 200))
# ax_ani.set_yticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)

# ax_ani.set_xlim([0, max_fov_xy])
# ax_ani.set_xlabel('Width [µm]')
# ax_ani.set_xticks(ticks = np.arange(0, max_fov_xy, 200))
# ax_ani.set_xticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)

# # invert y axis
# ax_ani.invert_yaxis()



# scatter_soma = ax_ani.scatter(cell_somacoordinates['X'], 
#                          cell_somacoordinates['Y'], 
#                          s = 25, 
#                          c= 'grey', 
#                          zorder = 2)



        # scatter_terminal.set_offsets([terminal_point_coor['X'], terminal_point_coor['Y']])  
# scatter_paths.set_offsets(list(zip(branch_terminal2soma['X'], branch_terminal2soma['Y'])))        
        
# # plt.scatter(terminal_point_coor['X'], terminal_point_coor['Y'], s = 10, c = 'r', zorder = 1)
# # plt.scatter(cell_somacoordinates['X'], cell_somacoordinates['Y'], s = 25, c = 'grey', zorder = 2)
# # plt.scatter(branch_terminal2soma['X'], branch_terminal2soma['Y'], s = 0.25, zorder = 0)
# # plt.xlim([0, 590.76])
# # plt.ylim([590.76, 0])
# plt.pause(0.001)

# # scatter_paths

# finpaths = finpaths + list(zip(branch_terminal2soma['X'], branch_terminal2soma['Y']))
# scatter_finpaths.set_offsets(finpaths)

# %matplotlib inline

# plt.show()


# %%


# # %% number of primary & terminal points & bifurcation ratio

# # initialize dataframes
# n_primary = pd.DataFrame(columns = [f'n_primary-{ntype}' for ntype in neurite_types],
#                           index = cell_IDs)
# n_terminal = pd.DataFrame(columns = [f'n_terminal-{ntype}' for ntype in neurite_types],
#                           index = cell_IDs)
# bifurcation_ratio = pd.DataFrame(columns = [f'bifurcation_ratio-{ntype}' for ntype in neurite_types],
#                                   index = cell_IDs)

# # rename index columns
# for df in [n_primary, n_terminal, bifurcation_ratio]:
#     df.index.name = 'cell_ID'

# # iterate through cells
# for cell_ID in cell_IDs:
    
#     # get cell coordinates
#     cell_pricoordinates = coordinates_dict[cell_ID]['pri_coor']
#     cell_endcoordinates = coordinates_dict[cell_ID]['end_coor']
    
#     for ntype in neurite_types:
        
#         # limit coordinates to specific neurite type
#         if ntype == 'neurites':
#             cell_pricoordinates_pertype = cell_pricoordinates
#             cell_endcoordinates_pertype = cell_endcoordinates
            
#         elif ntype == 'dendrites':
#             cell_pricoordinates_pertype = cell_pricoordinates[cell_pricoordinates['path_label'] != 'axons']
#             cell_endcoordinates_pertype = cell_endcoordinates[cell_endcoordinates['path_label'] != 'axons']
            
#         else:
#             cell_pricoordinates_pertype = cell_pricoordinates[cell_pricoordinates['path_label'] == ntype]
#             cell_endcoordinates_pertype = cell_endcoordinates[cell_endcoordinates['path_label'] == ntype]       
        
#         # get number of points
#         n_primary_pertype = cell_pricoordinates_pertype.shape[0]
#         n_terminal_pertype = cell_endcoordinates_pertype.shape[0]
        
#         # write to dataframe
#         n_primary.at[cell_ID, f'n_primary-{ntype}'] = n_primary_pertype
#         n_terminal.at[cell_ID, f'n_terminal-{ntype}'] = n_terminal_pertype
  
#         # bifurcation ratio
#         ### TODO: more than one origin
        
#         # check if ntype exsists
#         if (ntype in cell_endcoordinates['path_label'].to_list()) or (ntype in ['neurites', 'dendrites']):
#             bifurcation_ratio_pertype = n_terminal_pertype / n_primary_pertype
        
#         # check different possible configurations (axon)
#         elif ntype == 'axons':
            
#             # condition for dendritic axons
#             if n_primary_pertype == 0 and n_terminal_pertype > 0:     
#                 bifurcation_ratio_pertype = n_terminal_pertype / 1
            
#             # condition for somatic axons
#             elif n_primary_pertype == 1 and n_terminal_pertype > 0:     
#                 bifurcation_ratio_pertype = n_terminal_pertype / n_primary_pertype
                
#             else:
#                 raise ValueError('miscalculation of axons')
                
#         else:
#             bifurcation_ratio_pertype = np.nan
            
#         # write to dataframe
#         bifurcation_ratio.at[cell_ID, f'bifurcation_ratio-{ntype}'] = bifurcation_ratio_pertype


# # %% plot primary & terminal points & bifurcation ratio

# if vplots:   
#     # load plotting function    
#     from cellmorphology.AMC_analysis.plot_AMC_cellcoordinates_analysis import plot_endpoints
    
#     # plot
#     for cell_ID in tqdm(cell_IDs):  
#         plot_endpoints(cell_ID = cell_ID, 
#                         cell_coordinates = coordinates_dict[cell_ID],
#                         n_primary = n_primary,
#                         n_terminal = n_terminal,
#                         bifurcation_ratio = bifurcation_ratio)

# # %%

# cell_coordinates = coordinates_dict[cell_ID]

# # init plotting
# from cellmorphology.cellmorph_functions.initialize_AMC_cellmorph_plotting import *

# fig, ax = plt.subplots(nrows = 1,
#                         ncols = 1,
#                         layout = 'constrained',
#                         figsize = get_figure_size(width = 150, height = 100),
#                         dpi = 300)

# # set figure title
# fig.suptitle(f'{cell_ID} primary and terminal points', 
#               fontsize = 9)

# # set aspect ration of plot
# ax.set_aspect(1)

# # plot all cell coordinates
# ax.scatter(x = cell_coordinates['all_coor'].loc[:, 'X'],
#             y = cell_coordinates['all_coor'].loc[:, 'Y'],
#             color = 'gray',
#             s = 0.25)

# # plot primary points (end points of primary paths)
# for path_i in cell_coordinates['pri_coor'].index.to_list():
#     ax.scatter(x = cell_coordinates['pri_coor'].at[path_i, 'X'],
#                 y = cell_coordinates['pri_coor'].at[path_i, 'Y'],
#                 color = neurite_color_dict[cell_coordinates['pri_coor'].at[path_i, 'path_label']],
#                 s = 25,
#                 marker = 'x',
#                 linewidths=0.5)

# # plot terminal points
# for path_i in cell_coordinates['end_coor'].index.to_list():
#     ax.scatter(x = cell_coordinates['end_coor'].at[path_i, 'X'],
#                 y = cell_coordinates['end_coor'].at[path_i, 'Y'],
#                 color = neurite_color_dict[cell_coordinates['end_coor'].at[path_i, 'path_label']],
#                 s = 5)

# # plot soma on top
# ax.scatter(x = soma_coordinates.at[0, 'X'],
#             y = soma_coordinates.at[0, 'Y'],
#             color = neurite_color_dict[soma_coordinates.at[0, 'path_label']])


# # plane label
# ax.text(x = 10, 
#         y = 10, 
#         s = 'XY', 
#         ha = 'left', 
#         va = 'top', 
#         fontsize = 9)

# # edit axes
# ax.set_ylim([0, max_fov_xy])
# ax.set_ylabel('Height [µm]')
# ax.set_yticks(ticks = np.arange(0, max_fov_xy, 200))
# ax.set_yticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)

# ax.set_xlim([0, max_fov_xy])
# ax.set_xlabel('Width [µm]')
# ax.set_xticks(ticks = np.arange(0, max_fov_xy, 200))
# ax.set_xticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)

# # invert y axis
# ax.invert_yaxis()

# # soma inset

# # inset marker
# box_xmin   = cell_coordinates['soma_coor'].at[0,'X']-30
# box_width  = 60 
# box_ymin   = cell_coordinates['soma_coor'].at[0,'Y']-30
# box_height = 60
   
# # add rectangle marker
# ax.add_patch(Rectangle(xy = (box_xmin, box_ymin), 
#                         width = box_width, 
#                         height = box_height,
#                         fill = False,
#                         color = primecolor,
#                         linestyle = '--',
#                         lw = 0.5,
#                         alpha = 0.5))

# ## ([left, bottom, width, height]), percentages
# ax_inset = ax.inset_axes([1.05, 0.63, 0.35, 0.35],
#                           xlim=(box_xmin, box_xmin+box_width), 
#                           ylim=(box_ymin, box_ymin+box_height), 
#                           xticklabels=[], 
#                           yticklabels=[])

# # edit linewidth of inset axis and its ticks
# [ax_inset.spines[spine].set_linewidth(0.5) for spine in ['left', 'bottom']]
# ax_inset.tick_params(width=0.5)
# ax_inset.tick_params(which = 'minor', width=0.25)


# # plot inset
# # plot all cell coordinates
# ax_inset.scatter(x = cell_coordinates['all_coor'].loc[:, 'X'],
#                   y = cell_coordinates['all_coor'].loc[:, 'Y'],
#                   color = 'gray',
#                   s = 0.25)

# # plot primary points (end points of primary paths)
# for path_i in cell_coordinates['pri_coor'].index.to_list():
#     ax_inset.scatter(x = cell_coordinates['pri_coor'].at[path_i, 'X'],
#                       y = cell_coordinates['pri_coor'].at[path_i, 'Y'],
#                     color = neurite_color_dict[cell_coordinates['pri_coor'].at[path_i, 'path_label']],
#                     s = 25,
#                     marker = 'x',
#                     linewidths=0.5)

# # plot terminal points
# for path_i in cell_coordinates['end_coor'].index.to_list():
#     ax_inset.scatter(x = cell_coordinates['end_coor'].at[path_i, 'X'],
#                       y = cell_coordinates['end_coor'].at[path_i, 'Y'],
#                   color = neurite_color_dict[cell_coordinates['end_coor'].at[path_i, 'path_label']],
#                   s = 5)
    
# # plot soma on top
# ax_inset.scatter(x = soma_coordinates.at[0, 'X'],
#                   y = soma_coordinates.at[0, 'Y'],
#                   color = neurite_color_dict[soma_coordinates.at[0, 'path_label']])
   
# # x
# ax_inset.set_xticks(ticks = np.arange(0, max_fov_xy, 200), labels = [])
# ax_inset.set_xticks(ticks = np.arange(0, max_fov_xy, 25), labels = [], minor = True)
# ax_inset.set_xlim([box_xmin, box_xmin + box_width])

# # y
# ax_inset.set_yticks(ticks = np.arange(0, max_fov_xy, 200), labels = [])
# ax_inset.set_yticks(ticks = np.arange(0, max_fov_xy, 25), labels = [], minor = True)
# ax_inset.set_ylim([box_ymin, box_ymin + box_height])
# ax_inset.invert_yaxis()

# # # update measurements label      
# # m_label = 'n_primary [#] - n_terminal [#] - bifurcation_ratio'

# # for ntype in neurite_types:
# #     m_label = m_label + f'\n{ntype}: {"{:.2f}".format(n_primary.at[cell_ID, f"n_primary-{ntype}"])} - {"{:.2f}".format(n_terminal.at[cell_ID, f"n_terminal-{ntype}"])} - {"{:.2f}".format(bifurcation_ratio.at[cell_ID, f"bifurcation_ratio-{ntype}"])}'

# # # add measurements to plot
# # ax.text(x = 827, y = 580, 
# #         s = m_label, 
# #         ha = 'right', 
# #         va = 'bottom',
# #         size = 4)

# # create saving path and save
# from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_analysis_dir

# # save figure
# save_figures(fig, 
#               f'{cell_ID}-primary_terminal_bifurcation', 
#               join(AMCs_analysis_dir, 'primary_terminal_bifurcation_plots'), 
#               darkmode_bool, 
#               figure_format='png')
    

# # display plot
# plt.show()



