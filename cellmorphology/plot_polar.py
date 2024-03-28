# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 13:53:45 2024

@author: nesseler
"""

import pandas as pd
from os.path import join
import matplotlib.pyplot as plt
import numpy as np
import math

from parameters.directories_win import cell_morph_traces_coordinates_dir
from getter.get_onlyfiles_list import get_onlyfiles_list

# get onlyfiles list
onlyfiles = get_onlyfiles_list(cell_morph_traces_coordinates_dir)

# extract only filenames with all coordinates
onlyfiles_allcoor = [f for f in onlyfiles if 'all_coordinates' in f]
onlyfiles_endcoor = [f for f in onlyfiles if 'last_coordinates' in f]

# test 
cell_in_files = 0
all_coor_filename = onlyfiles_allcoor[cell_in_files]

# read all coordinates file
all_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, all_coor_filename)) 
end_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, onlyfiles_endcoor[cell_in_files])) 

# clean up dataframes
def clean_coordinates_df(coordinates_df):
    coordinates_df['path_ID'] = [int(s.replace("Path ", "").replace("(", "").replace(")", "").replace('1 [Single Point]', '1')) for s in coordinates_df['On Path']]
    coordinates_df.drop(columns = ['On Path'], inplace = True)

clean_coordinates_df(all_coordinates)

# end_coordinates['On Path'] = end_coordinates['On Path'].replace('Path (1) [Single Point]', '1')
clean_coordinates_df(end_coordinates)


# get total number of paths
n_paths = all_coordinates['path_ID'].max()



# %%

plt.scatter(all_coordinates['X'], all_coordinates['Y'], s = 0.5)
plt.scatter(end_coordinates['X'], end_coordinates['Y'], s = 0.75, color = 'r')
plt.scatter(end_coordinates['X'][0], end_coordinates['Y'][0], s = 2, color = 'm')
plt.gca().invert_yaxis()
plt.gca().set_box_aspect(1)
plt.xlim([0, 590.76])
plt.ylim([590.76, 0])


# get soma coordinates
soma_coordinates = end_coordinates[end_coordinates['path_ID'] == 1]

# define vector to compare to
reference_vector = pd.Series({'X' : 1., 'Y' : 0.})

# define dataframe for all vectors to endpoints
terminal_branches_df = pd.DataFrame(soma_coordinates.rename(columns = {'X': 'x_end', 'Y': 'y_end', 'Z': 'z_end'}).set_index('path_ID'))

# reindex all_coordinates
# all_coordinates.index += 1

# append soma coordinates to all coordinates
# all_coordinates = pd.concat([all_coordinates, soma_coordinates], axis = 0)

ls_deg = []


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
    
    plt.plot(vector_coordinates['X'], vector_coordinates['Y'], 'w', alpha = 0.5)
    
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


# %%

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

for i in terminal_branches_df.index:
    
    ax.plot([terminal_branches_df.at[i, 'angle_rad'], terminal_branches_df.at[i, 'angle_rad']], [0.1, 1], 'w', alpha = 0.5)
    ax.set_rticks([1])
    ax.set_xticks(np.arange(0, np.pi*2, np.pi/4))
    print(i, terminal_branches_df.at[i, 'angle_deg'])
    
ax.grid(False)

plt.show()

# %%


# skip first to avoid the soma
for idx_end_point in end_coordinates.index[2:3]:
    
    path_ID_end_point = end_coordinates.at[idx_end_point, 'path_ID']
    
    all_coordinates_path = all_coordinates[all_coordinates['path_ID'] == path_ID_end_point]
    
    # reset the index to start from 0
    all_coordinates_path.reset_index(drop=True, inplace=True)
    
    # get coordinates of first node
    first_node = all_coordinates_path.iloc[0, :]
    
    test_all_coor_parent = all_coordinates[all_coordinates['path_ID'] == 2]
    
    test_all_coor_non_parent = all_coordinates[all_coordinates['path_ID'] == 25]

    for path_idx in np.arange(1, n_paths + 1, 1):
        
        
        print(path_idx, len(all_coordinates[all_coordinates['path_ID'] == path_idx]))
        
        break
    
# %%


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

## list of all paths except current terminal path
path_IDs = np.arange(1, n_paths + 1, 1)

terminal_path_ID = 26

terminal_path_IDs = end_coordinates['path_ID'][end_coordinates['path_ID'] != 1]

for terminal_path_ID in terminal_path_IDs:
    
    print(terminal_path_ID)
    
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
        
        # set path_ID to parent for next step in while loop
        path_ID = parent_path_ID
    
    
    test = branch_coor_terminal_to_soma
    
    # test = test[test.index > 5000]
    
    # test['dist'] = test[['X', 'Y', 'Z']].apply(lambda row: np.linalg.norm((row.X, row.Y, row.Z)), axis=1)
    # test.sort_values('dist', ignore_index=True, inplace=True)
    plt.title(terminal_path_ID)
    plt.plot(test['X'], test['Y'])#, s = 0.5, alpha = 0.1)
    
    plt.scatter(end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['X'], 
                end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]['Y'], 
                s = 0.75, color = 'r')
    
    plt.scatter(end_coordinates['X'][0], end_coordinates['Y'][0], s = 2, color = 'm')
    
    end_point_coordinates = end_coordinates[end_coordinates['path_ID'] == terminal_path_ID]

    vector_coordinates = pd.concat([soma_coordinates, end_point_coordinates], axis = 0)
    vector_coordinates.drop(columns = ['path_ID'], inplace = True)
    
    plt.plot(vector_coordinates['X'], vector_coordinates['Y'], 'w', alpha = 0.5)
    
    plt.gca().invert_yaxis()
    plt.gca().set_box_aspect(1)
    plt.xlim([0, 590.76])
    plt.ylim([590.76, 0])
    plt.show()




# %%



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
    

    




