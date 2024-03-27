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
cell_in_files = 1
all_coor_filename = onlyfiles_allcoor[cell_in_files]

# read all coordinates file
all_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, all_coor_filename)) 
end_coordinates = pd.read_csv(join(cell_morph_traces_coordinates_dir, onlyfiles_endcoor[cell_in_files])) 

# clean up dataframes
def clean_coordinates_df(coordinates_df):
    coordinates_df['path_ID'] = coordinates_df['On Path'].str.replace('Path ', '').astype(int)
    coordinates_df.drop(columns = ['On Path'], inplace = True)

clean_coordinates_df(all_coordinates)

end_coordinates['On Path'] = end_coordinates['On Path'].replace('Path (1) [Single Point]', '1')
clean_coordinates_df(end_coordinates)


# get total number of paths
n_paths = all_coordinates['path_ID'].max()



# %%

plt.scatter(all_coordinates['X'], all_coordinates['Y'], s = 0.5)
plt.scatter(end_coordinates['X'], end_coordinates['Y'], s = 0.75, color = 'r')
plt.scatter(end_coordinates['X'][0], end_coordinates['Y'][0], s = 2, color = 'm')
plt.gca().invert_yaxis()
plt.xlim([0, 590.76])
plt.ylim([590.76, 0])


# get soma coordinates
soma_coordinates = end_coordinates[end_coordinates['path_ID'] == 1]

# define vector to compare to
reference_vector = pd.Series({'X' : 1., 'Y' : 0.})

# define dataframe for all vectors to endpoints
terminal_branches_df = pd.DataFrame(soma_coordinates.rename(columns = {'X': 'x_end', 'Y': 'y_end', 'Z': 'z_end'}).set_index('path_ID'))

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


# %%


# skip first to avoid the soma
for idx_end_point in end_coordinates.index[2:3]:
    
    path_ID_end_point = end_coordinates.at[idx_end_point, 'path_ID']
    
    all_coordinates_path = all_coordinates[all_coordinates['path_ID'] == path_ID_end_point]
    
    # reset the index to start from 0
    all_coordinates_path.reset_index(drop=True, inplace=True)
    
    # get coordinates of first node
    first_coordinates_path = all_coordinates_path.iloc[0, :]
    
    print(idx_end_point, path_ID_end_point)

    for path_idx in range()

# TODO:
    # polar plot with length-colorcoded lines
    # polar plot with histogram
        # bin size (quaters or )
    # histogram
    # calc length of terminal branches
        # get intersect points
    
    
    




