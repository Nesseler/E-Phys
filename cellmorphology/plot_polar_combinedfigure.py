# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 14:26:16 2024

@author: nesseler
"""

import pandas as pd
from os import mkdir
from os.path import join, exists
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np



from parameters.directories_win import table_file, cell_morph_descrip_dir, cell_morph_traces_coordinates_dir

from functions.functions_plotting import get_figure_size, get_colors
from cellmorphology.functions_cellmorph import clean_OnPath_column_to_path_ID_n_label

# settings
# get colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)


# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')


# %%

cell_ID = 'E-137'


# load dataframes

### terminal branches measurements
# diretory
terminal_branch_measurements_path = join(cell_morph_descrip_dir, 'terminal_branches_measurements', f'{cell_ID}-terminal_branches.xlsx')

# dataframe
terminal_branches_df = pd.read_excel(terminal_branch_measurements_path, index_col = 'path_ID')


### coordinates
# all coordinates
all_coordinates_path = join(cell_morph_traces_coordinates_dir, f'{cell_ID}-all_coordinates.csv')
cell_allcoordinates = pd.read_csv(all_coordinates_path)
clean_OnPath_column_to_path_ID_n_label(cell_allcoordinates)


# end / last / terminal coordinates
last_coordinates_path = join(cell_morph_traces_coordinates_dir, f'{cell_ID}-last_coordinates.csv')
cell_endcoordinates = pd.read_csv(last_coordinates_path)
clean_OnPath_column_to_path_ID_n_label(cell_endcoordinates)


# check if x coordinates need to be flipped
if MetaData.at[cell_ID, 'to_be_x_flipped']:
    print('flipping x coordinates')
    cell_allcoordinates['X'] = 590.76 - cell_allcoordinates['X']
    cell_endcoordinates['X'] = 590.76 - cell_endcoordinates['X']


# get soma coordinates
soma_coordinates = cell_endcoordinates[cell_endcoordinates['path_ID'] == 1]


# get all path_IDs
path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()


# %% combined figure for one cell

fig_all, ax_all = plt.subplots(nrows = 2,
                                ncols = 2,
                                layout = 'constrained',
                                height_ratios = [1, 1],
                                width_ratios = [1,1],
                                figsize = get_figure_size(width = 185.5))

# set figure title
fig_all.suptitle(cell_ID)




### subplot 1: XY coordinates ###

# loop through path
for path_ID in path_IDs:
    
    # get all coordinates
    path_all_coordinates = cell_allcoordinates[cell_allcoordinates['path_ID'] == path_ID]
    path_end_coordinates = cell_endcoordinates[cell_endcoordinates['path_ID'] == path_ID]           

    # get label of current path
    cur_path_label = path_all_coordinates['path_label'].iloc[0]
    
    # set colors for paths
    path_color_dict = {'dendrite' : {'path_color' : colors_dict['primecolor'], 'end_color' : 'r'},
                        'axon' :     {'path_color' : 'gray', 'end_color' : 'lightcoral'},
                        'soma' :     {'path_color' : colors_dict['color2'], 'end_color' : colors_dict['color2']}}

    path_color = path_color_dict[cur_path_label]['path_color']
    end_color = path_color_dict[cur_path_label]['end_color']


    
    ax_all.flat[0].scatter(path_all_coordinates['X'], path_all_coordinates['Y'],
                            s = 0.5, c = path_color, label = 'cell')
        
    ax_all.flat[0].scatter(path_end_coordinates['X'], path_end_coordinates['Y'], 
                            s = 0.75, c = end_color, label = 'end points')
    
ax_all.flat[0].scatter(soma_coordinates['X'], soma_coordinates['Y'], s = 10, color = colors_dict['color2'], label = 'soma')
ax_all.flat[0].text(x = 10, y = 10, s = 'XY', ha = 'left', va = 'top')
ax_all.flat[0].set_xlim([0, 590.76])
ax_all.flat[0].set_ylim([590.76, 0])

ax_all.flat[0].set_xticks(np.arange(0, 590, 200))
ax_all.flat[0].set_yticks(np.arange(0, 590, 200))



# define function to change projection type of subplot specific subplot
def change_projection(fig, axs, ax_tochange, projection = 'polar'):

    rows, cols, start, stop = ax_tochange.get_subplotspec().get_geometry()

    axs.flat[start].remove()
    axs.flat[start] = fig.add_subplot(rows, cols, start+1, projection=projection)

# change projection of subplots
[change_projection(fig_all, ax_all, ax, 'polar') for ax in ax_all.flat[1:]]



### all polar plot ### 
### color code for length of branches ###
# initialise color code
norm_min = 0
norm_max = 1000
cmap_str = 'viridis'
norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# colorbar
fig_all.colorbar(cmap, ax = ax_all.flat[3], label = 'Terminal branch length [µm]')


# initilize polar histogram
resul_n_bins = 8
binsize = (2 * np.pi) / resul_n_bins


def calc_polar_histo_binangles(n_bins):
    
    # step size is set bin number of resulting bins and 2*pi
    bin_stepsize = (2 * np.pi) / n_bins
    
    # start point: half of step size
    # because of rotated polar bins
    bin_start = bin_stepsize / 2
    
    # calc bin borders
    bin_angles = np.arange(bin_start, np.pi * 2, bin_stepsize)

    # roll bin borders
    bin_angles = np.roll(bin_angles, 1)    

    return bin_angles

bins_angles = calc_polar_histo_binangles(resul_n_bins)


            
                
                
                
                
                


def plot_colorcoded_polar(polar_occurances_df, ax):

    # define array with number of previouse numbers of branches in bin
    bottom = [0] * resul_n_bins
    
    # skip (drop) path 1, i.e. soma
    if 1 in polar_occurances_df.index.to_list():
        branch_idc = polar_occurances_df.drop(index = 1).index.to_list()
    else:
        branch_idc = polar_occurances_df.index.to_list()
    
    # loop through all branches to assign specific color for length            
    for branch_idx in branch_idc:
    
        # get angles of branches
        branch_length = polar_occurances_df.at[branch_idx, "length"]
        branch_bin = polar_occurances_df.at[branch_idx, "bin_id"].astype(int)
        
        # create empty bins and assign branch to bin
        hist_angles_occu = [0] * resul_n_bins
        hist_angles_occu[branch_bin] = 1
        
        # plot histogram as barplot
        ax.bar(bins_angles, hist_angles_occu, bottom = bottom,
                width = binsize, 
                align = 'edge',
                edgecolor = 'none',
                color = cmap.to_rgba(branch_length))
            
        # add to bottom list for next step
        bottom = np.add(bottom, hist_angles_occu)


# plot all
plot_colorcoded_polar(terminal_branches_df, ax_all.flat[1])

# plot dendrites
plot_colorcoded_polar(terminal_branches_df[terminal_branches_df['path_label'] == 'dendrite'], ax_all.flat[2])        

# plot dendrites
plot_colorcoded_polar(terminal_branches_df[terminal_branches_df['path_label'] == 'axon'], ax_all.flat[3])           

# axis for polar plots

for ax in ax_all.flat[1:]:
    
    # x axis
    ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
    ax.set_xticklabels(['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv'])
    
    # y axis          
    ax.set_yticks(ticks = np.arange(0, 15 + 1, 5))
    # ax.set_yticks(np.arange(0, round_to_base(max(bottom)+1, 5)+1, 1), minor = True)
    
    ax.grid(True, alpha = 0.5)
    

plt.show()