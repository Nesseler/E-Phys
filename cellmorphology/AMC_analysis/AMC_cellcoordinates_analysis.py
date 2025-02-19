# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 14:57:41 2025

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# get MetaData
from cellmorphology.AMC_analysis.AMC_analysis_import import get_cells_list
MetaData = get_cells_list()

# get cell_IDs to be analyzed
cell_IDs = MetaData.query('coordinates == "Yes" & paths_checked == "Yes"').index.to_list()

# set neurite types
neurite_types = ['neurites', 
                 'dendrites', 'glomerular_dendrites', 'basal_dendrites', 'LOT_dendrites', 'undefined_dendrites', 
                 'axons']

# field of view dimension
max_fov_xy = 590.76
max_fov_z = 300


# %% load coordinates

# all coordinates dict
coordinates_dict = dict.fromkeys(['all_coor', 'end_coor', 'soma_coor', 'path_IDs'])

# get directory of cell coordinates
from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_coordinates_dir

# import function that cleans up imported coordinates table
from cellmorphology.AMC_analysis.AMC_functions import clean_OnPath_column_to_path_ID_n_label_AMCs


# for cell_ID in cell_IDs:
    
cell_ID = 'Exp-162'
        
### coordinates
# all coordinates
all_coordinates_path = join(AMCs_coordinates_dir, f'{cell_ID}-all_coordinates.csv')
cell_allcoordinates = pd.read_csv(all_coordinates_path)
clean_OnPath_column_to_path_ID_n_label_AMCs(cell_allcoordinates)

# end / last / terminal coordinates
last_coordinates_path = join(AMCs_coordinates_dir, f'{cell_ID}-terminal_last_coordinates.csv')
cell_endcoordinates = pd.read_csv(last_coordinates_path)
clean_OnPath_column_to_path_ID_n_label_AMCs(cell_endcoordinates) 

# get soma coordinates
soma_coordinates = cell_endcoordinates[cell_endcoordinates['path_ID'] == 1]
 
# get all path_IDs
path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()

# set dict for all coordinates
cell_coordintes_dict = {'all_coor' : cell_allcoordinates,
                        'end_coor' : cell_endcoordinates,
                        'soma_coor': soma_coordinates,
                        'path_IDs' : path_IDs}
#                         # 'terminal_branches' : terminal_branches_df}

# # add to dict 
coordinates_dict[cell_ID] = cell_coordintes_dict


# %%

# initialize plotting packages
from cellmorphology.cellmorph_functions.initialize_AMC_cellmorph_plotting import *

scatter_plot_dict = {'s' : 0.25}


print('plotting ...')

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
fig.suptitle(f'{cell_ID} cell coordinates', 
             fontsize = 9)
    
# set coordinates data dict
cell_coordinates = coordinates_dict[cell_ID]

# get all path_IDs
cell_path_IDs = cell_coordinates['path_IDs']

# loop through all paths of cell
for path_ID in cell_path_IDs:
    
    # get path coordinates
    path_all_coordinates = cell_coordinates['all_coor'][cell_coordinates['all_coor']['path_ID'] == path_ID]
    # path_end_coordinates = cell_coordinates['end_coor'][cell_coordinates['end_coor']['path_ID'] == path_ID]           

    # get label of current path
    cur_path_label = path_all_coordinates['path_label'].iloc[0]
       
    # scatter plots
    for ax, dim1, dim2 in zip(axs[:3], ['X', 'Z', 'X'], ['Y', 'Y', 'Z']):
        
        # all coordinates
        ax.scatter(x = path_all_coordinates[dim1],
                   y = path_all_coordinates[dim2],
                   color = neurite_color_dict[cur_path_label],
                   label = cur_path_label,
                   **scatter_plot_dict)

# TODO: zorders, soma on top

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
axs[0].tick_params(axis = 'x', size = 0)
axs[0].tick_params(axis = 'x', which = 'minor', size = 0)
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
axs[1].tick_params(axis = 'y', size = 0)
axs[1].tick_params(axis = 'y', which = 'minor', size = 0)
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

# legend
h, l = axs[0].get_legend_handles_labels()

# get unique lables
unique_label = list(set(l))

# get and sort indices of unique lables
unique_label_indices = sorted([l.index(unique_label[u_i]) for u_i, u_l in enumerate(unique_label)])

# reassign handels
unique_handles = [h[u_li] for u_li in unique_label_indices]
unique_label = [l[u_li] for u_li in unique_label_indices]

axs[3].legend(unique_handles, unique_label, 
              title = 'neurite type',
              title_fontsize = 6,
              frameon = False, 
              ncol = 1, 
              loc = 'center',
              fontsize = 6)

axs[3].set_xticklabels(labels = [])
axs[3].tick_params(axis = 'y', size = 0)
axs[3].tick_params(axis = 'y', which = 'minor', size = 0)
axs[3].tick_params(axis = 'x', size = 0)
axs[3].tick_params(axis = 'x', which = 'minor', size = 0)

# remove spines
[axs[3].spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom', 'left']]

# align labels
fig.align_labels()

# display figure
plt.show()


