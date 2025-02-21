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
cell_IDs = ['Exp-160']

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

# initialize dataframes
n_primary = pd.DataFrame(columns = [f'n_primary-{ntype}' for ntype in neurite_types],
                         index = cell_IDs)
n_terminal = pd.DataFrame(columns = [f'n_terminal-{ntype}' for ntype in neurite_types],
                          index = cell_IDs)
bifurcation_ratio = pd.DataFrame(columns = [f'bifurcation_ratio-{ntype}' for ntype in neurite_types],
                                 index = cell_IDs)

# rename index columns
for df in [n_primary, n_terminal, bifurcation_ratio]:
    df.index.name = 'cell_ID'

for cell_ID in cell_IDs:
    
    # get cell coordinates
    cell_pricoordinates = coordinates_dict[cell_ID]['pri_coor']
    cell_endcoordinates = coordinates_dict[cell_ID]['end_coor']
    
    for ntype in neurite_types:
        # print(ntype)
        
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
        
        # write to dataframe
        n_primary.at[cell_ID, f'n_primary-{ntype}'] = cell_pricoordinates_pertype.shape[0]
        n_terminal.at[cell_ID, f'n_terminal-{ntype}'] = cell_endcoordinates_pertype.shape[0]
        
        # DENDRITIC AXONS!!!!
        
        # bifurcation ratio
        
        # check if ntype exsists
        if (ntype in cell_endcoordinates['path_label'].to_list()) or (ntype in ['neurites', 'dendrites', 'axons']):
            bifurcation_ratio.at[cell_ID, f'bifurcation_ratio-{ntype}'] = cell_endcoordinates_pertype.shape[0] / cell_pricoordinates_pertype.shape[0]
        else:
            bifurcation_ratio.at[cell_ID, f'bifurcation_ratio-{ntype}'] = np.nan

# %%

# init plotting
from cellmorphology.cellmorph_functions.initialize_AMC_cellmorph_plotting import *


# cell_ID
cell_coordinates = coordinates_dict[cell_ID]

fig, ax = plt.subplots(nrows = 1,
                       ncols = 1,
                       layout = 'constrained',
                       figsize = get_figure_size(width = 150, height = 100),
                       dpi = 300)

# set figure title
fig.suptitle(f'{cell_ID} primary and terminal points', 
             fontsize = 9)

# set aspect ration of plot
ax.set_aspect(1)

# plot all cell coordinates
ax.scatter(x = cell_coordinates['all_coor'].loc[:, 'X'],
           y = cell_coordinates['all_coor'].loc[:, 'Y'],
           color = 'gray',
           s = 0.25)

# plot primary points (end points of primary paths)
for path_i in cell_coordinates['pri_coor'].index.to_list():
    ax.scatter(x = cell_coordinates['pri_coor'].at[path_i, 'X'],
               y = cell_coordinates['pri_coor'].at[path_i, 'Y'],
               color = neurite_color_dict[cell_coordinates['pri_coor'].at[path_i, 'path_label']],
               s = 25,
               marker = 'x',
               linewidths=0.5)

# plot terminal points
for path_i in cell_coordinates['end_coor'].index.to_list():
    ax.scatter(x = cell_coordinates['end_coor'].at[path_i, 'X'],
               y = cell_coordinates['end_coor'].at[path_i, 'Y'],
               color = neurite_color_dict[cell_coordinates['end_coor'].at[path_i, 'path_label']],
               s = 5)

# plot soma on top
ax.scatter(x = soma_coordinates.at[0, 'X'],
           y = soma_coordinates.at[0, 'Y'],
           color = neurite_color_dict[soma_coordinates.at[0, 'path_label']])

max_fov_xy = 590.76

# plane label
ax.text(x = 10, 
        y = 10, 
        s = 'XY', 
        ha = 'left', 
        va = 'top', 
        fontsize = 9)

# edit axes
ax.set_ylim([0, max_fov_xy])
ax.set_ylabel('Height [µm]')
ax.set_yticks(ticks = np.arange(0, max_fov_xy, 200))
ax.set_yticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)

ax.set_xlim([0, max_fov_xy])
ax.set_xlabel('Width [µm]')
ax.set_xticks(ticks = np.arange(0, max_fov_xy, 200))
ax.set_xticks(ticks = np.arange(0, max_fov_xy, 25), minor = True)

# invert y axis
ax.invert_yaxis()

# soma inset

# inset marker
box_xmin   = cell_coordinates['soma_coor'].at[0,'X']-40
box_width  = 80 
box_ymin   = cell_coordinates['soma_coor'].at[0,'Y']-40
box_height = 80
   
# add rectangle marker
ax.add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                       width = box_width, 
                       height = box_height,
                       fill = False,
                       color = primecolor,
                       linestyle = '--',
                       lw = 0.5,
                       alpha = 0.5))

## ([left, bottom, width, height]), percentages
ax_inset = ax.inset_axes([1.05, 0.63, 0.35, 0.35],
                          xlim=(box_xmin, box_xmin+box_width), 
                          ylim=(box_ymin, box_ymin+box_height), 
                          xticklabels=[], 
                          yticklabels=[])

# edit linewidth of inset axis and its ticks
[ax_inset.spines[spine].set_linewidth(0.5) for spine in ['left', 'bottom']]
ax_inset.tick_params(width=0.5)
ax_inset.tick_params(which = 'minor', width=0.25)


# plot inset
# plot all cell coordinates
ax_inset.scatter(x = cell_coordinates['all_coor'].loc[:, 'X'],
                 y = cell_coordinates['all_coor'].loc[:, 'Y'],
                 color = 'gray',
                 s = 0.25)

# plot primary points (end points of primary paths)
for path_i in cell_coordinates['pri_coor'].index.to_list():
    ax_inset.scatter(x = cell_coordinates['pri_coor'].at[path_i, 'X'],
                     y = cell_coordinates['pri_coor'].at[path_i, 'Y'],
                   color = neurite_color_dict[cell_coordinates['pri_coor'].at[path_i, 'path_label']],
                   s = 25,
                   marker = 'x',
                   linewidths=0.5)

# plot terminal points
for path_i in cell_coordinates['end_coor'].index.to_list():
    ax_inset.scatter(x = cell_coordinates['end_coor'].at[path_i, 'X'],
                     y = cell_coordinates['end_coor'].at[path_i, 'Y'],
                  color = neurite_color_dict[cell_coordinates['end_coor'].at[path_i, 'path_label']],
                  s = 5)
    
# plot soma on top
ax_inset.scatter(x = soma_coordinates.at[0, 'X'],
                 y = soma_coordinates.at[0, 'Y'],
                 color = neurite_color_dict[soma_coordinates.at[0, 'path_label']])
   
# x
ax_inset.set_xticks(ticks = np.arange(0, max_fov_xy, 200), labels = [])
ax_inset.set_xticks(ticks = np.arange(0, max_fov_xy, 25), labels = [], minor = True)
ax_inset.set_xlim([box_xmin, box_xmin + box_width])

# y
ax_inset.set_yticks(ticks = np.arange(0, max_fov_xy, 200), labels = [])
ax_inset.set_yticks(ticks = np.arange(0, max_fov_xy, 25), labels = [], minor = True)
ax_inset.set_ylim([box_ymin, box_ymin + box_height])
ax_inset.invert_yaxis()

# update measurements label      
m_label = 'n_primary [#] - n_terminal [#] - bifurcation_ratio'

for ntype in neurite_types:
    m_label = m_label + f'\n{ntype}: {"{:.2f}".format(n_primary.at[cell_ID, f"n_primary-{ntype}"])} - {"{:.2f}".format(n_terminal.at[cell_ID, f"n_terminal-{ntype}"])} - {"{:.2f}".format(bifurcation_ratio.at[cell_ID, f"bifurcation_ratio-{ntype}"])}'

# add measurements to plot
ax.text(x = 580, y = 580, 
        s = m_label, 
        ha = 'right', 
        va = 'bottom',
        size = 4)

# display plot
plt.show()