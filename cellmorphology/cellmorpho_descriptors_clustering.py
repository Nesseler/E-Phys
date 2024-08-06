# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:22:45 2024

@author: nesseler
"""

from cellmorphology.cellmorph_packages import mtl, plt, sbn, pd, np, join

from parameters.directories_win import table_file, cell_morph_descrip_dir

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# get cell_IDs:
cell_IDs = MetaData[MetaData['reconstructed'] == 1].index.to_list()

# set list of neurite types
neurite_types = ['neurites', 'dendrites', 'axons']


# %%

# create dataframe that contains all parameters
cellmorph_descriptors = pd.DataFrame(index = cell_IDs)
cellmorph_descriptors.index.name = 'cell_ID'


# %% load height & width

# load
height_width = pd.read_excel(join(cell_morph_descrip_dir, 'height_width.xlsx'), index_col = 'cell_IDs')

# get list of columns to use
height_width_cols = [col for col in height_width.columns.to_list() if 'width' in col and 'neurites' not in col or 'height' in col and 'neurites' not in col]

# concat
cellmorph_descriptors = pd.concat([cellmorph_descriptors, height_width[height_width_cols]], axis = 1)


# %% load number of points

# load
n_primary = pd.read_excel(join(cell_morph_descrip_dir, 'n_primary_points.xlsx'), index_col = 'cell_ID')
n_terminal = pd.read_excel(join(cell_morph_descrip_dir, 'n_terminal_points.xlsx'), index_col = 'cell_ID')
n_bifurcation_ratios = pd.read_excel(join(cell_morph_descrip_dir, 'bifurcation_ratios.xlsx'), index_col = 'cell_ID')

# concat
cellmorph_descriptors = pd.concat([cellmorph_descriptors, 
                                   n_primary[['dendritic_primaries', 'axonic_primaries']],
                                   n_terminal[['dendritic_terminals', 'axonic_terminals']],
                                   n_bifurcation_ratios[['bifurcation_ratio_dendritic', 'bifurcation_ratio_axonic']]]
                                  , axis = 1)


# %% load total cable lengths

# load
total_cable_length = pd.read_excel(join(cell_morph_descrip_dir, 'total_cable_length.xlsx'), index_col = 'cell_ID')

# filter out neurites
total_cable_length_cols = [col for col in total_cable_length.columns.to_list() if 'neurites' not in col]

# concat
cellmorph_descriptors = pd.concat([cellmorph_descriptors, total_cable_length[total_cable_length_cols]], axis = 1)


# %% load critical and enclosing radius and max number of intersections

# load
# sholl metrics
sholl_metrics = {'neurites': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_neurites.xlsx'), index_col='cell_ID'),
                 'dendrites': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_dendrites.xlsx'), index_col='cell_ID'),
                 'axons': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_axons.xlsx'), index_col='cell_ID')}

for neurite_type in ['dendrites', 'axons']:
    
    # get sholl metrics
    sholl_metrics_pertype = sholl_metrics[neurite_type]

    # rename columns with neurite_type
    sholl_metrics_pertype.rename(columns = {'critical_radius'   : f'{neurite_type}-critical_radius', 
                                            'enclosing_radius'  : f'{neurite_type}-enclosing_radius', 
                                            'max_intersections' : f'{neurite_type}-max_intersections'},
                                 inplace = True)

    # concatenate to descriptors
    cellmorph_descriptors = pd.concat([cellmorph_descriptors, sholl_metrics_pertype], axis = 1)
    

# %% load AIS data

# load
axon_data = pd.read_excel(join(cell_morph_descrip_dir, 'axon_data.xlsx'), index_col='cell_ID')

# define columns to add
axon_cols = ['distance to soma']
axon_cols_ext = ['length(µm)', 'distance to soma', 'length + distance', 'source']

# rename columns with neurite_type
axon_data.rename(columns = {'distance to soma'   : 'AIS distance to soma'},
                inplace = True)

# concat
cellmorph_descriptors = pd.concat([cellmorph_descriptors, axon_data['AIS distance to soma']], axis = 1)


# %% load spine categorization

# load
spines_df = pd.read_excel(join(cell_morph_descrip_dir, 'spines.xlsx'), index_col='cell_ID')


# %% normalise cell descriptors

# min-max normalize cellmorph matrix
cellmorph_descriptors_minmax = (cellmorph_descriptors - cellmorph_descriptors.min()) / (cellmorph_descriptors.max() - cellmorph_descriptors.min())

# z-score cellmorph matrix
cellmorph_descriptors_zscored = (cellmorph_descriptors - cellmorph_descriptors.mean()) / cellmorph_descriptors.std()


# %% sort dataframe

# cellmorph_descriptors_zscored.sort_values(['total_cable_length-axons'], inplace = True)


# %% initialize plotting

from functions.functions_plotting import save_figures, get_colors, get_figure_size
from cellmorphology.cellmorph_parameters import cell_coordinates_field_of_view

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})





fig_heat, ax_heat = plt.subplots(nrows = 1,
                                  ncols = 1,
                                  layout = 'constrained',
                                  figsize = get_figure_size(width = 150, height = 125),
                                  dpi = 600)


sbn.heatmap(cellmorph_descriptors_zscored,
            vmin = -4,
            vmax = 4,
            square = False, 
            ax = ax_heat, 
            cmap="bwr", 
            yticklabels=False,
            linewidth = 0) 

    
# show plot
plt.show()