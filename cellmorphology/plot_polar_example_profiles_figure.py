# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:28:16 2024

@author: nesseler
"""

from cellmorphology.cellmorph_packages import plt, mtl, sbn, pd, join, np

# script specific directories / parameters / functions
from parameters.directories_win import cell_morph_descrip_dir, table_file, cell_morph_plots_dir, cell_morph_traces_coordinates_dir

# load cellmorpho parameters
# from cellmorphology.cellmorph_parameters import orientation_labels
orientation_labels = ['p', '', 'd', '', 'a', '', 'v', '']

# load cellmorhp functions
from cellmorphology.functions_cellmorph import calc_polar_histo_binangles, clean_OnPath_column_to_path_ID_n_label

# calc bins angles with functions
bins_angles, bin_size = calc_polar_histo_binangles(n_bins=8)

# set cell_IDs for figure
cell_IDs_dict = {'MeA' : ['E-065', 'E-111'],
                 'BAOT': ['E-089', 'E-130']}

# get cell_IDs list
cell_IDs = [cell_ID for region in cell_IDs_dict.values() for cell_ID in region]

# %% load

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# load histogram data
polar_plot_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_occurrances.xlsx'), index_col = 'orientation_rad')
polar_plot_dendrites_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_dendrites_occurrances.xlsx'), index_col = 'orientation_rad')
polar_plot_axons_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_axons_occurrances.xlsx'), index_col = 'orientation_rad')

# all coordinates dict
coordinates_dict = {}

# load cell coordinates
for cell_ID in cell_IDs:
        
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
    last_coordinates_path = join(cell_morph_traces_coordinates_dir, f'{cell_ID}-terminal_last_coordinates.csv')
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
    
    # set dict for all coordinates
    cell_coordintes_dict = {'all_coor' : cell_allcoordinates,
                            'end_coor' : cell_endcoordinates,
                            'soma_coor': soma_coordinates,
                            'path_IDs' : path_IDs,
                            'terminal_branches' : terminal_branches_df}
    
    # add to dict 
    coordinates_dict[cell_ID] = cell_coordintes_dict


# %% norm

### norm to neurites
# neurites
pp_occu_neurites_normed_to_neurites = polar_plot_occurrances.div(polar_plot_occurrances.sum(), axis = 1)

# dendrites
pp_occu_dendrites_normed_to_neurites = polar_plot_dendrites_occurrances.div(polar_plot_occurrances.sum(), axis = 1)

# axons
pp_occu_axons_normed_to_neurites = polar_plot_axons_occurrances.div(polar_plot_occurrances.sum(), axis = 1)


### norm to type
# dendrites
pp_occu_dendrites_normed_to_type = polar_plot_dendrites_occurrances.div(polar_plot_dendrites_occurrances.sum(), axis = 1)

# axons
pp_occu_axons_normed_to_type = polar_plot_axons_occurrances.div(polar_plot_axons_occurrances.sum(), axis = 1)


# %% initialize plotting

from functions.functions_plotting import save_figures, get_colors, get_figure_size, set_font_sizes, change_projection
from cellmorphology.cellmorph_colors import neurite_color_dict
from cellmorphology.cellmorph_parameters import cell_coordinates_field_of_view
from cellmorphology.functions_cellmorph import plot_colorcoded_polar_normed
from functions.functions_useful import round_to_base

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})

neurite_types = ['neurites', 'dendrites', 'axons']


# %% set plotting dicts

scatter_plot_dict = {'s' : 0.5}

path_color_dict = {'dendrite' : {'path_color' : colors_dict['primecolor'], 'end_color' : 'r'},
                    'axon' :     {'path_color' : 'b', 'end_color' : 'lightcoral'},
                    'soma' :     {'path_color' : colors_dict['color2'], 'end_color' : colors_dict['color2']}}

# set dict for titles of subplots
alpha_labels_dict = {cell_IDs_dict['MeA'][0]  : {'coordinates': '$\mathregular{A_{i}}$',
                                                 'neurites'   : '$\mathregular{A_{ii}}$',
                                                 'dendrites'  : '$\mathregular{A_{iii}}$',
                                                 'axons'      : '$\mathregular{A_{iv}}$'},
                     cell_IDs_dict['MeA'][1]  : {'coordinates': '$\mathregular{B_{i}}$',
                                                 'neurites'   : '$\mathregular{B_{ii}}$',
                                                 'dendrites'  : '$\mathregular{B_{iii}}$',
                                                 'axons'      : '$\mathregular{B_{iv}}$'},
                     cell_IDs_dict['BAOT'][0] : {'coordinates': '$\mathregular{C_{i}}$',
                                                 'neurites'   : '$\mathregular{C_{ii}}$',
                                                 'dendrites'  : '$\mathregular{C_{iii}}$',
                                                 'axons'      : '$\mathregular{C_{iv}}$'},
                     cell_IDs_dict['BAOT'][1] : {'coordinates': '$\mathregular{D_{i}}$',
                                                 'neurites'   : '$\mathregular{D_{ii}}$',
                                                 'dendrites'  : '$\mathregular{D_{iii}}$',
                                                 'axons'      : '$\mathregular{D_{iv}}$'},}


# %% figure

fig, axs = plt.subplots(nrows = 4,
                        ncols = 4,
                        layout = 'constrained',
                        figsize = get_figure_size(height = 150, width = 150),
                        dpi = 600)

# flatten axes array
axs = axs.flatten()

### cell coordinates ###
for c_idx, cell_ID in enumerate(cell_IDs):
    
    # set axis
    ax = axs[c_idx]
    
    # set coordinates data dict
    cell_coordinates = coordinates_dict[cell_ID]
    
    # get all path_IDs
    cell_path_IDs = cell_coordinates['path_IDs']
    
    # loop through all paths of cell
    for path_ID in cell_path_IDs:
        
        # get path coordinates
        path_all_coordinates = cell_coordinates['all_coor'][cell_coordinates['all_coor']['path_ID'] == path_ID]
        path_end_coordinates = cell_coordinates['end_coor'][cell_coordinates['end_coor']['path_ID'] == path_ID]           
    
        # get label of current path
        cur_path_label = path_all_coordinates['path_label'].iloc[0]
             
        # scatter plots
        # all coordinates
        ax.scatter(x = path_all_coordinates['X'],
                    y = path_all_coordinates['Y'],
                    color = path_color_dict[cur_path_label]['path_color'],
                    **scatter_plot_dict)
    
        # end coordinates
        ax.scatter(x = path_end_coordinates['X'],
                    y = path_end_coordinates['Y'],
                    color = path_color_dict[cur_path_label]['end_color'],
                    **scatter_plot_dict)

    # set aspect ratio
    ax.set_box_aspect(1)
    
    # add label
    ax.text(x = 10, y = 10, s = 'XY', ha = 'left', va = 'top', fontsize = 9)
    
    # edit axes
    # x
    ax.set_xlim([0, cell_coordinates_field_of_view])
    ax.set_xticks(np.arange(0, 590, 250))
    ax.set_xlabel('Width [µm]')
    
    # y
    ax.set_ylim([cell_coordinates_field_of_view, 0])
    ax.set_yticks(ticks = np.arange(0, 590, 250), labels = [])

axs[0].set_ylabel('Height [µm]')
    

### polar plots ###

# change projection for polar plots
for ax in axs[4:16]:
    change_projection(fig, axs, ax, projection = 'polar')
    
# get number of bins
n_bins = len(orientation_labels)

### color code for length of branches ###
# initialise color code
norm_min = 0
norm_max = 1000
cmap_str = 'viridis'
norm = mtl.colors.Normalize(norm_min, norm_max)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# initialise max dict for axes limits
max_dict = {'neurites' : 0, 'dendrites' : 0, 'axons' : 0}

for c_idx, cell_ID in enumerate(cell_IDs):
    
    # get histogram occurances
    cell_occu = pp_occu_neurites_normed_to_neurites[cell_ID]

    # set terminal branches dataframe
    terminal_branches_df = coordinates_dict[cell_ID]['terminal_branches']
    
    # get total number of neurites in bin
    total_n_neurites = terminal_branches_df.drop(index = 1).shape[0]
    total_n_dendrites = terminal_branches_df[terminal_branches_df['path_label'] == 'dendrite'].shape[0]
    total_n_axons = terminal_branches_df[terminal_branches_df['path_label'] == 'axon'].shape[0]
    
    # plot neurites
    neurites_occu = plot_colorcoded_polar_normed(terminal_branches_df,
                                                 max_n_neurites = total_n_neurites, 
                                                 ax = axs[c_idx+4],
                                                 n_bins = n_bins,
                                                 cmap = cmap)    

    # plot dendrites
    dendrites_occu = plot_colorcoded_polar_normed(terminal_branches_df[terminal_branches_df['path_label'] == 'dendrite'],
                                                  max_n_neurites = total_n_dendrites, 
                                                  ax = axs[c_idx+4+4],
                                                  n_bins = n_bins,
                                                  cmap = cmap)    
    
    # plot axons
    axons_occu = plot_colorcoded_polar_normed(terminal_branches_df[terminal_branches_df['path_label'] == 'axon'],
                                              max_n_neurites = total_n_axons, 
                                              ax = axs[c_idx+4+8],
                                              n_bins = n_bins,
                                              cmap = cmap)
    
    # get max of each type and update max dict
    def update_max(max_dict, neurite_type, occu):
        occu_max = occu.max()
        
        if occu_max > max_dict[neurite_type]:
            max_dict[neurite_type] = occu_max
            
    update_max(max_dict, 'neurites', neurites_occu)
    update_max(max_dict, 'dendrites', dendrites_occu)
    update_max(max_dict, 'axons', axons_occu)

   

# add axis titels
for r_idx, row in enumerate(['coordinates', 'neurites', 'dendrites', 'axons']):
    
    for c_idx, cell_ID in enumerate(cell_IDs):
    
        # set axis
        ax = axs[c_idx + r_idx * 4]
        
        if row == 'coordinates':
            titel = alpha_labels_dict[cell_ID][row] + ': ' + cell_ID
        else:
            titel = alpha_labels_dict[cell_ID][row]
        
        ax.set_title(titel, 
                     fontsize = 12,
                     loc='left',
                     x = -0.42)
        
# add row titles
for ax, row in zip(axs[4:16:4], ['Neurites', 'Dendrites', 'Axons']):
    ax.annotate(row, xy=(0, 0.5), 
                xytext=(-ax.yaxis.labelpad - 6, 0),
                xycoords = ax.yaxis.label, 
                textcoords='offset points',
                fontsize=9, 
                ha='right', 
                va='center',
                rotation = 90)


# edit axes
for n_idx, neurite_type in enumerate(neurite_types):
    
    for ax in axs[4*n_idx+4:4*n_idx+8]:
        
        # y
        ymin = 0
        ymax = round_to_base(max_dict[neurite_type], 0.1)
        ymax_label = "{:.0f}".format((ymax * 100)) + ' %'
        
        # apply
        ax.set_ylim([ymin, ymax])
        ax.set_yticks(ticks = np.arange(ymin, ymax+0.01, ymax/2), labels = ['', '', ymax_label])
        ax.tick_params(axis='x', which='major', pad=-3)

        # x axis
        ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
        ax.set_xticklabels(orientation_labels)
    
    
        # grid
        ax.grid(True, alpha = 0.5, color = 'gray')
        ax.set_axisbelow(True)
    

# colorbar
fig.colorbar(cmap, 
             ax = axs[12:16], 
             label = 'Terminal branch length [µm]', 
             orientation='horizontal',
             fraction=0.09,
             aspect=90)


# align labels
fig.align_labels()

# show figure
plt.show()

# save figure
save_figures(fig, figure_name = 'polar_plots-example_cells-figure', 
             save_dir = join(cell_morph_plots_dir, 'polar_plots'),
             darkmode_bool = darkmode_bool, 
             figure_format= 'both')