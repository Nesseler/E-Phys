# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:00:33 2024

@author: nesseler
"""

from cellmorphology.cellmorph_packages import mtl, plt, sbn, pd, np, join

from scipy.stats import circmean

from parameters.directories_win import cell_morph_descrip_dir, table_file, cell_morph_traces_coordinates_dir, cell_morph_plots_dir

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# get cell_IDs:
cell_IDs = MetaData[MetaData['reconstructed'] == 1].index.to_list()

# set dict for all terminal branches measurements files
terminal_branches_dict = {}

# load all terminal branches measurements files
for cell_ID in cell_IDs:
    # set path to terminal branches file
    terminal_branches_measurements_path = join(cell_morph_descrip_dir, 'terminal_branches_measurements', f'{cell_ID}-terminal_branches.xlsx')
    
    # load file
    cell_terminal_branches = pd.read_excel(terminal_branches_measurements_path, index_col='path_ID')
    
    # write to dict
    terminal_branches_dict[cell_ID] = cell_terminal_branches



# set list of neurite types
neurite_types = ['neurites', 'dendrites', 'axons']

# initialise dataframe for circmeans
circular_means = pd.DataFrame(columns = ['neurites-circmean', 'dendrites-circmean', 'axons-circmean'], index = cell_IDs)
circular_means.index.name = 'cell_ID'


for cell_ID in cell_IDs:

    # get cell terminal branches
    cell_terminal_branches = terminal_branches_dict[cell_ID]
    
    # drop soma
    cell_terminal_branches = cell_terminal_branches[cell_terminal_branches['path_label'] != 'soma']
    
    for neurite_type in neurite_types:
        
        # filter paths for label
        if neurite_type == 'neurites':
            cell_terminal_branches_pertype = cell_terminal_branches
        else:
            cell_terminal_branches_pertype = cell_terminal_branches[cell_terminal_branches['path_label'] == neurite_type[:-1]]
        
        # get list of angles in rad
        cell_angles_pertype = cell_terminal_branches_pertype['angle_rad'].to_list()
    
        # calc circmean
        cell_circmean_pertype = circmean(cell_angles_pertype)
        
        # write to dataframe
        circular_means.at[cell_ID, f'{neurite_type}-circmean'] = cell_circmean_pertype
        
# save circ means dataframe to excel
circular_means.to_excel(join(cell_morph_descrip_dir, 'circular_means.xlsx'), index_label = 'cell_ID')
        
        
# %% load cell coordinates for plot

from cellmorphology.functions_cellmorph import clean_OnPath_column_to_path_ID_n_label
    
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


# %% initialize plotting

from functions.functions_plotting import save_figures, get_colors, get_figure_size

from cellmorphology.cellmorph_parameters import cell_coordinates_field_of_view


# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})

# set dict for color of paths in scatter plot of coordinates
path_color_dict = {'dendrite' : {'path_color' : colors_dict['primecolor'], 'end_color' : 'r'},
                    'axon' :     {'path_color' : 'b', 'end_color' : 'lightcoral'},
                    'soma' :     {'path_color' : colors_dict['color2'], 'end_color' : colors_dict['color2']}}

# set dict for scatter plot of coordinates
scatter_plot_dict = {'s' : 0.5}

# set colors for neurite types
neurite_type_color_dict = {'neurites' : 'gray',
                           'dendrites': 'k',
                           'axons'    : 'b'}


for cell_ID in cell_IDs:
    
    # initialise figure
    fig_circmean, ax_circmean = plt.subplots(nrows = 2,
                                             ncols = 1,
                                             layout = 'constrained',
                                             figsize = get_figure_size(width = 75, height = 150),
                                             dpi = 600)
    
    # set title
    fig_circmean.suptitle(f'{cell_ID} circular means')
    
    
    ### cell coordinates ###
    
    # set axis
    ax = ax_circmean[0]
    
    # set axis title
    ax_circmean[0].set_title('A: Cell coordinates', 
                             fontsize = 9,
                             loc = 'left',
                             x = -0.15)
    
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
    
    ax.set_ylabel('Height [µm]')
    
    
    ### circular mean marker ###
    
    # set axis title
    ax_circmean[1].set_title('B: Circular mean on unit circle', 
                             fontsize = 9,
                             loc = 'left',
                             x= -0.15)
    
    
    # plot unit circle
    ax_circmean[1].plot(np.cos(np.linspace(0, 2*np.pi, 500)), np.sin(np.linspace(0, 2*np.pi, 500)), 
                     c=colors_dict['primecolor'], 
                     zorder= 0,
                     label = '_nolegend_')
    
    # plot soma marker
    ax_circmean[1].scatter([0], [0], 
                        c = neurite_type_color_dict['neurites'], 
                        label = 'soma',
                        marker = 'o',
                        s = 60,
                        zorder = 5)
    
    # get cell terminal branches
    cell_terminal_branches = terminal_branches_dict[cell_ID]
    
    # set aspect ratio
    ax_circmean[1].set_box_aspect(1)
    
    # plot neurites per type
    for n_idx, neurite_type in enumerate(neurite_types):
        
        # filter paths for label
        if neurite_type == 'neurites':
            cell_terminal_branches_pertype = cell_terminal_branches
        else:
            cell_terminal_branches_pertype = cell_terminal_branches[cell_terminal_branches['path_label'] == neurite_type[:-1]]
        
        # get list of angles in rad
        cell_angles_pertype = cell_terminal_branches_pertype['angle_rad'].to_list()
        
        # plot angles on circle
        ax_circmean[1].scatter(np.cos(cell_angles_pertype), np.sin(cell_angles_pertype), 
                            c = neurite_type_color_dict[neurite_type],
                            zorder = n_idx+1,
                            label = neurite_type,
                            s = 20)
        
        # set circ mean for cell and type
        cell_circmean_pertype = circular_means.at[cell_ID, f'{neurite_type}-circmean']
        
        # ax_circmean.scatter(np.cos(cell_circmean_pertype), np.sin(cell_circmean_pertype), 
        #                     c = neurite_type_color_dict[neurite_type], 
        #                     label = f'{neurite_type} circular mean',
        #                     marker = 'x',
        #                     s = 40)
        
        # plot circular means
        ax_circmean[1].plot([0, np.cos(cell_circmean_pertype)], [0, np.sin(cell_circmean_pertype)], 
                         c = neurite_type_color_dict[neurite_type], 
                         lw = 1,
                         label = f'{neurite_type} circular mean',
                         marker = 'x',
                         markersize = 6)
        
    
        
    # edit axis
    # y
    ymin = -1
    ymax = 1
    ypad = 0.4
    
    ax_circmean[1].set_ylim([ymin - ypad, ymax + ypad])
    ax_circmean[1].set_yticks(ticks = [ymin, ymax], labels = [])
    ax_circmean[1].spines['left'].set_bounds([ymin, ymax])
    
    ax_circmean[1].set_xlim([ymin - ypad, ymax + ypad])
    ax_circmean[1].set_xticks(ticks = [ymin, ymax], labels = [])
    ax_circmean[1].spines['bottom'].set_bounds([ymin, ymax])
    
    # remove spines
    [ax_circmean[1].spines[spine].set_visible(False) for spine in ['top', 'right']]
        
    # set legend
    # legend = plt.figlegend(loc='upper center', 
    #                        bbox_to_anchor=(0.5, 0),
    #                        ncol=2, 
    #                        prop={'size': 6})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [1,3,5,0,2,4,6]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
               loc='upper center', 
              bbox_to_anchor=(0.5, -0.05),
              ncol=2, 
              prop={'size': 6})
        
    # show plot
    plt.show()
    
    # save figure
    fig_dir = join(cell_morph_plots_dir, 'circular_means')
    save_figures(fig_circmean, 
                 f'{cell_ID}-circular_means', 
                 fig_dir,
                 darkmode_bool = darkmode_bool,
                 figure_format='both')


# %%




import numpy as np

from scipy.stats import circmean

import matplotlib.pyplot as plt

from astropy import stats as astro_stats

angles = np.deg2rad(np.arange(90, 180, 10))

print(astro_stats.rayleightest(angles))

circmean = circmean(angles)

np.rad2deg(circmean)


mean = angles.mean()

np.rad2deg(mean)


plt.plot(np.cos(np.linspace(0, 2*np.pi, 500)),

         np.sin(np.linspace(0, 2*np.pi, 500)),

         c='k')

plt.scatter(np.cos(angles), np.sin(angles), c='k')

plt.scatter(np.cos(circmean), np.sin(circmean), c='b',

            label='circmean')

# plt.scatter(np.cos(mean), np.sin(mean), c='r', label='mean')

plt.legend()

plt.axis('equal')

plt.show()