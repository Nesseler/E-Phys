# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 15:22:59 2024

@author: nesseler
"""

from cellmorphology.cellmorph_packages import mtl, plt, pd, np, join

from parameters.directories_win import table_file, cell_morph_traces_coordinates_dir, cell_morph_descrip_dir, cell_morph_plots_dir

from cellmorphology.functions_cellmorph import clean_OnPath_column_to_path_ID_n_label

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# get cell_IDs:
cell_IDs = MetaData[MetaData['reconstructed'] == 1].index.to_list()

# set list of neurite types
neurite_types = ['neurites', 'dendrites', 'axons']


# %% load cell coordinates

# all coordinates dict
coordinates_dict = {}

# create list of columns for dataframe
columns = [neurite_type + '-' + col for neurite_type in neurite_types for col in ['x_min', 'y_min', 'width', 'height']]

# initialise dataframe for height, width measurements
height_width_df = pd.DataFrame(columns=columns, index = cell_IDs)
height_width_df.index.name = 'cell_IDs'


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
    
    # get coordinates min and max
    for neurite_type in neurite_types:
        
        # limit coordinates to specific neurite type
        if neurite_type == 'neurites':
            cell_allcoordinates_pertype = cell_allcoordinates
        else:
            cell_allcoordinates_pertype = cell_allcoordinates[cell_allcoordinates['path_label'] == neurite_type[:-1]]
        
        # get x and ymin
        xmin = cell_allcoordinates_pertype['X'].min()
        ymin = cell_allcoordinates_pertype['Y'].min()
        xmax = cell_allcoordinates_pertype['X'].max()
        ymax = cell_allcoordinates_pertype['Y'].max()
        
        # write to dataframe
        # mins
        height_width_df.at[cell_ID, neurite_type + '-' + 'x_min'] = xmin
        height_width_df.at[cell_ID, neurite_type + '-' + 'y_min'] = ymin
        
        # calc width and height
        height_width_df.at[cell_ID, neurite_type + '-' + 'width'] = xmax - xmin
        height_width_df.at[cell_ID, neurite_type + '-' + 'height'] = ymax - ymin
    




# %% initialize plotting

from functions.functions_plotting import save_figures, get_colors, get_figure_size
from cellmorphology.cellmorph_parameters import cell_coordinates_field_of_view

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})



cell_ID = 'E-087'

vplots = True


# %% set plotting dicts

scatter_plot_dict = {'s' : 0.5}

path_color_dict = {'dendrite' : {'path_color' : colors_dict['primecolor'], 'end_color' : 'r'},
                   'axon' :     {'path_color' : 'b', 'end_color' : 'lightcoral'},
                   'soma' :     {'path_color' : colors_dict['color2'], 'end_color' : colors_dict['color2']}}

neurite_color_dict = {'neurites'  : 'grey',
                      'dendrites' : colors_dict['primecolor'],
                      'axons'     : 'b'}


# %% plot

if vplots:

    for cell_ID in cell_IDs: #['E-137']:   
        
        fig, ax = plt.subplots(nrows = 1,
                               ncols = 1,
                               layout = 'constrained',
                               figsize = get_figure_size(width = 100, height = 100),
                               dpi = 600)
        
        # set figure title
        fig.suptitle(f'{cell_ID}')
        
            
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
            
            
        # max points
        enclosing_points_x = [np.min(cell_coordinates['all_coor']['X']), np.max(cell_coordinates['all_coor']['X'])] 
        enclosing_points_y = [np.min(cell_coordinates['all_coor']['Y']), np.max(cell_coordinates['all_coor']['Y'])] 
        
        # get height, width of cell
        cell_height_width = height_width_df.loc[cell_ID,:]

        # create measurements label
        m_label = ''
        
        # loop through neurite types and plot enclosing rectangles
        for neurite_type in neurite_types:   
            # initialize rectanlge
            rectangle = mtl.patches.Rectangle(xy = (cell_height_width[f'{neurite_type}-x_min'], cell_height_width[f'{neurite_type}-y_min']),
                                              width = cell_height_width[f'{neurite_type}-width'],
                                              height = cell_height_width[f'{neurite_type}-height'],
                                              facecolor = 'None',
                                              edgecolor = neurite_color_dict[neurite_type])
        
            # add rectangle as patch
            ax.add_patch(rectangle)
            
            # update measurements label
            m_label = m_label + f'\n{neurite_type} width [µm] = {"{:.2f}".format(cell_height_width[f"{neurite_type}-width"])}\n{neurite_type} height [µm] = {"{:.2f}".format(cell_height_width[f"{neurite_type}-height"])}\n'
            
            
        # add measurements to plot
        ax.text(x = 580, y = 0, 
                s = m_label, 
                ha = 'right', 
                va = 'top',
                size = 6)
        
        # set aspect ratio
        ax.set_box_aspect(1)
        
        # add label
        ax.text(x = 10, y = 10, s = 'XY', ha = 'left', va = 'top', fontsize = 9)
        
        # edit axes
        # x
        ax.set_xlim([0, cell_coordinates_field_of_view])
        ax.set_xticks(np.arange(0, 590, 250))
        ax.set_xticks(np.arange(0, 590, 50), minor = True)
        ax.set_xlabel('Width [µm]')
        
        # y
        ax.set_ylim([cell_coordinates_field_of_view, 0])
        ax.set_yticks(ticks = np.arange(0, 590, 250))
        ax.set_yticks(np.arange(0, 590, 50), minor = True)
        ax.set_ylabel('Height [µm]')
        
        # show plot
        plt.show()
        
        # save figure
        fig_dir = join(cell_morph_plots_dir, 'height_width', 'height_width_per_cell')
        save_figures(fig, 
                     f'{cell_ID}-width_height', 
                     fig_dir,
                     darkmode_bool = darkmode_bool,
                     figure_format='png')


# %% save dataframe

height_width_df.to_excel(join(cell_morph_descrip_dir, 'height_width.xlsx'), index_label='cell_ID')