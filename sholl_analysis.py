# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 11:34:34 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import seaborn as sbn
import pandas as pd
from os.path import join
import numpy as np

from parameters.directories_win import cell_morph_traces_sholl_dir, table_file, cell_morph_plots_dir, cell_morph_descrip_dir

from functions.functions_import import get_onlyfiles_list
from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size


# %% parameters

# define sholl radius step size
sholl_step_size = 1 # um

# directory of table files
sholl_neurites_tables_dir = join(cell_morph_traces_sholl_dir, 'sholl_tables_neurites')
sholl_dendrites_tables_dir = join(cell_morph_traces_sholl_dir, 'sholl_tables_dendrites')
# sholl_axons_tables_dir = join(cell_morph_traces_sholl_dir, 'sholl_tables_axons')


# %% 

# get files in directory
sholl_neurites_onlyfiles = get_onlyfiles_list(sholl_neurites_tables_dir)
sholl_dendrites_onlyfiles = get_onlyfiles_list(sholl_dendrites_tables_dir)
# sholl_axons_onlyfiles = get_onlyfiles_list(sholl_axons_tables_dir)

# remove summary table from filenames
sholl_neurites_onlyfiles = [fname for fname in sholl_neurites_onlyfiles if '_profile' in fname]
sholl_dendrites_onlyfiles = [fname for fname in sholl_dendrites_onlyfiles if '_profile' in fname]
# sholl_axons_onlyfiles = [fname for fname in sholl_axons_onlyfiles if '_profile' in fname]

# get cell_IDs
cell_IDs = ['E' + fname[16:20] for fname in sholl_neurites_onlyfiles]
# cell_IDs_axons = ['E' + fname[16:20] for fname in sholl_axons_onlyfiles]

# create dataframe with neurites filenames
sholl_filenames_neurites = pd.DataFrame({'neurites' : sholl_neurites_onlyfiles},
                               index = cell_IDs)

# create dataframe with axons filenames
sholl_filenames_dendrites = pd.DataFrame({'dendrites' : sholl_dendrites_onlyfiles},
                                         index = cell_IDs)

# concatenate both filenames dataframe
sholl_filenames = pd.concat([sholl_filenames_neurites, sholl_filenames_dendrites], axis = 1)

# quality control (loading filenames)
for cell_ID, row in zip(sholl_filenames.index, sholl_filenames.values):
    # check index and neurites
    if not cell_ID[1:] == row[0][16:20]:
        raise ValueError(f'Filenames and cell_IDs do not match! {cell_ID}')
    # check index and axons
    if not type(row[1]) == float:
        if not cell_ID[1:] == row[1][16:20]:
            raise ValueError(f'Filenames and cell_IDs do not match! {cell_ID}')

# prepare dataframes for all sholl profiles
## neurites
sholl_neurites = pd.DataFrame(columns = cell_IDs,
                              index = np.arange(0, 500 + sholl_step_size, sholl_step_size))
sholl_neurites.index.name = 'radius'

## axons
sholl_dendrites = pd.DataFrame(columns = sholl_filenames_dendrites.index,
                           index = np.arange(0, 500 + sholl_step_size, sholl_step_size))
sholl_dendrites.index.name = 'radius'


# cell_ID = 'E-137'

for cell_ID in cell_IDs:
    ## neurites
    # load sholl profiles
    cell_sholl_neurites = pd.read_csv(join(sholl_neurites_tables_dir, sholl_filenames.at[cell_ID, 'neurites']))
    cell_sholl_neurites.set_index('Radius', inplace = True)
    
    # write to dataframe
    sholl_neurites[cell_ID] = cell_sholl_neurites['Inters.']
    
    ## dendrites
    # load sholl profiles
    cell_sholl_dendrites = pd.read_csv(join(sholl_dendrites_tables_dir, sholl_filenames.at[cell_ID, 'dendrites']))
    cell_sholl_dendrites.set_index('Radius', inplace = True)
    
    # write to dataframe
    sholl_dendrites[cell_ID] = cell_sholl_dendrites['Inters.']
    
    # ## axons
    # # load sholl profiles if available
    # if not type(sholl_filenames.at[cell_ID, 'axons']) == float:
    #     cell_sholl_axons = pd.read_csv(join(sholl_axons_tables_dir, sholl_filenames.at[cell_ID, 'axons']))
    #     cell_sholl_axons.set_index('Radius', inplace = True)
    
    #     # write to dataframe
    #     sholl_axons[cell_ID] = cell_sholl_axons['Inters.']


# %% remove first row (SNT: sholl analysis error: sometimes by radius of 0)

sholl_neurites = sholl_neurites.iloc[1:]
sholl_dendrites = sholl_dendrites.iloc[1:]


# %% calculate dendrite sholl profile as difference between neurites and axons

sholl_axons = sholl_neurites.sub(sholl_dendrites, fill_value = 0)

# drop columns (cells) with no axons, i.e. all zeros in column
not_null_cell_IDs = sholl_axons.sum() != 0
sholl_axons = sholl_axons.loc[:, not_null_cell_IDs]

# set cell_IDs with axons from sholl tabel
cell_IDs_axons = sholl_axons.columns.to_list()


# %% get sholl metrics

## neurites
sholl_metrics_neurites = pd.DataFrame({'critical_radius' : sholl_neurites.idxmax(),
                                       'max_intersections' : sholl_neurites.max(),
                                       'enclosing_radius' : [sholl_neurites[cell_ID].dropna().index[-1] for cell_ID in cell_IDs]})

## dendrites
sholl_metrics_dendrites = pd.DataFrame({'critical_radius' : sholl_dendrites.idxmax(),
                                        'max_intersections' : sholl_dendrites.max(),
                                        'enclosing_radius' : [sholl_dendrites[cell_ID].dropna().index[-1] for cell_ID in cell_IDs]})

## neurites
sholl_metrics_axons = pd.DataFrame({'critical_radius' : sholl_axons.idxmax(),
                                    'max_intersections' : sholl_axons.max(),
                                    'enclosing_radius' : [sholl_axons[cell_ID].dropna().index[-1] for cell_ID in cell_IDs_axons]})

# fill nan values with zeros
sholl_neurites.fillna(0, inplace = True)
sholl_dendrites.fillna(0, inplace = True)
sholl_axons.fillna(0, inplace = True)



# %% save sholl profiles & metrics dataframes

if True:
    # save profiles metrics
    sholl_neurites.to_excel(join(cell_morph_descrip_dir, 'sholl_profiles_neurites.xlsx'), index_label = 'Radius')
    sholl_dendrites.to_excel(join(cell_morph_descrip_dir, 'sholl_profiles_dendrites.xlsx'), index_label = 'Radius')
    sholl_axons.to_excel(join(cell_morph_descrip_dir, 'sholl_profiles_axons.xlsx'), index_label = 'Radius')

    # save sholl metrics
    sholl_metrics_neurites.to_excel(join(cell_morph_descrip_dir, 'sholl_metrics_neurites.xlsx'), index_label = 'cell_ID')
    sholl_metrics_dendrites.to_excel(join(cell_morph_descrip_dir, 'sholl_metrics_dendrites.xlsx'), index_label = 'cell_ID')
    sholl_metrics_axons.to_excel(join(cell_morph_descrip_dir, 'sholl_metrics_axons.xlsx'), index_label = 'cell_ID')


# %% single cell sholl profiles

from functions.functions_useful import round_up_to_base



# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)
from cellmorphology.cellmorph_colors import neurite_color_dict

# set font size 
mtl.rcParams.update({'font.size': 9})

# set export directory for cell sholl plots
cell_sholl_plots_dir = join(cell_morph_plots_dir, 'sholl_plots', 'sholl_profile')

# cell_ID = 'E-137'

for cell_ID in cell_IDs:

    fig_cell, axs_cell = plt.subplots(nrows = 1,
                                      ncols = 1,
                                      figsize = get_figure_size(100, 50),
                                      dpi = 600)
    
    fig_cell.suptitle(f'{cell_ID} sholl profile')
    
    # plot sholl profile
    axs_cell.plot(sholl_neurites[cell_ID], neurite_color_dict['all']['neurites'], label = 'neurites')
    
    axs_cell.plot(sholl_dendrites[cell_ID], neurite_color_dict['all']['dendrites'], label = 'dendrites')
    
    # plot axon profile if cell has axon
    if cell_ID in cell_IDs_axons:
        axs_cell.plot(sholl_axons[cell_ID], neurite_color_dict['all']['axons'], label = 'axons')
    
    # plot legend
    plt.legend()
    
    
    # x axis
    axs_cell.set_xlabel('Radius [µm]')
    xmin = 0
    xmax = round_up_to_base(sholl_metrics_dendrites.at[cell_ID, 'enclosing_radius'], 100)
    xpad = (xmax - xmin) * 0.025
    
    axs_cell.set_xlim(xmin - xpad, xmax)
    axs_cell.set_xticks(np.arange(0, xmax + 1, 50))
    axs_cell.spines['bottom'].set_bounds([0, xmax])
    
    # y axis
    axs_cell.set_ylabel('Number of intersections [#]')
    ymin = 0
    ymax = round_up_to_base(sholl_metrics_dendrites.at[cell_ID, 'max_intersections'], 10)
    ypad = (ymax - ymin) * 0.025
    
    axs_cell.set_ylim(0 - ypad, ymax + ypad)
    axs_cell.set_yticks(np.arange(0, ymax + 1, 5))
    axs_cell.spines['left'].set_bounds([0, ymax])
    
    # remove spines
    [axs_cell.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # save figure
    save_figures(fig_cell, f'{cell_ID}-sholl_profile', cell_sholl_plots_dir, figure_format= 'both')
    
    # show figure
    plt.show()












