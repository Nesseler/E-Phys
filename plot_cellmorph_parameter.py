# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:52:58 2024

@author: nesseler
"""

import pandas as pd
import seaborn as sbn
import matplotlib.pyplot as plt
from os.path import join

from getter.get_onlyfiles_list import get_onlyfiles_list
from parameters.directories_win import cell_morph_measures_dir, table_file, cell_morph_figures_dir

from functions.functions_plotting import set_font_sizes, get_colors, get_figure_size, save_figures


# get list of all files in directory
onlyfiles = get_onlyfiles_list(cell_morph_measures_dir)

# get cell_IDs from filenames
cell_IDs = [f_str.replace('.csv', '') for f_str in onlyfiles]

# initialize common dataframe
morph_parameters = pd.DataFrame()

for i, cell_ID in enumerate(cell_IDs):
    cell_parameters = pd.read_csv(join(cell_morph_measures_dir, cell_ID + '.csv'),
                                   encoding='unicode_escape')
    
    cell_parameters.rename(index = {0 : cell_ID}, inplace = True)
    
    morph_parameters = pd.concat([morph_parameters, cell_parameters], axis = 0)
    

morph_parameters.drop(columns = ['\\'], inplace = True)

morph_parameters.index.name = 'cell_ID'


# %% import meta data sheet

MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]


# %% combine morph parameters and MetaData

morph_parameters = pd.concat([morph_parameters, MetaData], axis = 1)

# export to combined list
morph_parameters.T.to_excel(join(cell_morph_measures_dir, 'morph_parameters.xlsx'), index = 'cell_ID')

# %% trail plot

# set_font_sizes()

parameters = ['No. of branches [Single value]',
              'No. of terminal branches [Single value]',
              'No. of primary branches [Single value]',
              'Convex hull: Size (µm³) [Single value]',
              'Convex hull: Centroid-root distance (µm) [Single value]',
              'Convex hull: Roundness [Single value]',
              'Convex hull: Elongation (µm) [Single value]',
              'Complexity index [Single value]',
              'Horton-Strahler bifurcation ratio [Single value]'
              ]


darkmode_bool = True

colors_dict, regions_c = get_colors(darkmode_bool)



for parameter in parameters:

    fig_sep_regions, axs_sep_regions = plt.subplots(nrows = 1,
                                                    ncols = 1,
                                                    layout = 'constrained')
    
    violin = sbn.violinplot(data = morph_parameters,
                            x = 'Region',
                            y = parameter,
                            inner = 'quart',
                            linewidth = 1,
                            hue = 'put. morph category',
                            ax = axs_sep_regions,
                            palette = 'deep')
    
    for l in violin.lines:
        l.set_color(colors_dict['primecolor'])
    
    for violin in violin.collections:
        violin.set_edgecolor(colors_dict['primecolor'])
        violin.set_facecolor('None')
    
        
    swarm = sbn.swarmplot(data = morph_parameters,
                          x = 'Region',
                          y = parameter, 
                          ax = axs_sep_regions,
                          size = 7,
                          color = colors_dict['primecolor'],
                          hue = 'put. morph category',
                          palette = 'deep',
                          dodge = True)
    
    
    plt.grid(False)
    
    # save_figures(fig_sep_regions, f'{parameter.replace(" ", "_").replace(":", "_")}', cell_morph_figures_dir, darkmode_bool)

    plt.show()






