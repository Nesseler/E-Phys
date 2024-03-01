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
from parameters.directories_win import cell_morph_measures_dir, table_file

from functions.functions_plotting import set_font_sizes, get_colors, get_figure_size


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

# %% trail plot

# set_font_sizes()

parameter = 'No. of terminal branches [Single value]'

darkmode_bool = True

colors_dict, regions_c = get_colors(darkmode_bool)


fig, axs = plt.subplots(nrows = 1,
                        ncols = 1,
                        layout = 'constrained')

violin = sbn.violinplot(data = morph_parameters,
                        y = parameter,
                        inner = 'quart',
                        linewidth = 1,
                        ax = axs)

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

for violin in violin.collections:
    violin.set_edgecolor(colors_dict['primecolor'])
    violin.set_facecolor('None')

    
swarm = sbn.swarmplot(data = morph_parameters,
                      y = parameter, 
                      ax = axs,
                      size = 7,
                      color = colors_dict['primecolor'],
                      hue = 'Region',
                      palette = regions_c)

plt.show()


# plot with separate violins


fig_sep_regions, axs_sep_regions = plt.subplots(nrows = 1,
                                                ncols = 1,
                                                layout = 'constrained')

violin = sbn.violinplot(data = morph_parameters,
                        x = 'Region',
                        y = parameter,
                        inner = 'quart',
                        linewidth = 1,
                        ax = axs_sep_regions)

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
                      hue = 'Region',
                      palette = regions_c)












