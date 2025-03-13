# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:34:01 2024

@author: nesseler
"""

import matplotlib as mtl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sbn

from functions.functions_plotting import get_figure_size, save_figures


spines_table_file = '//fileserver2/AG Spehr BigData/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA/Plots/spines_plots' + '/' + 'spines.xlsx'

# load dataframe
spines_absolute_df = pd.read_excel(spines_table_file, 
                                   sheet_name = 'barplot',
                                   index_col = 'Region')

spines_relative_df = pd.read_excel(spines_table_file, 
                                   sheet_name = 'barplot_percentage',
                                   index_col = 'Region')

# get list of regions
regions = spines_absolute_df.index.to_list()

# spinyness categories list
spinyness = list(reversed(spines_absolute_df.columns.to_list()))

# set a list for x coordinates for regions
x = [0, 1, 2]

# initialise 
spines_color_dict = {'both' : {'high' : '#776B5D' , 'moderate' : '#B0A695', 'low': '#EBE3D5'},
                     'MeA' :  {'high' : '#A02334' , 'moderate' : '#EB5B00', 'low': '#FFB200'},
                     'BAOT' : {'high' : '#201E43' , 'moderate' : '#134B70', 'low': '#508C9B'}}

# set font size
mtl.rcParams.update({'font.size': 9})

# initialise figure

fig, axs = plt.subplots(nrows = 1,
                        ncols = 2, 
                        layout = 'tight',
                        dpi = 600,
                        figsize = get_figure_size(width = 150, height = 90))



for r_idx, region in enumerate(regions):
    
    # define bottom values before loop
    bottom_abs = [0]
    bottom_rel = [0]
    
    # set region label
    if region == 'both':
        region_label = 'Both regions'
    else:
        region_label = region
    
    for spinyness_cat in spinyness:

        # plot relative number
        axs[0].bar(x = x[r_idx],
                   height = spines_relative_df.at[region, spinyness_cat],
                   bottom = bottom_rel,
                   color = spines_color_dict[region][spinyness_cat],
                   label = region_label + ': ' + spinyness_cat)
        
        # plot absolute numbers
        axs[1].bar(x = x[r_idx],
                   height = spines_absolute_df.at[region, spinyness_cat],
                   bottom = bottom_abs,
                   color = spines_color_dict[region][spinyness_cat])
        
        # update and add height to bottom list
        bottom_abs = bottom_abs + spines_absolute_df.at[region, spinyness_cat]
        bottom_rel = bottom_rel + spines_relative_df.at[region, spinyness_cat]



legend =plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 0),
                      ncol=3)



# edit axis

# loop through panes
for ax in axs:
    
    # x
    ax.set_xticks(ticks = x,
                  labels = ['Both\nregions', 'MeA', 'BAOT'])
    
    ax.set_xlabel('Regions')
    ax.spines['bottom'].set_bounds([0, 2])
    
    
# y
ymin = 0
ymax = 100
ystep = 20
ystepminor = 5
ypad = 2

axs[0].set_ylabel('Number of cells [%]')
axs[0].set_ylim([ymin - ypad, ymax])
axs[0].set_yticks(ticks = np.arange(ymin, ymax+1, ystep))
axs[0].set_yticks(ticks = np.arange(ymin, ymax+1, ystepminor), minor = True)
axs[0].spines['left'].set_bounds([ymin, ymax])


ymax = 80

axs[1].set_ylabel('Number of cells [#]')
axs[1].set_ylim([ymin - ypad, ymax])
axs[1].set_yticks(ticks = np.arange(ymin, ymax+1, ystep))
axs[1].set_yticks(ticks = np.arange(ymin, ymax+1, ystepminor), minor = True)
axs[1].spines['left'].set_bounds([ymin, ymax])

# remove spines
for ax in axs:
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)


# set figure directory
fig_dir = "//fileserver2/AG Spehr BigData/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA/Plots/spines_plots/"

# save figure
fig.savefig(fig_dir + 'spines_barplots_combined.png', bbox_extra_artists=[legend], bbox_inches='tight', format = 'png')
fig.savefig(fig_dir + 'spines_barplots_combined.svg', bbox_extra_artists=[legend], bbox_inches='tight', format = 'svg')