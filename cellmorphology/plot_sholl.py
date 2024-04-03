# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 17:36:31 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import seaborn as sbn
import pandas as pd
from os.path import join
import numpy as np

from parameters.directories_win import cell_morph_traces_sholl_dir, table_file, cell_morph_plots_dir

from getter.get_onlyfiles_list import get_onlyfiles_list

from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size

# define sholl radius step size
sholl_step_size = 1 # um

# define directory for figures
sholl_plots_dir = join(cell_morph_plots_dir, 'sholl_plots')

# get onlyfiles list
onlyfiles = get_onlyfiles_list(cell_morph_traces_sholl_dir)

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# extract only filenames with all coordinates
onlyfiles_sholl_profile = [f for f in onlyfiles if '_profile' in f]

# create combined sholl profile dataframe
all_sholl_profiles = pd.DataFrame(index = np.arange(0, 590 + 1, sholl_step_size))

# loop through all files
for filename in onlyfiles_sholl_profile:
    # get cell_ID
    cell_ID = 'E' + filename[16:20]
    
    # read sholl profile from csv
    sholl_csv = pd.read_csv(join(cell_morph_traces_sholl_dir, filename))
    
    # get only intersection counts
    intersections = sholl_csv['Inters.']
    
    # rename column to cell_ID
    intersections.name = cell_ID
    
    # concatenate to all sholl profiles
    all_sholl_profiles = pd.concat([all_sholl_profiles, intersections], axis = 1)

# get all cell_IDs
cell_IDs = all_sholl_profiles.columns.to_list()
    
# rename index
all_sholl_profiles.index.name = 'Radius'

# dataframe with additional metrics
all_sholl_metrics = pd.DataFrame(index = all_sholl_profiles.columns.to_list())

all_sholl_metrics['max_intersections'] = all_sholl_profiles.max(axis = 0)
all_sholl_metrics['idxmax_intersections'] = all_sholl_profiles.idxmax(axis = 0)

all_sholl_metrics['end_radius'] = [len(all_sholl_profiles[c][all_sholl_profiles[c].notnull()])-1 for c in all_sholl_profiles.columns]

# concatenate metadata with all_sholl_metrics
all_sholl_metrics = pd.concat([all_sholl_metrics, MetaData.loc[cell_IDs, 'Region']], axis = 1)

# calculate mean, median and std
all_sholl_profiles_metrics = pd.DataFrame(index = all_sholl_profiles.index)
all_sholl_profiles_metrics.index.name = 'Radius'
all_sholl_profiles_metrics['mean_intersections'] = all_sholl_profiles.mean(axis = 1)
all_sholl_profiles_metrics['median_intersections'] = all_sholl_profiles.median(axis = 1)
all_sholl_profiles_metrics['std_intersections'] = all_sholl_profiles.std(axis = 1)

# %% plot

darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)
set_font_sizes()

fig_sholl, axs_sholl = plt.subplots(nrows = 3, 
                                    ncols = 1,
                                    height_ratios = [5, 1, 1],
                                    sharex = True,
                                    layout = 'constrained')

# sholl profiles
lines = axs_sholl[0].plot(all_sholl_profiles, c = 'gray', zorder = 0, label = 'Sholl profiles')

std_shade = axs_sholl[0].fill_between(x = all_sholl_profiles_metrics.index.to_list(), 
                          y1 = all_sholl_profiles_metrics['mean_intersections'] - all_sholl_profiles_metrics['std_intersections'],
                          y2 = all_sholl_profiles_metrics['mean_intersections'] + all_sholl_profiles_metrics['std_intersections'], 
                          color = 'blue',
                          alpha = 0.25,
                          edgecolor = None,
                          zorder = 1,
                          label = 'std')



mean_line = axs_sholl[0].plot(all_sholl_profiles_metrics.index.to_list(), all_sholl_profiles_metrics['mean_intersections'],
                              color = colors_dict['primecolor'],
                              linewidth = 1,
                              label = 'mean')

axs_sholl[0].legend(handles = [lines[0], std_shade, mean_line[0]])

# y axis
axs_sholl[0].set_ylim([0, 40])
axs_sholl[0].set_yticks(np.arange(0, 40 + 1, 10))
axs_sholl[0].set_yticks(np.arange(0, 40 + 1, 5), minor = True)

# radius of max
axs_sholl[1].text(x = 450, y = -0.5, s = 'Critical radius', ha = 'right', va = 'top', fontsize = 8)

violin = sbn.violinplot(data = all_sholl_metrics,
                        x = 'idxmax_intersections',
                        ax = axs_sholl[1],
                        inner = 'quart',
                        linewidth = 1)

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

for violin in violin.collections:
    violin.set_edgecolor(colors_dict['primecolor'])
    violin.set_facecolor('None')

swarm = sbn.swarmplot(data = all_sholl_metrics,
                      x = 'idxmax_intersections',
                      ax = axs_sholl[1], 
                      color = 'gray',
                      size = 4)

axs_sholl[1].set_xlabel('')


# end radius
axs_sholl[2].text(x = 450, y = -0.5, s = 'Enclosing radius', ha = 'right', va = 'top', fontsize = 8)

violin = sbn.violinplot(data = all_sholl_metrics,
                        x = 'end_radius',
                        ax = axs_sholl[2],
                        inner = 'quart',
                        linewidth = 1)

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

for violin in violin.collections:
    violin.set_edgecolor(colors_dict['primecolor'])
    violin.set_facecolor('None')

swarm = sbn.swarmplot(x = all_sholl_metrics['end_radius'], 
                      ax = axs_sholl[2],
                      color = 'gray',
                      size = 4)


# x axis
xmin = 0
xmax = 450
axs_sholl[2].set_xlim([xmin-5, xmax+5])
axs_sholl[2].set_xticks(np.arange(xmin, xmax + 1, 50))
axs_sholl[2].spines['bottom'].set_bounds([xmin, xmax])
axs_sholl[2].set_xlabel('Radius [µm]')

# despine bottom panels
[ax.spines['left'].set_visible(False) for ax in axs_sholl[1::]]

for ax in axs_sholl:
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    ax.spines['bottom'].set_bounds([xmin, xmax])

# save figure
save_figures(fig_sholl, 'Sholl_profiles_figure-all', sholl_plots_dir, darkmode_bool)

# %% plot with color coded regions

fig_cc, axs_cc = plt.subplots(nrows = 4, 
                              ncols = 1,
                              height_ratios = [3, 3, 1, 1],
                              sharex = True,
                              layout = 'constrained',
                              figsize = get_figure_size(width = 165.5))

set_font_sizes()

# sholl profiles

regions = ['MeA', 'BAOT']

for i_region, region in enumerate(regions):
    
    # get dataframe for specified region in loo
    region_sholl_metrics = all_sholl_metrics[all_sholl_metrics['Region'] == region]
    
    # get cell IDs in region
    region_cell_IDs = region_sholl_metrics.index.to_list()
    
    # get sholl profiles with cell IDs
    region_sholl_profils = all_sholl_profiles[region_cell_IDs]
    
    # recalculate profil metrics
    region_sholl_profil_metrics = pd.DataFrame(index = region_sholl_profils.index)
    region_sholl_profil_metrics.index.name = 'Radius'
    region_sholl_profil_metrics['mean_intersections'] = region_sholl_profils.mean(axis = 1)
    region_sholl_profil_metrics['median_intersections'] = region_sholl_profils.median(axis = 1)
    region_sholl_profil_metrics['std_intersections'] = region_sholl_profils.std(axis = 1)
    
    # define axis
    ax = axs_cc[i_region]
    
    ax.set_title(region)

    lines = ax.plot(region_sholl_profils, 
                              c = region_colors[region], 
                              zorder = 0, 
                              label = 'Sholl profiles')

    std_shade = ax.fill_between(x = region_sholl_profil_metrics.index.to_list(), 
                                          y1 = region_sholl_profil_metrics['mean_intersections'] - region_sholl_profil_metrics['std_intersections'],
                                          y2 = region_sholl_profil_metrics['mean_intersections'] + region_sholl_profil_metrics['std_intersections'], 
                                          color = 'gray',
                                          alpha = 0.5,
                                          edgecolor = None,
                                          zorder = 1,
                                          label = 'std')



    mean_line = ax.plot(region_sholl_profil_metrics.index.to_list(), region_sholl_profil_metrics['mean_intersections'],
                                  color = colors_dict['primecolor'],
                                  linewidth = 1,
                                  label = 'mean')

    ax.legend(handles = [lines[0], std_shade, mean_line[0]], prop={'size': 8})

    # y axis
    ax.set_ylim([0, 40])
    ax.set_yticks(np.arange(0, 40 + 1, 10))
    ax.set_yticks(np.arange(0, 40 + 1, 5), minor = True)


# radius of max
axs_cc[2].text(x = 450, y = -0.5, s = 'Critical radius', ha = 'right', va = 'top', fontsize = 8)

violin = sbn.violinplot(data = all_sholl_metrics,
                        x = 'idxmax_intersections',
                        y = 'Region',
                        ax = axs_cc[2],
                        inner = 'quart',
                        linewidth = 1)

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

for violin in violin.collections:
    violin.set_edgecolor(colors_dict['primecolor'])
    violin.set_facecolor('None')

swarm = sbn.swarmplot(data = all_sholl_metrics,
                      x = 'idxmax_intersections',
                      y = 'Region',
                      hue = 'Region',
                      ax = axs_cc[2], 
                      color = 'gray',
                      size = 4,
                      palette = region_colors)

axs_cc[2].set_xlabel('')
axs_cc[2].set_ylabel('')
axs_cc[2].legend().set_visible(False)


# end radius
axs_cc[3].text(x = 450, y = -0.5, s = 'Enclosing radius', ha = 'right', va = 'top', fontsize = 8)

violin = sbn.violinplot(data = all_sholl_metrics,
                        x = 'end_radius',
                        y = 'Region',
                        ax = axs_cc[3],
                        inner = 'quart',
                        linewidth = 1)

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

for violin in violin.collections:
    violin.set_edgecolor(colors_dict['primecolor'])
    violin.set_facecolor('None')

swarm = sbn.swarmplot(data = all_sholl_metrics,
                      x = 'end_radius',
                      y = 'Region',
                      hue = 'Region',
                      ax = axs_cc[3],
                      color = 'gray',
                      size = 4,
                      palette = region_colors)

axs_cc[3].legend().set_visible(False)
axs_cc[3].set_ylabel('')

# x axis
xmin = 0
xmax = 450
axs_cc[3].set_xlim([xmin-5, xmax+5])
axs_cc[3].set_xticks(np.arange(xmin, xmax + 1, 50))
axs_cc[3].spines['bottom'].set_bounds([xmin, xmax])
axs_cc[3].set_xlabel('Radius [µm]')

# despine bottom panels
[axs_cc[2].spines[spine].set_visible(False) for spine in ['left', 'bottom']]
axs_cc[2].tick_params(axis = 'x', size = 0)
axs_cc[3].spines['left'].set_visible(False)
 
for ax in axs_cc:
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    ax.spines['bottom'].set_bounds([xmin, xmax])

# save figure
save_figures(fig_cc, 'Sholl_profiles_figure-cc_regions', sholl_plots_dir, darkmode_bool)


# %% enclosing vs critical radius

# both regions combined    
fig_cr_er_c_com, ax_cv_er_com = plt.subplots(nrows = 1, 
                                   ncols = 1,
                                   layout = 'constrained',
                                   figsize = get_figure_size(width = 165.5))



sbn.scatterplot(data = all_sholl_metrics,
                x = 'end_radius',
                y = 'idxmax_intersections',
                hue = 'Region', 
                palette = region_colors,
                ax = ax_cv_er_com)


# plt.gca().legend(title = 'Region', prop={'size': 8})


for cell_ID in all_sholl_metrics.index:

    ax_cv_er_com.text(x = all_sholl_metrics.at[cell_ID, 'end_radius'],
                      y = all_sholl_metrics.at[cell_ID, 'idxmax_intersections']-2,
                      s = cell_ID,
                      ha = 'center',
                      va = 'top',
                      fontsize = 8,
                      color = region_colors[all_sholl_metrics.at[cell_ID, 'Region']])
    
# x axis
xmin = 0
xmax = 450
ax_cv_er_com.set_xlim([xmin-5, xmax+5])
ax_cv_er_com.set_xticks(np.arange(xmin, xmax + 1, 100))
ax_cv_er_com.set_xticks(np.arange(xmin, xmax + 1, 25), minor = True)
ax_cv_er_com.spines['bottom'].set_bounds([xmin, xmax])
ax_cv_er_com.set_xlabel('Enclosing radius [µm]')

# y axis
ymin = 0
ymax = 200
ax_cv_er_com.set_ylim([ymin-5, ymax+5])
ax_cv_er_com.set_yticks(np.arange(ymin, ymax + 1, 50))
ax_cv_er_com.set_yticks(np.arange(ymin, ymax + 1, 25), minor = True)
ax_cv_er_com.spines['left'].set_bounds([ymin, ymax])
ax_cv_er_com.set_ylabel('Critical radius [µm]')

# despine
[ax_cv_er_com.spines[spine].set_visible(False) for spine in ['top', 'right']]

# save figure
save_figures(fig_cr_er_c_com, 'Enclosing_v_critical_radius-regions', sholl_plots_dir, darkmode_bool)


# %% enclosing vs critical radius for each region (color coded max number of intersections)

for i_region, region in enumerate(regions):

    region_sholl_metrics = all_sholl_metrics[all_sholl_metrics['Region'] == region]
    
    fig_cr_er_cc, ax_cv_er_cc = plt.subplots(nrows = 1, 
                                               ncols = 1,
                                               layout = 'constrained',
                                               figsize = get_figure_size(width = 185.5))
    
    fig_cr_er_cc.suptitle(region, color = region_colors[region])
    
    ### color code for length of branches ###
    # initialise color code
    norm_min = 0
    norm_max = 40
    cmap_str = 'hot'
    norm = mtl.colors.Normalize(norm_min, norm_max)
    cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
    
    # colorbar
    fig_cr_er_cc.colorbar(cmap, ax = ax_cv_er_cc, label = 'Max. number of intersections')
    
    ax_cv_er_cc.scatter(x = region_sholl_metrics['end_radius'],
                        y = region_sholl_metrics['idxmax_intersections'],
                        color = cmap.to_rgba(region_sholl_metrics['max_intersections']),
                        edgecolor = colors_dict['primecolor'],
                        lw = 0.5,
                        s = 50)
    
    
    for cell_ID in region_sholl_metrics.index:
    
        ax_cv_er_cc.text(x = region_sholl_metrics.at[cell_ID, 'end_radius'],
                          y = region_sholl_metrics.at[cell_ID, 'idxmax_intersections']-2,
                          s = cell_ID,
                          ha = 'center',
                          va = 'top',
                          fontsize = 8,
                          color = colors_dict['primecolor'])
        
    # legend
    ax_cv_er_cc.legend().set_visible(False)
        
    # x axis
    xmin = 0
    xmax = 450
    ax_cv_er_cc.set_xlim([xmin-5, xmax+5])
    ax_cv_er_cc.set_xticks(np.arange(xmin, xmax + 1, 100))
    ax_cv_er_cc.set_xticks(np.arange(xmin, xmax + 1, 25), minor = True)
    ax_cv_er_cc.spines['bottom'].set_bounds([xmin, xmax])
    ax_cv_er_cc.set_xlabel('Enclosing radius [µm]')
    
    # y axis
    ymin = 0
    ymax = 200
    ax_cv_er_cc.set_ylim([ymin-5, ymax+5])
    ax_cv_er_cc.set_yticks(np.arange(ymin, ymax + 1, 50))
    ax_cv_er_cc.set_yticks(np.arange(ymin, ymax + 1, 25), minor = True)
    ax_cv_er_cc.spines['left'].set_bounds([ymin, ymax])
    ax_cv_er_cc.set_ylabel('Critical radius [µm]')
    
    # despine
    [ax_cv_er_cc.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    plt.show()
    
    save_figures(fig_cr_er_cc, f'Enclosing_v_critical_radius-cc_max_inters-{region}', sholl_plots_dir, darkmode_bool)