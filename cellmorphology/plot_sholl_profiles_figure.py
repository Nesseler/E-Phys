# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:17:07 2024

@author: nesseler
"""

# general packages
from cellmorphology.cellmorph_packages import *

# script specific directories / parameters / functions
from parameters.directories_win import cell_morph_descrip_dir, table_file, cell_morph_plots_dir

from cellmorphology.cellmorph_parameters import sholl_step_size


# load files from sholl analysis script
## sholl profiles
sholl_profiles = {'neurites' : pd.read_excel(join(cell_morph_descrip_dir, 'sholl_profiles_neurites.xlsx'), index_col= 'Radius'),
                  'dendrites' : pd.read_excel(join(cell_morph_descrip_dir, 'sholl_profiles_dendrites.xlsx'), index_col= 'Radius'),
                  'axons' : pd.read_excel(join(cell_morph_descrip_dir, 'sholl_profiles_axons.xlsx'), index_col= 'Radius')}


## sholl metrics
sholl_metrics = {'neurites' : pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_neurites.xlsx'), index_col= 'cell_ID'),
                 'dendrites' : pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_dendrites.xlsx'), index_col= 'cell_ID'),
                 'axons' : pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_axons.xlsx'), index_col= 'cell_ID')}


# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# get all cell_IDs
sholl_cell_IDs = sholl_profiles['neurites'].columns.to_list()

# get cell_IDs
cell_IDs_dict = {'all' : sholl_cell_IDs,
                 'MeA' : [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in sholl_cell_IDs],
                 'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in sholl_cell_IDs],
                 'BAOT/MeA' : [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list() if cell_ID in sholl_cell_IDs]}

# get cell_IDs for cells with axon
cell_IDs_w_axon_dict = {'all' : sholl_profiles['axons'].columns.to_list(),
                        'MeA' : [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in sholl_profiles['axons'].columns.to_list()],
                        'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in sholl_profiles['axons'].columns.to_list()],
                        'BAOT/MeA' : [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list() if cell_ID in sholl_profiles['axons'].columns.to_list()]}


# %% calc means, stds, etc...

# initialise dataframe
neurite_types = ['neurites', 'dendrites', 'axons']
regions = ['all', 'MeA', 'BAOT']
measurement_columns = []

for region in regions:
    for neurite_type in neurite_types:
        measurement_columns.append(region + '-' + neurite_type)
        
mean_columns = ['mean_' + col for col in measurement_columns]    
std_columns = ['mean_' + col for col in measurement_columns]

# average
mean_sholl_profiles = pd.DataFrame(columns = mean_columns,
                                   index = sholl_profiles['neurites'].index)
mean_sholl_profiles.index.name = 'Radius'

# standard deviation
std_sholl_profiles = pd.DataFrame(columns = mean_columns,
                                  index = sholl_profiles['neurites'].index)
std_sholl_profiles.index.name = 'Radius'

# fill dataframes
for region in regions:
    for neurite_type in neurite_types:
        
        if neurite_type == 'axons':
            mean_sholl_profiles['mean_' + region + '-' + neurite_type] = sholl_profiles[neurite_type][cell_IDs_w_axon_dict[region]].mean(axis = 1)
            std_sholl_profiles['std_' + region + '-' + neurite_type] = sholl_profiles[neurite_type][cell_IDs_w_axon_dict[region]].std(axis = 1)
        
        else:
            mean_sholl_profiles['mean_' + region + '-' + neurite_type] = sholl_profiles[neurite_type][cell_IDs_dict[region]].mean(axis = 1)
            std_sholl_profiles['std_' + region + '-' + neurite_type] = sholl_profiles[neurite_type][cell_IDs_dict[region]].std(axis = 1)

       



# %% initialize plotting

from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)
from cellmorphology.cellmorph_colors import neurite_color_dict

# set font size 
mtl.rcParams.update({'font.size': 9})


# %% create figure of sholl profiles

# initialise figure
fig_sholl, axs_sholl = plt.subplot_mosaic('AB;AB;AB;CD;CE;CF',
                                          figsize = get_figure_size(height = 150, width = 150),
                                          layout = 'constrained',
                                          dpi = 600)

# flatten nd array of axes
# axs_sholl = axs_sholl.flatten()
axs_keys = {'all' : 'A', 'MeA' : 'B', 'BAOT' : 'C'}

# specify line
line_dict = {'lw' : 1, 'alpha' : 1}

### plot all profiles together
for region in regions:
    for neurite_type in neurite_types:
        
        # set axis
        ax = axs_sholl[axs_keys[region]]
        
        # ax.plot(sholl_profiles[neurite_type], 
        #         color = neurite_color_dict[region][neurite_type],
        #         **line_dict,
        #         zorder = 0)
    
        ax.fill_between(x = sholl_profiles[neurite_type].index.to_list(),
                        y1 = mean_sholl_profiles['mean_' + region + '-' + neurite_type] - std_sholl_profiles['std_' + region + '-' + neurite_type],
                        y2 = mean_sholl_profiles['mean_' + region + '-' + neurite_type] + std_sholl_profiles['std_' + region + '-' + neurite_type],
                        color = neurite_color_dict[region][neurite_type],
                        alpha = 0.5,
                        linewidth = 0,
                        zorder = 1,
                        label = '_nolegend_')
    
        ax.plot(mean_sholl_profiles['mean_' + region + '-' + neurite_type],
                color = neurite_color_dict[region][neurite_type],
                zorder = 3,
                label = neurite_type)

# edit sholl profile axes
for region in regions:
    
    # set axis
    ax = axs_sholl[axs_keys[region]]
    
    # subplot title
    ax.set_title(region, fontsize = 9, loc = 'left')
    
    # legend
    ax.legend(prop={'size': 7})

    # y
    ymin = 0
    ymax = 20
    ystep = 5
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks = np.arange(ymin, ymax +ystep , ystep))
    ax.set_yticks(ticks = np.arange(ymin, ymax +ystepminor , ystepminor), minor = True)
    
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


### comparisons of neurite types between regions
axs_keys_neurites = {'neurites' : 'D', 'dendrites' : 'E', 'axons' : 'F'}

# add color for 'all' to region_colors
region_colors['all'] = 'gray'

for neurite_type in neurite_types:
    
    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]
    
    for region in ['BAOT', 'MeA']:
        
        ax.set_title(neurite_type, fontsize = 9, loc = 'left')
        
        ax.plot(mean_sholl_profiles['mean_' + region + '-' + neurite_type],
                color = neurite_color_dict[region][neurite_type], #region_colors[region],
                zorder = 3,
                label = region)
        
    # legend
    ax.legend(prop={'size': 7})


# configure axis of neurites and dendrites
for neurite_type in ['neurites', 'dendrites']:
    
    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]  
    
    # y
    ymin = 0
    ymax = 15
    ystep = 15
    ystepminor = 5
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks = np.arange(ymin, ymax +ystep , ystep))
    ax.set_yticks(ticks = np.arange(ymin, ymax +ystepminor , ystepminor), minor = True)
    
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


for neurite_type in ['axons']:
    
    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]  
    
    # y
    ymin = 0
    ymax = 2
    ystep = 2
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks = np.arange(ymin, ymax +ystep , ystep))
    ax.set_yticks(ticks = np.arange(ymin, ymax +ystepminor , ystepminor), minor = True)
    
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# set all x axis
for key in ['A', 'B', 'C', 'D', 'E', 'F']:
    # set axis
    ax = axs_sholl[key]
    
    # x
    xmin = 0
    xmax = 450
    xstep = 100
    xstepminor = 25
    ax.set_xlim(xmin-10, xmax-10)
    ax.set_xticks(ticks = np.arange(xmin, xmax +1, xstep))
    ax.set_xticks(ticks = np.arange(xmin, xmax +ystepminor , xstepminor), minor = True)
    ax.spines['bottom'].set_bounds([xmin, xmax])


# axis labels
[axs_sholl[key].set_ylabel('Number of intersections [#]') for key in ['A', 'C']]
[axs_sholl[key].set_xlabel('Radius  [µm]') for key in ['C', 'F']]

fig_sholl.align_labels()


plt.show()


save_figures(fig_sholl, 'sholl_profiles_figure', join(cell_morph_plots_dir, 'sholl_plots'),
             darkmode_bool = darkmode_bool, figure_format = 'both')


# %% sholl metrics figure

sholl_metrics_plot_df = pd.DataFrame()

# reorder dataframe for plotting
for neurite_type in ['dendrites', 'axons']:
    
    # get dataframe
    sholl_metrics_pertype = sholl_metrics[neurite_type]
    
    # set neuritype in dataframe
    sholl_metrics_pertype['neurite_type'] = neurite_type
    
    # set region in dataframe
    sholl_metrics_pertype['Region'] = MetaData['Region'][sholl_metrics_pertype.index]
    
    # combine 
    sholl_metrics_plot_df = pd.concat([sholl_metrics_plot_df, sholl_metrics_pertype], axis = 0)
   
# remove non categorised cells
sholl_metrics_plot_df = sholl_metrics_plot_df[sholl_metrics_plot_df['Region'] != 'BAOT/MeA']
    

# %% create figure

# initialise figure
fig_metrics, axs_metrics = plt.subplot_mosaic('AABBCC;DDDEEE;FFFGGG',
                                          figsize = get_figure_size(width = 150, height = 225),
                                          layout = 'constrained',
                                          dpi = 600)




for axs_key, metric in {'A' : 'critical_radius', 'B' : 'enclosing_radius', 'C' : 'max_intersections'}.items():
    
    # set axis
    ax = axs_metrics[axs_key]

    # violin plots
    violin = sbn.violinplot(data = sholl_metrics_plot_df,
                            x = 'neurite_type',
                            y = metric,
                            hue = 'Region',
                            #split = True,
                            inner = 'quart',
                            ax = ax,
                            linewidth = 1,
                            palette = region_colors,
                            gap = 0.1)
    
    for l in violin.lines:
        l.set_color(colors_dict['primecolor'])
    
    for violin in violin.collections:
        violin.set_edgecolor(colors_dict['primecolor'])
        # violin.set_facecolor('None')
        
    # swarmplots
    swarm = sbn.swarmplot(data = sholl_metrics_plot_df,
                          x = 'neurite_type',
                          y = metric,
                          hue = 'Region',
                          ax = ax,
                          size = 1.5,
                          color = colors_dict['primecolor'],
                          dodge = True)
    
    # edit seaborn legend
    ax.legend().set_visible(False)
    
    # x
    ax.set_xlim(0-0.5, 1+0.5)
    ax.set_xticks(ticks = np.arange(0, 1 +0.5, 1), labels = ['dendrites', 'axons'], rotation = 45)
    ax.spines['bottom'].set_bounds([0, 1])
    ax.set_xlabel('')
    
    # y
    if metric == 'critical_radius':
        ax.set_ylabel('Critical radius [µm]')
        ymin = 0
        ymax = 300
        ystep = 100
        ystepminor = 25
        ypad = 10
        
    elif metric == 'enclosing_radius':
        ax.set_ylabel('Enclosing radius [µm]')
        ymin = 0
        ymax = 500
        ystep = 100
        ystepminor = 25
        ypad = 10
        
    elif metric == 'max_intersections':
        ax.set_ylabel('Max. number of intersections [#]')
        ymin = 0
        ymax = 50
        ystep = 25
        ystepminor = 5
        ypad = 5
        
        
    ax.set_ylim(ymin -ypad, ymax)
    ax.set_yticks(ticks = np.arange(ymin, ymax +ystep , ystep))
    ax.set_yticks(ticks = np.arange(ymin, ymax +ystepminor , ystepminor), minor = True)
    ax.spines['left'].set_bounds([ymin, ymax])
    
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]



# critical vs enclosing vs max_intersection (color-coding)


# initialise color-coding of max intersections
# specify color map
cmap = 'magma'

# min max normalize time for color-code
norm = mtl.colors.Normalize(vmin = 0, 
                            vmax = 40)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap)

# colorbar
fig_metrics.colorbar(cmap, ax = [axs_metrics['F'], axs_metrics['G']],
                     orientation = 'horizontal', 
                     label = 'Max. number of intersections [#]',
                     fraction = 0.05,
                     aspect = 70)


ax = axs_metrics['D']


axs_keys = {'MeA' : {'dendrites' : 'D', 'axons' : 'F'},
            'BAOT' : {'dendrites' : 'E', 'axons' : 'G'}}

## ([left, bottom, width, height]), percentages
inset_sort = {'MeA' : True,
              'BAOT' : False}


for region in ['MeA', 'BAOT']:
    
    for neurite_type in ['dendrites', 'axons']:

        # set axis
        ax = axs_metrics[axs_keys[region][neurite_type]]
        
        # subplot title
        ax.set_title(region + ' ' + neurite_type, fontsize = 9, loc = 'left')
        
        # set cell_IDs
        if neurite_type == 'dendrites':
            plt_cell_IDs = cell_IDs_dict
        elif neurite_type == 'axons':
            plt_cell_IDs = cell_IDs_w_axon_dict

        ax.scatter(x = sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['enclosing_radius'],
                   y = sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['critical_radius'],
                   color = cmap.to_rgba(sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['max_intersections']),
                   s = 10)
        
        
        
        # edit seaborn legend
        ax.legend().set_visible(False)
        
        
        
        # x
        xmin = 0
        xmax = 500
        xstep = 100
        xstepminor = 25
        ax.set_xlim(xmin, xmax)
        ax.set_xticks(ticks = np.arange(xmin, xmax +1, xstep))
        ax.set_xticks(ticks = np.arange(xmin, xmax +ystepminor , xstepminor), minor = True)
        ax.spines['bottom'].set_bounds([xmin, xmax])
        
        # y
        ymin = 0
        ymax = 400
        ystep = 100
        ystepminor = 25
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(ticks = np.arange(ymin, ymax +1, ystep))
        ax.set_yticks(ticks = np.arange(ymin, ymax +ystepminor , ystepminor), minor = True)
        ax.spines['left'].set_bounds([ymin, ymax])
        
        
        
        
        
        
        # add inset
        ## ([left, bottom, width, height]), percentages
        ax_inset = ax.inset_axes([0.1,0.65,0.25,0.25])
        
        # plot in inset
        scatter = sbn.scatterplot(data = sholl_metrics_plot_df[sholl_metrics_plot_df['neurite_type'] == neurite_type].sort_values(by = ['Region'], ascending = inset_sort[region]),
                                  x = 'enclosing_radius',
                                  y = 'critical_radius',
                                  hue = 'Region', 
                                  palette = region_colors,
                                  ax = ax_inset,
                                  s = 5,
                                  linewidth = 0)
        
        
        ax_inset.legend().set_visible(False)
        
        ax_inset.set_xlim(xmin, xmax)
        ax_inset.set_xticks(ticks = np.arange(xmin, xmax +1, xstep), labels = [])
        ax_inset.set_xlabel('')
        ax_inset.tick_params('both', length=1.5, which='major')
        
        ax_inset.set_ylim(ymin, ymax)
        ax_inset.set_yticks(ticks = np.arange(ymin, ymax +1, ystep), labels = [])
        ax_inset.set_ylabel('')
        
        
        
# axis labels
[axs_metrics[key].set_ylabel('Critical radius [µm]') for key in ['D', 'F']]
[axs_metrics[key].set_xlabel('Enclosing radius  [µm]') for key in ['F', 'G']]

fig_metrics.align_labels()


plt.show()

save_figures(fig_metrics, 'sholl_metrics_figure', join(cell_morph_plots_dir, 'sholl_plots'),
             darkmode_bool = darkmode_bool, figure_format = 'both')


