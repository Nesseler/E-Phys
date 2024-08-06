# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:17:07 2024

@author: nesseler
"""

# general packages
from cellmorphology.cellmorph_colors import neurite_color_dict
from functions.functions_plotting import get_colors, save_figures, get_figure_size
from cellmorphology.cellmorph_packages import mtl, plt, sbn, pd, np, join

# script specific directories / parameters / functions
from parameters.directories_win import cell_morph_descrip_dir, table_file, cell_morph_plots_dir

from cellmorphology.cellmorph_parameters import sholl_step_size


# load files from sholl analysis script
# sholl profiles
sholl_profiles = {'neurites': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_profiles_neurites.xlsx'), index_col='Radius'),
                  'dendrites': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_profiles_dendrites.xlsx'), index_col='Radius'),
                  'axons': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_profiles_axons.xlsx'), index_col='Radius')}


# sholl metrics
sholl_metrics = {'neurites': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_neurites.xlsx'), index_col='cell_ID'),
                 'dendrites': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_dendrites.xlsx'), index_col='cell_ID'),
                 'axons': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_axons.xlsx'), index_col='cell_ID')}


# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# get all cell_IDs
sholl_cell_IDs = sholl_profiles['neurites'].columns.to_list()

# get cell_IDs
cell_IDs_dict = {'all': sholl_cell_IDs,
                 'MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in sholl_cell_IDs],
                 'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in sholl_cell_IDs],
                 'BAOT/MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list() if cell_ID in sholl_cell_IDs]}

# get cell_IDs for cells with axon
cell_IDs_w_axon_dict = {'all': sholl_profiles['axons'].columns.to_list(),
                        'MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in sholl_profiles['axons'].columns.to_list()],
                        'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in sholl_profiles['axons'].columns.to_list()],
                        'BAOT/MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list() if cell_ID in sholl_profiles['axons'].columns.to_list()]}


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
mean_sholl_profiles = pd.DataFrame(columns=mean_columns,
                                   index=sholl_profiles['neurites'].index)
mean_sholl_profiles.index.name = 'Radius'

# standard deviation
std_sholl_profiles = pd.DataFrame(columns=mean_columns,
                                  index=sholl_profiles['neurites'].index)
std_sholl_profiles.index.name = 'Radius'

# fill dataframes
for region in regions:
    for neurite_type in neurite_types:

        if neurite_type == 'axons':
            mean_sholl_profiles['mean_' + region + '-' +
                                neurite_type] = sholl_profiles[neurite_type][cell_IDs_w_axon_dict[region]].mean(axis=1)
            std_sholl_profiles['std_' + region + '-' +
                               neurite_type] = sholl_profiles[neurite_type][cell_IDs_w_axon_dict[region]].std(axis=1)

        else:
            mean_sholl_profiles['mean_' + region + '-' +
                                neurite_type] = sholl_profiles[neurite_type][cell_IDs_dict[region]].mean(axis=1)
            std_sholl_profiles['std_' + region + '-' +
                               neurite_type] = sholl_profiles[neurite_type][cell_IDs_dict[region]].std(axis=1)


# %% initialize plotting

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% create figure of sholl profiles

# initialise figure
fig_sholl, axs_sholl = plt.subplot_mosaic('AB;AB;AB;CD;CE;CF',
                                          figsize=get_figure_size(
                                              height=150, width=150),
                                          layout='tight',
                                          dpi=600)

# flatten nd array of axes
# axs_sholl = axs_sholl.flatten()
axs_keys = {'all': 'A', 'MeA': 'B', 'BAOT': 'C'}

# specify line
line_dict = {'lw': 1, 'alpha': 1}

# plot all profiles together
for region in regions:
    for neurite_type in neurite_types:

        # set axis
        ax = axs_sholl[axs_keys[region]]

        # ax.plot(sholl_profiles[neurite_type],
        #         color = neurite_color_dict[region][neurite_type],
        #         **line_dict,
        #         zorder = 0)

        ax.fill_between(x=sholl_profiles[neurite_type].index.to_list(),
                        y1=mean_sholl_profiles['mean_' + region + '-' + neurite_type] -
                        std_sholl_profiles['std_' +
                                           region + '-' + neurite_type],
                        y2=mean_sholl_profiles['mean_' + region + '-' + neurite_type] +
                        std_sholl_profiles['std_' +
                                           region + '-' + neurite_type],
                        color=neurite_color_dict[region][neurite_type],
                        alpha=0.5,
                        linewidth=0,
                        zorder=1,
                        label='_nolegend_')

        ax.plot(mean_sholl_profiles['mean_' + region + '-' + neurite_type],
                color=neurite_color_dict[region][neurite_type],
                zorder=3,
                label=neurite_type)

# edit sholl profile axes
for region in regions:

    # set axis
    ax = axs_sholl[axs_keys[region]]

    # subplot title
    ax.set_title(region, fontsize=12, loc='left')

    # legend
    ax.legend(prop={'size': 7})

    # y
    ymin = 0
    ymax = 20
    ystep = 5
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# comparisons of neurite types between regions
axs_keys_neurites = {'neurites': 'D', 'dendrites': 'E', 'axons': 'F'}

# add color for 'all' to region_colors
region_colors['all'] = 'gray'

for neurite_type in neurite_types:

    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]

    for region in ['BAOT', 'MeA']:

        ax.set_title(neurite_type, fontsize=12, loc='left')

        ax.plot(mean_sholl_profiles['mean_' + region + '-' + neurite_type],
                # region_colors[region],
                color=neurite_color_dict[region][neurite_type],
                zorder=3,
                label=region)

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
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


for neurite_type in ['axons']:

    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]

    # y
    ymin = 0
    ymax = 4
    ystep = 4
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

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
    ax.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep))
    ax.set_xticks(ticks=np.arange(xmin, xmax + xstepminor, xstepminor), minor=True)
    ax.spines['bottom'].set_bounds([xmin, xmax])


# axis labels
[axs_sholl[key].set_ylabel('Number of intersections [#]') for key in ['A', 'C']]
[axs_sholl[key].set_xlabel('Radius  [µm]') for key in ['C', 'F']]

# align all axis labels
fig_sholl.align_labels()

# show figure
plt.show()

# save figure
save_figures(fig_sholl, 'sholl_profiles_figure', join(cell_morph_plots_dir, 'sholl_plots'),
             darkmode_bool=darkmode_bool, figure_format='both')


# %% sholl profiles figure v2!

# initialise figure
fig_sholl, axs_sholl = plt.subplot_mosaic('BD;BD;BE;CE;CF;CF',
                                          figsize=get_figure_size(height=100, width=150),
                                          layout='tight',
                                          dpi=600)

# flatten nd array of axes
# axs_sholl = axs_sholl.flatten()
axs_keys = {'MeA': 'B', 'BAOT': 'C'}

# specify line
line_dict = {'lw': 1, 'alpha': 1}

# plot all profiles together
for region in ['MeA', 'BAOT']:
    for neurite_type in neurite_types:

        # set axis
        ax = axs_sholl[axs_keys[region]]

        # ax.plot(sholl_profiles[neurite_type],
        #         color = neurite_color_dict[region][neurite_type],
        #         **line_dict,
        #         zorder = 0)

        ax.fill_between(x=sholl_profiles[neurite_type].index.to_list(),
                        y1=mean_sholl_profiles['mean_' + region + '-' + neurite_type] -
                        std_sholl_profiles['std_' +
                                           region + '-' + neurite_type],
                        y2=mean_sholl_profiles['mean_' + region + '-' + neurite_type] +
                        std_sholl_profiles['std_' +
                                           region + '-' + neurite_type],
                        color=neurite_color_dict[region][neurite_type],
                        alpha=0.5,
                        linewidth=0,
                        zorder=1,
                        label='_nolegend_')

        ax.plot(mean_sholl_profiles['mean_' + region + '-' + neurite_type],
                color=neurite_color_dict[region][neurite_type],
                zorder=3,
                label=neurite_type)

# edit sholl profile axes
for region in ['MeA', 'BAOT']:

    # set axis
    ax = axs_sholl[axs_keys[region]]

    # subplot title
    ax.set_title(region, fontsize=12, loc='left')

    # legend
    ax.legend(prop={'size': 7})

    # y
    ymin = 0
    ymax = 20
    ystep = 5
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# comparisons of neurite types between regions
axs_keys_neurites = {'neurites': 'D', 'dendrites': 'E', 'axons': 'F'}

# add color for 'all' to region_colors
region_colors['all'] = 'gray'

for neurite_type in neurite_types:

    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]

    for region in ['BAOT', 'MeA']:

        ax.set_title(neurite_type, fontsize=12, loc='left')

        ax.plot(mean_sholl_profiles['mean_' + region + '-' + neurite_type],
                # region_colors[region],
                color=neurite_color_dict[region][neurite_type],
                zorder=3,
                label=region)

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
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


for neurite_type in ['axons']:

    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]

    # y
    ymin = 0
    ymax = 4
    ystep = 4
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# set all x axis
for key in ['B', 'C', 'D', 'E', 'F']:
    # set axis
    ax = axs_sholl[key]

    # x
    xmin = 0
    xmax = 450
    xstep = 100
    xstepminor = 25
    ax.set_xlim(xmin-10, xmax-10)
    ax.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep))
    ax.set_xticks(ticks=np.arange(
        xmin, xmax + ystepminor, xstepminor), minor=True)
    ax.spines['bottom'].set_bounds([xmin, xmax])


# axis labels
[axs_sholl[key].set_ylabel('Number of\nintersections [#]') for key in ['B', 'C']]
[axs_sholl[key].set_xlabel('Radius  [µm]') for key in ['C', 'F']]

# align all axis labels
fig_sholl.align_labels()

# show figure
plt.show()

# save figure
save_figures(fig_sholl, 'sholl_profiles_figure_v2', join(cell_morph_plots_dir, 'sholl_plots'),
             darkmode_bool=darkmode_bool, figure_format='both')



# %% sholl metrics figure

sholl_metrics_plot_df = pd.DataFrame()

# reorder dataframe for plotting
for neurite_type in neurite_types:

    # get dataframe
    sholl_metrics_pertype = sholl_metrics[neurite_type]

    # set neuritype in dataframe
    sholl_metrics_pertype['neurite_type'] = neurite_type

    # set region in dataframe
    sholl_metrics_pertype['Region'] = MetaData['Region'][sholl_metrics_pertype.index]

    # combine
    sholl_metrics_plot_df = pd.concat(
        [sholl_metrics_plot_df, sholl_metrics_pertype], axis=0)

# remove non categorised cells
sholl_metrics_plot_df = sholl_metrics_plot_df[sholl_metrics_plot_df['Region'] != 'BAOT/MeA']


# %% create figure: violins

# initialise figure
fig_metrics_violins, axs_metrics_violins = plt.subplots(nrows=2,
                                                        ncols=2,
                                                        figsize=get_figure_size(
                                                            width=150, height=150),
                                                        layout='constrained',
                                                        dpi=600)

# flatten axes array
axs_metrics_violins = axs_metrics_violins.flatten()

# create plotting dict for plotting order
order_dict = {1: 'critical_radius',
              2: 'enclosing_radius',
              3: 'max_intersections'}

# create list for subplot alphabetical labels
alpha_labels = ['', 'B', 'C', 'D']

# plot different sholl metrics
for axs_idx, metric in order_dict.items():

    # set axis
    ax = axs_metrics_violins[axs_idx]

    # set label for subplot
    ax.set_title(alpha_labels[axs_idx], fontsize=12, loc='left')

    # violin plots
    violins = sbn.violinplot(data=sholl_metrics_plot_df,
                             x='neurite_type',
                             y=metric,
                             hue='Region',
                             bw=0.4,
                             inner='quart',
                             ax=ax,
                             linewidth=1,
                             palette=region_colors,
                             zorder = 1)

    # edit lines of quarts
    for v_idx in np.arange(len(order_dict.items()) * len(sholl_metrics_plot_df['Region'].drop_duplicates())):

        for l_idx in np.arange(0, 3):
            all_l_idx = v_idx * 3 + l_idx

            if v_idx % 2 == 0:
                violins.lines[all_l_idx].set_color(region_colors['MeA'])
            else:
                violins.lines[all_l_idx].set_color(region_colors['BAOT'])

    # edit color of edges and faces
    for v_idx, violin in enumerate(violins.collections):
        # print(violin)
        if v_idx % 2 == 0:
            violin.set_edgecolor(region_colors['MeA'])
        else:
            violin.set_edgecolor(region_colors['BAOT'])

        # violin.set_edgecolor(colors_dict['primecolor'])
        violin.set_facecolor('None')

    # swarmplots
    swarm = sbn.swarmplot(data=sholl_metrics_plot_df,
                          x='neurite_type',
                          y=metric,
                          hue='Region',
                          ax=ax,
                          size=2,
                          color=colors_dict['primecolor'],
                          dodge=True,
                          zorder = 2)
    
    # errorbar
    for neurite_idx, neurite_type in enumerate(neurite_types):
        for region_x, region in zip([-0.2, +0.2], ['MeA', 'BAOT']):
            
            
            if neurite_type == 'axons':
                cur_cell_IDs_dict = cell_IDs_w_axon_dict
            else:
                cur_cell_IDs_dict = cell_IDs_dict
            
            # data set per type and region for mean and std calculation
            sholl_metric_per_type_n_region = sholl_metrics[neurite_type].loc[cur_cell_IDs_dict[region], :]
            sholl_metric_mean = sholl_metric_per_type_n_region[metric].mean()
            sholl_metric_std = sholl_metric_per_type_n_region[metric].std()
            
            ax.errorbar(x = neurite_idx + region_x,
                        y = sholl_metric_mean,
                        yerr = sholl_metric_std,
                        fmt='_', 
                        markersize = 6,
                        markerfacecolor = 'none',
                        capsize = 2,
                        color=region_colors[region],
                        linewidth = 1,
                        label = '_nolegend_',
                        zorder = 3)
    

    # edit seaborn legend
    ax.legend().set_visible(False)

    # x
    ax.set_xlim(0-0.5, 2+0.5)
    ax.set_xticks(ticks=np.arange(0, 2 + 0.5, 1),
                  labels=neurite_types, rotation=45)
    ax.spines['bottom'].set_bounds([0, 2])
    ax.set_xlabel('')

    # y
    if metric == 'critical_radius':
        ax.set_ylabel('Critical radius [µm]')
        ymin = 0
        ymax = 300
        ystep = 100
        ystepminor = 25

    elif metric == 'enclosing_radius':
        ax.set_ylabel('Enclosing radius [µm]')
        ymin = 0
        ymax = 500
        ystep = 100
        ystepminor = 25

    elif metric == 'max_intersections':
        ax.set_ylabel('Max. number of intersections [#]')
        ymin = 0
        ymax = 50
        ystep = 25
        ystepminor = 5

    # define ypad relative to range
    ypad = (ymax - ymin) * 0.03

    # apply y axis settings
    ax.set_ylim(ymin - ypad, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)
    ax.spines['left'].set_bounds([ymin, ymax])

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# delete first subplot
fig_metrics_violins.delaxes(axs_metrics_violins[0])

# align labels
fig_metrics_violins.align_labels()

# show plot
plt.show()

# save plot
save_figures(fig_metrics_violins, 'sholl_metrics_violins_figure', join(cell_morph_plots_dir, 'sholl_plots'),
             darkmode_bool=darkmode_bool, figure_format='both')


# %% figure for critical vs enclosing radius

# inititalise figure
fig_sholl_scatter, axs_sholl_scatter = plt.subplots(nrows=3,
                                                    ncols=2,
                                                    figsize=get_figure_size(width=150, height=225),
                                                    layout='constrained',
                                                    dpi=600)

# flatten axes array
axs_sholl_scatter = axs_sholl_scatter.flatten()

# create list for subplot alphabetical labels
alpha_labels = ['A', 'B', 'C', 'D', 'E', 'F']

# initialise color-coding of max intersections
cmap_str = 'viridis'

# min max normalize time for color-code
norm = mtl.colors.Normalize(vmin=0,
                            vmax=40)

# create mappable color-map object
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# add colorbar
fig_sholl_scatter.colorbar(cmap, ax=axs_sholl_scatter[-2::],
                           orientation='horizontal',
                           label='Max. number of intersections [#]',
                           fraction=0.05,
                           aspect=70)

# create plotting dict for plotting order
order_dict = {'MeA' : {'neurites' : 0, 'dendrites': 2, 'axons': 4},
              'BAOT': {'neurites' : 1, 'dendrites': 3, 'axons': 5}}

# create dictionary for sorting of inset plots
inset_sort = {'MeA': True, 'BAOT': False}

# create scatter plots splited by region and neurite_type
# loop through region = columns
for region in ['MeA', 'BAOT']:
    
    # loop through neurite_types = rows
    for neurite_type in neurite_types:
        
        # set axis
        ax = axs_sholl_scatter[order_dict[region][neurite_type]]
    
        # subplot title
        ax.set_title(alpha_labels[order_dict[region][neurite_type]] + ': ' + region + ' ' + neurite_type, 
                     fontsize=12, 
                     loc='left')
        
        # set cell_IDs
        if neurite_type == 'dendrites' or neurite_type == 'neurites':
            plt_cell_IDs = cell_IDs_dict
        elif neurite_type == 'axons':
            plt_cell_IDs = cell_IDs_w_axon_dict
        
        # plot scatter plots
        ax.scatter(x=sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['enclosing_radius'],
                   y=sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['critical_radius'],
                   color=cmap.to_rgba(sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['max_intersections']),
                   s=15)
        
        
        # edit main plot axes
        # x
        xmin = 0
        xmax = 500
        xstep = 100
        xstepminor = 25
        xpad = 10
        
        ax.set_xlim(xmin-xpad, xmax)
        ax.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep))
        ax.set_xticks(ticks=np.arange(
            xmin, xmax + ystepminor, xstepminor), minor=True)
        ax.spines['bottom'].set_bounds([xmin, xmax])

        # y
        ymin = 0
        ymax = 400
        ystep = 100
        ystepminor = 25
        ypad = 10
        
        ax.set_ylim(ymin-ypad, ymax)
        ax.set_yticks(ticks=np.arange(ymin, ymax + 1, ystep))
        ax.set_yticks(ticks=np.arange(
            ymin, ymax + ystepminor, ystepminor), minor=True)
        ax.spines['left'].set_bounds([ymin, ymax])
        
        # remove top and right spines
        [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
        
        
        ### add insets

        # ([left, bottom, width, height]), percentages
        ax_inset = ax.inset_axes([0.1, 0.7, 0.25, 0.25])
        
        # plot in inset
        scatter = sbn.scatterplot(data=sholl_metrics_plot_df[sholl_metrics_plot_df['neurite_type'] == neurite_type].sort_values(by=['Region'], ascending=inset_sort[region]),
                                  x='enclosing_radius',
                                  y='critical_radius',
                                  hue='Region',
                                  palette=region_colors,
                                  ax=ax_inset,
                                  s=6,
                                  linewidth=0)
        
        # remove seaborn legend
        ax_inset.legend().set_visible(False)
        
        # edit inset axis
        ax_inset.set_xlim(xmin-xpad, xmax)
        ax_inset.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep), labels=[])
        ax_inset.set_xlabel('')
        ax_inset.tick_params('both', length=1.5, which='major')
        
        ax_inset.set_ylim(ymin-ypad, ymax)
        ax_inset.set_yticks(ticks=np.arange(ymin, ymax + 1, ystep), labels=[])
        ax_inset.set_ylabel('')
        
        # remove top and right spines
        [ax_inset.spines[spine].set_visible(False) for spine in ['top', 'right']]

# axis labels
[ax.set_ylabel('Critical radius [µm]') for ax in axs_sholl_scatter[::2]] 
[ax.set_xlabel('Enclosing radius  [µm]') for ax in axs_sholl_scatter[-2:]]

# align labels
fig_sholl_scatter.align_labels()

# show figure
plt.show()

# save figure
save_figures(fig_sholl_scatter, 'sholl_metrics_crit_v_enclos_figure', join(cell_morph_plots_dir, 'sholl_plots'),
             darkmode_bool = darkmode_bool, figure_format = 'both')

