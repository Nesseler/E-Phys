# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 15:22:34 2024

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


# %% drop cells E-126 & E-158

from cellmorphology.cellmorph_parameters import cell_IDs_toDrop

for neurite_type in sholl_metrics.keys():
    
    # iterate through cell_IDs to drop to check wether their contain in the dataframes
    for cell_ID_toDrop in cell_IDs_toDrop:
        
        # test if cell_ID is in index
        if cell_ID_toDrop in sholl_metrics[neurite_type].index:
            sholl_metrics[neurite_type].drop(index = cell_ID_toDrop, inplace = True)
            
        # test if cell_ID is in columns
        if cell_ID_toDrop in sholl_profiles[neurite_type].columns:
            sholl_profiles[neurite_type].drop(columns = cell_ID_toDrop, inplace = True)


# %%

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
std_columns = ['std_' + col for col in measurement_columns]

# average
mean_sholl_profiles = pd.DataFrame(columns=mean_columns,
                                   index=sholl_profiles['neurites'].index)
mean_sholl_profiles.index.name = 'Radius'

# standard deviation
std_sholl_profiles = pd.DataFrame(columns=std_columns,
                                  index=sholl_profiles['neurites'].index)
std_sholl_profiles.index.name = 'Radius'

# fill dataframes
for region in regions:
    for neurite_type in neurite_types:

        if neurite_type == 'axons':
            mean_sholl_profiles['mean_' + region + '-' + neurite_type] = sholl_profiles[neurite_type][cell_IDs_w_axon_dict[region]].mean(axis=1)
            std_sholl_profiles['std_' + region + '-' + neurite_type] = sholl_profiles[neurite_type][cell_IDs_w_axon_dict[region]].std(axis=1)

        else:
            mean_sholl_profiles['mean_' + region + '-' + neurite_type] = sholl_profiles[neurite_type][cell_IDs_dict[region]].mean(axis=1)
            std_sholl_profiles['std_' + region + '-' + neurite_type] = sholl_profiles[neurite_type][cell_IDs_dict[region]].std(axis=1)


# %% initialize plotting

# set colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

# set color dict for neurite_types and regions
neurites_regions_color_dict = {'all'  : {'dendrites' : 'grey',                'axons' : 'lightgrey'},
                               'MeA'  : {'dendrites' : region_colors['MeA'],  'axons' : colors_dict['MeA_lighter']},
                               'BAOT' : {'dendrites' : region_colors['BAOT'], 'axons' : colors_dict['BAOT_lighter']}}

# set font size
mtl.rcParams.update({'font.size': 12})



# %% sholl profiles figure v2!

# initialise figure
fig_sholl, axs_sholl = plt.subplots(nrows = 2,
                                    ncols = 1,
                                    figsize=get_figure_size(height=100, width=87.48),
                                    layout='constrained',
                                    dpi=600)

# flatten nd array of axes
# axs_sholl = axs_sholl.flatten()
axs_keys = {'MeA': 0, 'BAOT': 1}

# specify line
line_dict = {'lw': 1, 'alpha': 1}

# set labels for subplots
axs_titles = {'MeA' : 'MeA',
              'BAOT': 'BAOT'}

# plot all profiles together
for region in ['MeA', 'BAOT']:
    
    for neurite_type in ['dendrites', 'axons']:

        # set axis
        ax = axs_sholl[axs_keys[region]]

        # ax.plot(sholl_profiles[neurite_type],
        #         color = neurite_color_dict[region][neurite_type],
        #         **line_dict,
        #         zorder = 0)

        ax.fill_between(x = sholl_profiles[neurite_type].index.to_list(),
                        y1 = mean_sholl_profiles['mean_' + region + '-' + neurite_type] - std_sholl_profiles['std_' + region + '-' + neurite_type],
                        y2 = mean_sholl_profiles['mean_' + region + '-' + neurite_type] + std_sholl_profiles['std_' + region + '-' + neurite_type],
                        color = neurites_regions_color_dict[region][neurite_type],
                        alpha=0.5,
                        linewidth=0,
                        zorder=1,
                        label='_nolegend_')

        ax.plot(mean_sholl_profiles['mean_' + region + '-' + neurite_type],
                color = neurites_regions_color_dict[region][neurite_type],
                zorder=3,
                label=neurite_type.title())

# edit sholl profile axes
for region in ['MeA', 'BAOT']:

    # set axis
    ax = axs_sholl[axs_keys[region]]

    # subplot title
    ax.set_title(axs_titles[region], fontsize=12, loc='left')

    # legend
    ax.legend(prop={'size': 9})

    # y
    ymin = 0
    ymax = 20
    ystep = 10
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]



# set all x axis
for ax in axs_sholl:

    # x
    xmin = 0
    xmax = 400
    xstep = 100
    xstepminor = 25
    ax.set_xlim(xmin-10, xmax-10)
    ax.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep))
    ax.set_xticks(ticks=np.arange(
        xmin, xmax + ystepminor, xstepminor), minor=True)
    ax.spines['bottom'].set_bounds([xmin, xmax])


# axis labels
[ax.set_ylabel('') for ax in axs_sholl] 
[ax.set_xlabel('') for ax in axs_sholl]

# set axis title for figure
fig_sholl.supylabel('Number of intersections [#]', fontsize = 12)
fig_sholl.supxlabel('Radius  [µm]', fontsize = 12)


# align all axis labels
fig_sholl.align_labels()

# show figure
plt.show()

# save figure
sholl_profiles_fig_dir = "C:/Users/nesseler/Desktop/Poster_iBehave"
save_figures(fig_sholl, 
              figure_name = 'sholl_profiles_figure_v2', 
              save_dir = sholl_profiles_fig_dir,
              darkmode_bool=darkmode_bool, 
              figure_format='both')


# %% collect data and write to dataframe to save

# concat mean and std of sholl profiles 
mean_std_sholl_profiles = pd.concat([mean_sholl_profiles, std_sholl_profiles], axis = 1)

# save dataframe
mean_std_sholl_profiles.to_excel(join(sholl_profiles_fig_dir, 'sholl_profiles-means_n_stds.xlsx'), index_label = 'Radius')



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



# %% figure for critical vs enclosing radius

# inititalise figure
fig_sholl_scatter, axs_sholl_scatter = plt.subplots(nrows=2,
                                                    ncols=2,
                                                    figsize=get_figure_size(width=135, height=100),
                                                    layout='constrained',
                                                    dpi=600)


# adjust constrained layout
fig_sholl_scatter.get_layout_engine().set(wspace = 0.05, hspace = 0.1)

# flatten axes array
axs_sholl_scatter = axs_sholl_scatter.flatten()

# create list for subplot alphabetical labels
# alpha_labels = ['A', 'B', 'C', 'D', 'E', 'F']

# create dict for axis titles 
axs_titles = {'MeA' : {'neurites' : 'MeA neurites',
                       'dendrites': 'MeA dendrites',
                       'axons'    : 'MeA axons'}, 
              'BAOT': {'neurites' : 'BAOT neurites',
                       'dendrites': 'BAOT dendrites',
                       'axons'    : 'BAOT axons'}}

# initialise color-coding of max intersections
cmap_str = 'viridis'

# min max normalize time for color-code
norm = mtl.colors.Normalize(vmin=0,
                            vmax=40)

# create mappable color-map object
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# add colorbar
cbar = fig_sholl_scatter.colorbar(cmap, ax=axs_sholl_scatter[::],
                                  orientation='vertical',
                                  label='Max. number of intersections [#]',
                                  fraction=0.1,
                                  aspect=60)

# adjust ticks on colorbar
cbar.set_ticks(np.arange(0, 40 + 1, 10))

# create plotting dict for plotting order
order_dict = {'MeA' : {'neurites' : 0, 'dendrites': 0, 'axons': 1},
              'BAOT': {'neurites' : 1, 'dendrites': 2, 'axons': 3}}

# create dictionary for sorting of inset plots
inset_sort = {'MeA': True, 'BAOT': False}

# create scatter plots splited by region and neurite_type
# loop through region = columns
for region in ['MeA', 'BAOT']:
    
    # loop through neurite_types = rows
    for neurite_type in ['dendrites', 'axons']:
        
        # set axis
        ax = axs_sholl_scatter[order_dict[region][neurite_type]]
    
        # subplot title
        ax.set_title(axs_titles[region][neurite_type], 
                     fontsize=12, 
                     loc='left')
        
        # set cell_IDs
        if neurite_type == 'dendrites' or neurite_type == 'neurites':
            plt_cell_IDs = cell_IDs_dict
        elif neurite_type == 'axons':
            plt_cell_IDs = cell_IDs_w_axon_dict
            
        # plot diagonal line in plot
        ax.plot([0, 400], [0, 400],
                color = colors_dict['primecolor'],
                linewidth = 1,
                linestyle = 'dashed',
                alpha = 0.5, 
                zorder = 0)
        
        # plot scatter plots
        ax.scatter(x=sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['enclosing_radius'],
                   y=sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['critical_radius'],
                   color=cmap.to_rgba(sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['max_intersections']),
                   s=10,
                   zorder = 1)
        
        
        # edit main plot axes
        # x
        xmin = 0
        xmax = 500
        xstep = 250
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
        ystep = 200
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
        ax_inset.plot([0, 400], [0, 400],
                      color = colors_dict['primecolor'],
                      linewidth = 0.5,
                      linestyle = 'dashed',
                      alpha = 0.5, 
                      zorder = 0)
        
        
        # scatterplot in inset
        scatter = sbn.scatterplot(data=sholl_metrics_plot_df[sholl_metrics_plot_df['neurite_type'] == neurite_type].sort_values(by=['Region'], ascending=inset_sort[region]),
                                  x='enclosing_radius',
                                  y='critical_radius',
                                  hue='Region',
                                  palette=region_colors,
                                  ax=ax_inset,
                                  s=2,
                                  linewidth=0, 
                                  zorder = 1)
        
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
[ax.set_ylabel('') for ax in axs_sholl_scatter[::2]] 
[ax.set_xlabel('') for ax in axs_sholl_scatter]

# set axis title for figure
fig_sholl_scatter.supylabel('Critical radius [µm]', fontsize = 12)
fig_sholl_scatter.supxlabel('Enclosing radius  [µm]', fontsize = 12)

# align labels
fig_sholl_scatter.align_labels()

# show figure
plt.show()

# save figure
crit_v_enclos_figure_dir = "C:/Users/nesseler/Desktop/Poster_iBehave"

save_figures(fig_sholl_scatter, 'sholl_metrics-crit_v_enclos-figure', 
              save_dir = crit_v_enclos_figure_dir,
              darkmode_bool = darkmode_bool, 
              figure_format = 'both')

# %% collected data to save

# create dataframe
sholl_metrics_df = pd.DataFrame(index = sholl_cell_IDs)

# iterate through neurite_type
for neurite_type in neurite_types:
    
    # iterate through regions
    for region in ['MeA', 'BAOT']:
        
        # get sholl metric per Type
        sholl_metric_perType = sholl_metrics[neurite_type]

        # limit to current region
        sholl_metric_perType_n_region = sholl_metric_perType[sholl_metric_perType['Region'] == region]
        
        # get current cell_IDs
        sholl_metric_region_cell_IDs = sholl_metric_perType_n_region.index.to_list()
        
        # iterate through metrics
        for sholl_metric in ['enclosing_radius', 'critical_radius', 'max_intersections']:

            # write to dataframe
            sholl_metrics_df.loc[sholl_metric_region_cell_IDs, f'{sholl_metric}-{neurite_type}-{region}'] = sholl_metric_perType_n_region[sholl_metric]

# save dataframe
# sholl_metrics_df.to_excel(join(crit_v_enclos_figure_dir, 'sholl_metrics_sorted.xlsx'), index_label = 'cell_ID')