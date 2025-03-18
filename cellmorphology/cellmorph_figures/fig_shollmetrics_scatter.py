# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:17:07 2024

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# init plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *

# get MetaData
from functions.functions_import import get_MetaData
MetaData = get_MetaData()

# define regions
regions = ['BAOT', 'MeA']

# set neurite types
neurite_types = ['neurites', 'dendrites', 'axons']

# define sholl radius step size
sholl_step_size = 1 # µm


# %% load data

from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_analysis_dir, cellmorph_metrics_dir, cellmorph_shollfigs_dir

# get mean profiles
avg_sholl_profiles = pd.read_excel(join(cellmorph_analysis_dir, 'sholl-combined_profiles', 'sholl_profiles-avg.xlsx'),
                                   index_col = 'Radius')

# get sholl metrics
sholl_metrics = pd.read_excel(join(cellmorph_metrics_dir, 'sholl_metrics.xlsx'),
                              index_col = 'cell_ID')


# %% cell_IDs dicts

# get all cell_IDs
cell_IDs = sholl_metrics.index.to_list()

# get cell_IDs
cell_IDs_dict = {'all': cell_IDs,
                 'MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in cell_IDs],
                 'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in cell_IDs],
                 'BAOT/MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list() if cell_ID in cell_IDs]}

# get cell_IDs for cells with axon
cell_IDs_w_axon_dict = {'all': sholl_metrics['critical_radius-axons'].dropna().index.to_list(),
                        'MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in sholl_metrics['critical_radius-axons'].dropna().index.to_list()],
                        'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in sholl_metrics['critical_radius-axons'].dropna().index.to_list()],
                        'BAOT/MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list() if cell_ID in sholl_metrics['critical_radius-axons'].dropna().index.to_list()]}

sholl_metrics_plot_df = pd.concat([sholl_metrics, MetaData.loc[cell_IDs, 'Region']])

# %% figure for critical vs enclosing radius

# inititalise figure
fig, axs = plt.subplots(nrows=3,
                        ncols=2,
                        figsize=get_figure_size(width = 160, height = 225), #(width=260.334, height=225),
                        layout='constrained',
                        dpi=600)

# flatten axes array
axs = axs.flatten()

# create dict for axis titles 
axs_titles = {'MeA' : {'neurites' : '$\mathregular{A_{i}}$: MeA neurites',
                       'dendrites': '$\mathregular{B_{i}}$: MeA dendrites',
                       'axons'    : '$\mathregular{C_{i}}$: MeA axons'}, 
              'BAOT': {'neurites' : '$\mathregular{A_{ii}}$: BAOT neurites',
                       'dendrites': '$\mathregular{B_{ii}}$: BAOT dendrites',
                       'axons'    : '$\mathregular{C_{ii}}$: BAOT axons'}}

# initialise color-coding of max intersections
cmap_str = 'viridis'

# min max normalize time for color-code
norm = mtl.colors.Normalize(vmin=0,
                            vmax=40)

# create mappable color-map object
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# add colorbar
fig.colorbar(cmap, ax=axs[-2::],
             orientation = 'horizontal',
             label = 'Max. number of intersections [#]',
             fraction = 0.05,
             aspect = 70)

# create plotting dict for plotting order
order_dict = {'MeA' : {'neurites' : 0, 'dendrites': 2, 'axons': 4},
              'BAOT': {'neurites' : 1, 'dendrites': 3, 'axons': 5}}

# create dictionary for sorting of inset plots
inset_sort = {'MeA': True, 'BAOT': False}

# create scatter plots splited by region and neurite_type
# loop through region = columns
for region in ['MeA', 'BAOT']:
    
    # get the other region
    other_region = [r for r in regions if r != region][0]
    
    # loop through neurite_types = rows
    for ntype in neurite_types:
        
        # set axis
        ax = axs[order_dict[region][ntype]]
    
        # subplot title
        ax.set_title(axs_titles[region][ntype], 
                     fontsize=12, 
                     loc='left')
        
        # # set cell_IDs
        # if ntype == 'dendrites' or ntype == 'neurites':
        #     plt_cell_IDs = cell_IDs_dict
        # elif ntype == 'axons':
        #     plt_cell_IDs = cell_IDs_w_axon_dict
            
        # plot diagonal line in plot
        ax.plot([0, 400], [0, 400],
                color = colors_dict['primecolor'],
                linewidth = 1,
                linestyle = 'dashed',
                alpha = 0.5, 
                zorder = 0)
        
        # plot scatter plots
        ax.scatter(x = sholl_metrics.loc[cell_IDs_dict[region], f'enclosing_radius-{ntype}'],
                   y = sholl_metrics.loc[cell_IDs_dict[region], f'critical_radius-{ntype}'],
                   color = cmap.to_rgba(sholl_metrics.loc[cell_IDs_dict[region], f'max_intersections-{ntype}']),
                   s = 15,
                   zorder = 1)
        
        
        # edit main plot axes
        # x
        xdict = {'ax_min' : 0,
                 'ax_max' : 500,
                 'pad' : None,
                 'step' : 100,
                 'stepminor' : 25,
                 'label' : 'Enclosing radius [µm]'}
        
        apply_axis_settings(ax, axis = 'x', **xdict)

        # y
        ydict = {'ax_min' : 0,
                 'ax_max' : 400,
                 'pad' : None,
                 'step' : 100,
                 'stepminor' : 25,
                 'label' : 'Enclosing radius [µm]'}
        
        apply_axis_settings(ax, axis = 'y', **ydict)
        
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
        # scatter = sbn.scatterplot(data=sholl_metrics_plot_df.loc[cell_IDs_dict[region], :]),
        #                           x = f'enclosing_radius-{ntype}',
        #                           y = f'critical_radius-{ntype}',
        #                           hue = 'Region',
        #                           palette  =region_colors,
        #                           ax=ax_inset,
        #                           s=6,
        #                           linewidth=0, 
        #                           zorder = 2)
    
        scatter = sbn.scatterplot(data=sholl_metrics_plot_df.loc[cell_IDs_dict[other_region], :]),
                                  x = f'enclosing_radius-{ntype}',
                                  y = f'critical_radius-{ntype}',
                                  hue = 'Region',
                                  palette = region_colors,
                                  ax = ax_inset,
                                  s = 6,
                                  linewidth = 0, 
                                  zorder = 1,
                                  alpha = 0.5)
        
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
[ax.set_ylabel('Critical radius [µm]') for ax in axs[::2]] 
[ax.set_xlabel('Enclosing radius  [µm]') for ax in axs]

# align labels
fig.align_labels()

# show figure
plt.show()

# save figure
# crit_v_enclos_figure_dir = join(cell_morph_plots_dir, 'figure-sholl_metrics-enclosing_v_critical_radius')

# save_figures(fig_sholl_scatter, 'sholl_metrics-crit_v_enclos-figure', 
#              save_dir = crit_v_enclos_figure_dir,
#              darkmode_bool = darkmode_bool, 
#              figure_format = 'both')


# %% collected data to save

