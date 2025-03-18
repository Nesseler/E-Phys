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


# %% sholl profiles figure 

# initialise figure
fig_sholl, axs_sholl = plt.subplot_mosaic('BBBCCC;DDEEFF',#'BD;BD;BE;CE;CF;CF',
                                          figsize=get_figure_size(height=125.103, width=260.334),
                                          layout='tight',
                                          height_ratios = [2.5,1],
                                          dpi=600)

# flatten nd array of axes
# axs_sholl = axs_sholl.flatten()
axs_keys = {'MeA': 'B', 'BAOT': 'C'}

# specify line
line_dict = {'lw': 1, 'alpha': 1}

# set labels for subplots
axs_titles = {'MeA' : 'A: MeA',
              'BAOT': 'B: BAOT',
              'neurites' : '$\mathregular{C_{i}}$: Neurites',
              'dendrites' : '$\mathregular{C_{ii}}$: Dendrites',
              'axons' : '$\mathregular{C_{iii}}$: Axons'}

# plot all profiles together
for region in ['MeA', 'BAOT']:
    for neurite_type in ['dendrites', 'axons']:

        # set axis
        ax = axs_sholl[axs_keys[region]]
        
        # calc upper and lower border
        y_bottom = avg_sholl_profiles[f'mean-sholl_profile-{neurite_type}-{region}'] - avg_sholl_profiles[f'std-sholl_profile-{neurite_type}-{region}']
        y_top = avg_sholl_profiles[f'mean-sholl_profile-{neurite_type}-{region}'] + avg_sholl_profiles[f'std-sholl_profile-{neurite_type}-{region}']

        # plot std as shade
        ax.fill_between(x = avg_sholl_profiles.index.to_list(),
                        y1 = y_bottom,
                        y2 = y_top,
                        color = neurite_color_dict[region][neurite_type],
                        alpha = 0.5,
                        linewidth = 0,
                        zorder = 1,
                        label = '_nolegend_')

        # plot mean
        ax.plot(avg_sholl_profiles[f'mean-sholl_profile-{neurite_type}-{region}'],
                color = neurite_color_dict[region][neurite_type],
                zorder = 3,
                label = neurite_type.title())

# edit sholl profile axes
for region in ['MeA', 'BAOT']:

    # set axis
    ax = axs_sholl[axs_keys[region]]

    # subplot title
    ax.set_title(axs_titles[region], fontsize=12, loc='left')

    # legend
    ax.legend(fontsize = 9,
              frameon = False,
              title = 'Neurite types',
              title_fontsize = 9)

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

        ax.set_title(axs_titles[neurite_type], fontsize=12, loc='left')

        ax.plot(avg_sholl_profiles[f'mean-sholl_profile-{neurite_type}-{region}'],
                color = neurite_color_dict[region][neurite_type],
                zorder = 3,
                label  =region)

    # legend
    ax.legend(fontsize = 9,
              frameon = False,
              title = 'Region',
              title_fontsize = 9)


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
    xmax = 400
    xstep = 100
    xstepminor = 25
    ax.set_xlim(xmin-10, xmax-10)
    ax.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep))
    ax.set_xticks(ticks=np.arange(
        xmin, xmax + ystepminor, xstepminor), minor=True)
    ax.spines['bottom'].set_bounds([xmin, xmax])


# axis labels
[axs_sholl[key].set_ylabel('Number of\nintersections [#]') for key in ['B', 'C', 'D']]
[axs_sholl[key].set_xlabel('Radius  [µm]') for key in ['B', 'C', 'D','E', 'F']]

# align all axis labels
fig_sholl.align_labels()

# show figure
plt.show()

# save figure
# sholl_profiles_fig_dir = join(cellmorph_shollfigs_dir, 'figure-sholl_profiles')
save_figures(fig_sholl, 
              figure_name = 'sholl_profiles_figure', 
              save_dir = cellmorph_shollfigs_dir,
              darkmode_bool=darkmode_bool, 
              figure_format='both')


# # %% sholl metrics figure

# sholl_metrics_plot_df = pd.DataFrame()

# # # reorder dataframe for plotting
# # for neurite_type in neurite_types:

# #     # get dataframe
# #     sholl_metrics_pertype = sholl_metrics[neurite_type]

# #     # set neuritype in dataframe
# #     sholl_metrics_pertype['neurite_type'] = neurite_type

# #     # set region in dataframe
# #     sholl_metrics_pertype['Region'] = MetaData['Region'][sholl_metrics_pertype.index]

# #     # combine
# #     sholl_metrics_plot_df = pd.concat(
# #         [sholl_metrics_plot_df, sholl_metrics_pertype], axis=0)

# # remove non categorised cells
# sholl_metrics_plot_df = pd.concat([sholl_metrics, MetaData.loc[cell_IDs, 'Region']])


# # %% create figure: violins

# # # initialise figure
# # fig, axs = plt.subplots(nrows=2,
# #                         ncols=2,
# #                         figsize=get_figure_size(width=160, height=160),
# #                         layout='constrained',
# #                         dpi=600)

# # # flatten axes array
# # axs_metrics_violins = axs.flatten()

# # # create plotting dict for plotting order
# # order_dict = {1: 'enclosing_radius',
# #               2: 'critical_radius',
# #               3: 'max_intersections'}

# # # create list for subplot alphabetical labels
# # alpha_labels = ['', 'B', 'C', 'D']

# # # plot different sholl metrics
# # for axs_idx, metric in order_dict.items():

# #     # set axis
# #     ax = axs[axs_idx]

# #     # set label for subplot
# #     ax.set_title(alpha_labels[axs_idx], fontsize=12, loc='left')

# #     # violin plots
# #     violins = sbn.violinplot(data = sholl_metrics_plot_df,
# #                              x = 'neurite_type',
# #                              y = metric,
# #                              hue='Region',
# #                              hue_order = ['MeA', 'BAOT'],
# #                              bw=0.4,
# #                              inner=None,
# #                              split = True,
# #                              gap = 0.15,
# #                              ax=ax,
# #                              linewidth=1,
# #                              palette=region_colors,
# #                              zorder = 1,
# #                              density_norm= 'width')

# #     # edit lines of quarts
# #     # for v_idx in np.arange(len(order_dict.items()) * len(sholl_metrics_plot_df['Region'].drop_duplicates())):

# #     #     for l_idx in np.arange(0, 3):
# #     #         all_l_idx = v_idx * 3 + l_idx

# #     #         if v_idx % 2 == 0:
# #     #             violins.lines[all_l_idx].set_color(region_colors['MeA'])
# #     #         else:
# #     #             violins.lines[all_l_idx].set_color(region_colors['BAOT'])

# #     # edit color of edges and faces
# #     for n_idx, neurite_type in enumerate(neurite_types):
        
# #         for r_idx, region in enumerate(['MeA', 'BAOT']):
# #             v_idx = n_idx * 2 + r_idx
            
# #             violins.collections[v_idx].set_edgecolor(region_colors[region])
# #             violins.collections[v_idx].set_facecolor(region_colors[region])

# #     # swarmplots
# #     swarm = sbn.swarmplot(data=sholl_metrics_plot_df,
# #                           x='neurite_type',
# #                           y=metric,
# #                           hue='Region',
# #                           hue_order = ['MeA', 'BAOT'],
# #                           ax=ax,
# #                           size=2,
# #                           color=colors_dict['primecolor'],
# #                           dodge=True,
# #                           zorder = 2)
    
# #     # errorbar
# #     for neurite_idx, neurite_type in enumerate(neurite_types):
# #         for region_x, region in zip([-0.1, +0.1], ['MeA', 'BAOT']):
            
            
# #             if neurite_type == 'axons':
# #                 cur_cell_IDs_dict = cell_IDs_w_axon_dict
# #             else:
# #                 cur_cell_IDs_dict = cell_IDs_dict
            
# #             # data set per type and region for mean and std calculation
# #             sholl_metric_per_type_n_region = sholl_metrics[neurite_type].loc[cur_cell_IDs_dict[region], :]
# #             sholl_metric_mean = sholl_metric_per_type_n_region[metric].mean()
# #             sholl_metric_std = sholl_metric_per_type_n_region[metric].std()
# #             sholl_metric_median = sholl_metric_per_type_n_region[metric].median()
            
# #             ax.errorbar(x = neurite_idx + region_x,
# #                         y = sholl_metric_mean,
# #                         yerr = sholl_metric_std,
# #                         fmt='_', 
# #                         markersize = 6,
# #                         markerfacecolor = 'none',
# #                         capsize = 2,
# #                         color= colors_dict['primecolor'],
# #                         linewidth = 1,
# #                         label = '_nolegend_',
# #                         zorder = 3)
            
# #             # plot median
# #             ax.scatter(x = neurite_idx + region_x,
# #                        y = sholl_metric_median,
# #                        marker='D', 
# #                        s = 5,
# #                        color=colors_dict['primecolor'],
# #                        linewidth = 1,
# #                        label = '_nolegend_',
# #                        zorder = 4)
        
    
# #     # get handles and labels of legend
# #     legend_handles, legend_labels = ax.get_legend_handles_labels()

# #     # remove seaborn legend
# #     if axs_idx != 3:
# #         ax.legend().set_visible(False)
# #     else:
# #         ax.legend(handles = legend_handles[:2],
# #                   labels = legend_labels[:2],
# #                   loc='upper right',
# #                   ncol=1)


# #     # x
# #     ax.set_xlim(0-0.5, 2+0.5)
# #     ax.set_xticks(ticks=np.arange(0, 2 + 0.5, 1),
# #                   labels=[label.title() for label in neurite_types], rotation=45)
# #     ax.spines['bottom'].set_bounds([0, 2])
# #     ax.set_xlabel('Type of neurites')

# #     # y
# #     if metric == 'critical_radius':
# #         ax.set_ylabel('Critical radius\nper cell [µm]')
# #         ymin = 0
# #         ymax = 300
# #         ystep = 100
# #         ystepminor = 25

# #     elif metric == 'enclosing_radius':
# #         ax.set_ylabel('Enclosing radius\nper cell [µm]')
# #         ymin = 0
# #         ymax = 500
# #         ystep = 100
# #         ystepminor = 25

# #     elif metric == 'max_intersections':
# #         ax.set_ylabel('Max. number of intersections\nper cell [#]')
# #         ymin = 0
# #         ymax = 50
# #         ystep = 25
# #         ystepminor = 5

# #     # define ypad relative to range
# #     ypad = (ymax - ymin) * 0.03

# #     # apply y axis settings
# #     ax.set_ylim(ymin - ypad, ymax)
# #     ax.set_yticks(ticks=np.arange(ymin, ymax + ystepminor, ystep))
# #     ax.set_yticks(ticks=np.arange(
# #         ymin, ymax + ystepminor, ystepminor), minor=True)
# #     ax.spines['left'].set_bounds([ymin, ymax])

# #     [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# # # delete first subplot
# # fig_metrics_violins.delaxes(axs_metrics_violins[0])

# # # align labels
# # fig_metrics_violins.align_labels()

# # # show plot
# # plt.show()

# # # save plot
# # sholl_metrics_violins_fig_dir = join(cell_morph_plots_dir, 'figure-sholl_metrics-violins')


# # # %% sholl metrics statistics

# # from scipy.stats import normaltest, mannwhitneyu, ks_1samp
# # from scipy import stats
# # from statsmodels.stats.multitest import multipletests # bonferroni correction

# # # create dataframe for statistics measures
# # sholl_metrics_normaltest = pd.DataFrame()
# # sholl_metrics_mannwhitneyu = pd.DataFrame()

# # # iterate through metric
# # for sholl_metric in ['enclosing_radius', 'critical_radius', 'max_intersections']:
    
# #     # iterate through neurite type
# #     for neurite_type in neurite_types:
        
# #         # set cell_IDs per neurite types
# #         # get cell_IDs per neurite_type
# #         if neurite_type != 'axons':
# #             cell_IDs_perNeuriteType = cell_IDs_dict
# #         elif neurite_type == 'axons':
# #             cell_IDs_perNeuriteType = cell_IDs_w_axon_dict
        
# #         # iterate through regions
# #         for region in ['MeA', 'BAOT']:
            
# #             # get cell_IDs per neurite_type and region
# #             cell_IDs_perNeuriteType_n_region = cell_IDs_perNeuriteType[region]

# #             # run normal test
# #             normaltest_res = normaltest(sholl_metrics[neurite_type].loc[cell_IDs_perNeuriteType_n_region, sholl_metric])

# #             # write to dataframe
# #             sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-{region}', 'normaltest_statistic'] = normaltest_res.statistic
# #             sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-{region}', 'normaltest_p_value'] = normaltest_res.pvalue
            
# #             # write boolean
# #             if normaltest_res.pvalue < 0.05 :
# #                 sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-{region}', 'normaltest-normally_distributed'] = True
# #             else:
# #                 sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-{region}', 'normaltest-normally_distributed'] = False
            
# #             # run rest
# #             kstest_res = ks_1samp(sholl_metrics[neurite_type].loc[cell_IDs_perNeuriteType_n_region, sholl_metric],
# #                                   stats.norm.cdf,
# #                                   nan_policy='omit')
            
# #             # write to dataframe
# #             sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-{region}', 'kstest_statistic'] = kstest_res.statistic
# #             sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-{region}', 'kstest_pvalue'] = kstest_res.pvalue
            
# #             # write boolean
# #             if kstest_res.pvalue > 0.05 :  
# #                 sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-{region}', 'kstest-normally_distributed'] = True
# #             else:
# #                 sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-{region}', 'kstest-normally_distributed'] = False

  
# #         # compare both regions
# #         sholl_metric_perNeurite_MeA = sholl_metrics[neurite_type].loc[cell_IDs_perNeuriteType['MeA'], sholl_metric]
# #         sholl_metric_perNeurite_BAOT = sholl_metrics[neurite_type].loc[cell_IDs_perNeuriteType['BAOT'], sholl_metric]

# #         if not (sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-MeA', 'kstest-normally_distributed'] and sholl_metrics_normaltest.at[f'{sholl_metric}-{neurite_type}-BAOT', 'kstest-normally_distributed']):
# #             # run test
# #             mannwhitneyu_res = mannwhitneyu(sholl_metric_perNeurite_MeA, sholl_metric_perNeurite_BAOT, 
# #                                             alternative = 'two-sided', 
# #                                             nan_policy='omit')

# #             # write to dataframe
# #             sholl_metrics_mannwhitneyu.at[f'{sholl_metric}-{neurite_type}', 'mannwhitneyu_statistic'] = mannwhitneyu_res.statistic
# #             sholl_metrics_mannwhitneyu.at[f'{sholl_metric}-{neurite_type}', 'mannwhitneyu_pvalue'] = mannwhitneyu_res.pvalue
            
# #             # write boolean
# #             if mannwhitneyu_res.pvalue < 0.05:
# #                 sholl_metrics_mannwhitneyu.at[f'{sholl_metric}-{neurite_type}', 'mannwhitneyu-statistically_different'] = True
# #             else:
# #                 sholl_metrics_mannwhitneyu.at[f'{sholl_metric}-{neurite_type}', 'mannwhitneyu-statistically_different'] = False
        
        
# #     # Bonferroni correction
# #     # get rows of pvalues
# #     rows = [sholl_metric + '-' + neurite_type for neurite_type in neurite_types]
    
# #     # get pvalues
# #     mannwhitneyu_pvalues = sholl_metrics_mannwhitneyu.loc[rows, 'mannwhitneyu_pvalue']
    
# #     # apply bonferroni correction
# #     rejects, pvals_corrected, _, _ = multipletests(mannwhitneyu_pvalues,
# #                                                    alpha=0.05, 
# #                                                    method='bonferroni')
    
# #     # write to dataframe
# #     sholl_metrics_mannwhitneyu.loc[rows, 'bonferroni_pvalue'] = pvals_corrected
    
# #     # write boolean
# #     sholl_metrics_mannwhitneyu.loc[rows, 'bonferroni-statistical_difference'] = sholl_metrics_mannwhitneyu.loc[rows, 'bonferroni_pvalue'] < 0.05


# # # save statistics dataframes
# # # sholl_metrics_normaltest.to_excel(join(sholl_metrics_violins_fig_dir, 'sholl_metrics_normaltest.xlsx'), index_label='sholl_metric-neurite_type-region')
# # # sholl_metrics_mannwhitneyu.to_excel(join(sholl_metrics_violins_fig_dir, 'sholl_metrics_mannwhitneyu.xlsx'), index_label='sholl_metric-neurite_type')


# # %% figure for critical vs enclosing radius

# # inititalise figure
# fig_sholl_scatter, axs_sholl_scatter = plt.subplots(nrows=3,
#                                                     ncols=2,
#                                                     figsize=get_figure_size(width=160, height=225),
#                                                     layout='constrained',
#                                                     dpi=600)

# # flatten axes array
# axs_sholl_scatter = axs_sholl_scatter.flatten()

# # create list for subplot alphabetical labels
# # alpha_labels = ['A', 'B', 'C', 'D', 'E', 'F']

# # create dict for axis titles 
# axs_titles = {'MeA' : {'neurites' : '$\mathregular{A_{i}}$: MeA neurites',
#                        'dendrites': '$\mathregular{B_{i}}$: MeA dendrites',
#                        'axons'    : '$\mathregular{C_{i}}$: MeA axons'}, 
#               'BAOT': {'neurites' : '$\mathregular{A_{ii}}$: BAOT neurites',
#                        'dendrites': '$\mathregular{B_{ii}}$: BAOT dendrites',
#                        'axons'    : '$\mathregular{C_{ii}}$: BAOT axons'}}

# # initialise color-coding of max intersections
# cmap_str = 'viridis'

# # min max normalize time for color-code
# norm = mtl.colors.Normalize(vmin=0,
#                             vmax=40)

# # create mappable color-map object
# cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# # add colorbar
# fig_sholl_scatter.colorbar(cmap, ax=axs_sholl_scatter[-2::],
#                            orientation='horizontal',
#                            label='Max. number of intersections [#]',
#                            fraction=0.05,
#                            aspect=70)

# # create plotting dict for plotting order
# order_dict = {'MeA' : {'neurites' : 0, 'dendrites': 2, 'axons': 4},
#               'BAOT': {'neurites' : 1, 'dendrites': 3, 'axons': 5}}

# # create dictionary for sorting of inset plots
# inset_sort = {'MeA': True, 'BAOT': False}

# # create scatter plots splited by region and neurite_type
# # loop through region = columns
# for region in ['MeA', 'BAOT']:
    
#     # loop through neurite_types = rows
#     for neurite_type in neurite_types:
        
#         # set axis
#         ax = axs_sholl_scatter[order_dict[region][neurite_type]]
    
#         # subplot title
#         ax.set_title(axs_titles[region][neurite_type], 
#                      fontsize=12, 
#                      loc='left')
        
#         # set cell_IDs
#         if neurite_type == 'dendrites' or neurite_type == 'neurites':
#             plt_cell_IDs = cell_IDs_dict
#         elif neurite_type == 'axons':
#             plt_cell_IDs = cell_IDs_w_axon_dict
            
#         # plot diagonal line in plot
#         ax.plot([0, 400], [0, 400],
#                 color = colors_dict['primecolor'],
#                 linewidth = 1,
#                 linestyle = 'dashed',
#                 alpha = 0.5, 
#                 zorder = 0)
        
#         # plot scatter plots
#         ax.scatter(x=sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['enclosing_radius'],
#                    y=sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['critical_radius'],
#                    color=cmap.to_rgba(sholl_metrics[neurite_type].loc[plt_cell_IDs[region]]['max_intersections']),
#                    s=15,
#                    zorder = 1)
        
        
#         # edit main plot axes
#         # x
#         xmin = 0
#         xmax = 500
#         xstep = 100
#         xstepminor = 25
#         xpad = 10
        
#         ax.set_xlim(xmin-xpad, xmax)
#         ax.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep))
#         ax.set_xticks(ticks=np.arange(
#             xmin, xmax + ystepminor, xstepminor), minor=True)
#         ax.spines['bottom'].set_bounds([xmin, xmax])

#         # y
#         ymin = 0
#         ymax = 400
#         ystep = 100
#         ystepminor = 25
#         ypad = 10
        
#         ax.set_ylim(ymin-ypad, ymax)
#         ax.set_yticks(ticks=np.arange(ymin, ymax + 1, ystep))
#         ax.set_yticks(ticks=np.arange(
#             ymin, ymax + ystepminor, ystepminor), minor=True)
#         ax.spines['left'].set_bounds([ymin, ymax])
        
#         # remove top and right spines
#         [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
        
        
#         ### add insets

#         # ([left, bottom, width, height]), percentages
#         ax_inset = ax.inset_axes([0.1, 0.7, 0.25, 0.25])
        
#         # plot in inset
#         ax_inset.plot([0, 400], [0, 400],
#                       color = colors_dict['primecolor'],
#                       linewidth = 0.5,
#                       linestyle = 'dashed',
#                       alpha = 0.5, 
#                       zorder = 0)
        
        
#         # scatterplot in inset
#         scatter = sbn.scatterplot(data=sholl_metrics_plot_df[sholl_metrics_plot_df['neurite_type'] == neurite_type].sort_values(by=['Region'], ascending=inset_sort[region]),
#                                   x='enclosing_radius',
#                                   y='critical_radius',
#                                   hue='Region',
#                                   palette=region_colors,
#                                   ax=ax_inset,
#                                   s=6,
#                                   linewidth=0, 
#                                   zorder = 1)
        
#         # remove seaborn legend
#         ax_inset.legend().set_visible(False)
        
#         # edit inset axis
#         ax_inset.set_xlim(xmin-xpad, xmax)
#         ax_inset.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep), labels=[])
#         ax_inset.set_xlabel('')
#         ax_inset.tick_params('both', length=1.5, which='major')
        
#         ax_inset.set_ylim(ymin-ypad, ymax)
#         ax_inset.set_yticks(ticks=np.arange(ymin, ymax + 1, ystep), labels=[])
#         ax_inset.set_ylabel('')
        
#         # remove top and right spines
#         [ax_inset.spines[spine].set_visible(False) for spine in ['top', 'right']]

# # axis labels
# [ax.set_ylabel('Critical radius [µm]') for ax in axs_sholl_scatter[::2]] 
# [ax.set_xlabel('Enclosing radius  [µm]') for ax in axs_sholl_scatter]

# # align labels
# fig_sholl_scatter.align_labels()

# # show figure
# plt.show()

# # save figure
# crit_v_enclos_figure_dir = join(cell_morph_plots_dir, 'figure-sholl_metrics-enclosing_v_critical_radius')

# # save_figures(fig_sholl_scatter, 'sholl_metrics-crit_v_enclos-figure', 
# #              save_dir = crit_v_enclos_figure_dir,
# #              darkmode_bool = darkmode_bool, 
# #              figure_format = 'both')

# # %% collected data to save

# # create dataframe
# sholl_metrics_df = pd.DataFrame(index = sholl_cell_IDs)

# # iterate through neurite_type
# for neurite_type in neurite_types:
    
#     # iterate through regions
#     for region in ['MeA', 'BAOT']:
        
#         # get sholl metric per Type
#         sholl_metric_perType = sholl_metrics[neurite_type]

#         # limit to current region
#         sholl_metric_perType_n_region = sholl_metric_perType[sholl_metric_perType['Region'] == region]
        
#         # get current cell_IDs
#         sholl_metric_region_cell_IDs = sholl_metric_perType_n_region.index.to_list()
        
#         # iterate through metrics
#         for sholl_metric in ['enclosing_radius', 'critical_radius', 'max_intersections']:

#             # write to dataframe
#             sholl_metrics_df.loc[sholl_metric_region_cell_IDs, f'{sholl_metric}-{neurite_type}-{region}'] = sholl_metric_perType_n_region[sholl_metric]

# # save dataframe
# # sholl_metrics_df.to_excel(join(crit_v_enclos_figure_dir, 'sholl_metrics_sorted.xlsx'), index_label = 'cell_ID')