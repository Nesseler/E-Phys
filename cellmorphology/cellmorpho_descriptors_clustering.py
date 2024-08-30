# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:22:45 2024

@author: nesseler
"""

from cellmorphology.cellmorph_packages import mtl, plt, sbn, pd, np, join

from parameters.directories_win import table_file, cell_morph_descrip_dir, cell_morph_plots_dir

from cellmorphology.cellmorph_parameters import cell_IDs_toDrop

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# get cell_IDs:
cell_IDs = MetaData[MetaData['reconstructed'] == 1].index.to_list()

# set list of neurite types
neurite_types = ['neurites', 'dendrites', 'axons']

# set directory for clustering figures
cellmorph_clustering_fig_dir = join(cell_morph_plots_dir, 'clustering')


# %%

# create dataframe that contains all parameters
cellmorph_descriptors = pd.DataFrame(index = cell_IDs)
cellmorph_descriptors.index.name = 'cell_ID'


# %% load height & width

# load
height_width = pd.read_excel(join(cell_morph_descrip_dir, 'height_width.xlsx'), index_col = 'cell_ID')

# get list of columns to use
height_width_cols = [col for col in height_width.columns.to_list() if 'width' in col and 'neurites' not in col or 'height' in col and 'neurites' not in col]

# concat
cellmorph_descriptors = pd.concat([cellmorph_descriptors, height_width[height_width_cols]], axis = 1)


# %% load number of points

# load
n_primary = pd.read_excel(join(cell_morph_descrip_dir, 'n_primary_points.xlsx'), index_col = 'cell_ID')
n_terminal = pd.read_excel(join(cell_morph_descrip_dir, 'n_terminal_points.xlsx'), index_col = 'cell_ID')
n_bifurcation_ratios = pd.read_excel(join(cell_morph_descrip_dir, 'bifurcation_ratios.xlsx'), index_col = 'cell_ID')

# concat
cellmorph_descriptors = pd.concat([cellmorph_descriptors, 
                                   n_primary[['dendritic_primaries', 'axonic_primaries']],
                                   n_terminal[['dendritic_terminals', 'axonic_terminals']],
                                   n_bifurcation_ratios[['bifurcation_ratio_dendritic', 'bifurcation_ratio_axonic']]]
                                  , axis = 1)


# %% load total cable lengths

# load
total_cable_length = pd.read_excel(join(cell_morph_descrip_dir, 'total_cable_length.xlsx'), index_col = 'cell_ID')

# filter out neurites
total_cable_length_cols = [col for col in total_cable_length.columns.to_list() if 'neurites' not in col]

# concat
cellmorph_descriptors = pd.concat([cellmorph_descriptors, total_cable_length[total_cable_length_cols]], axis = 1)


# rename columns to unify naming scheme
for neurite_prefix, neurite_type in zip(['dendritic', 'axonic'], neurite_types[1:]):

    cellmorph_descriptors.rename(columns = {f'{neurite_prefix}_primaries'         : f'{neurite_type}-n_primaries', 
                                            f'{neurite_prefix}_terminals'         : f'{neurite_type}-n_terminals',
                                            f'bifurcation_ratio_{neurite_prefix}' : f'{neurite_type}-bifurcation_ratio',
                                            f'total_cable_length-{neurite_type}'  : f'{neurite_type}-total_cable_length'},
                                            inplace = True)


# %% load critical and enclosing radius and max number of intersections

# load
# sholl metrics
sholl_metrics = {'neurites': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_neurites.xlsx'), index_col='cell_ID'),
                 'dendrites': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_dendrites.xlsx'), index_col='cell_ID'),
                 'axons': pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics_axons.xlsx'), index_col='cell_ID')}

for neurite_type in ['dendrites', 'axons']:
    
    # get sholl metrics
    sholl_metrics_pertype = sholl_metrics[neurite_type]

    # rename columns with neurite_type
    sholl_metrics_pertype.rename(columns = {'critical_radius'   : f'{neurite_type}-critical_radius', 
                                            'enclosing_radius'  : f'{neurite_type}-enclosing_radius', 
                                            'max_intersections' : f'{neurite_type}-max_intersections'},
                                 inplace = True)

    # concatenate to descriptors
    cellmorph_descriptors = pd.concat([cellmorph_descriptors, sholl_metrics_pertype], axis = 1)
 

# %% load circular means data

# load
circular_means = pd.read_excel(join(cell_morph_descrip_dir, 'circular_means.xlsx'), index_col = 'cell_ID')

# filter out neurites
circular_means_cols = [col for col in circular_means.columns.to_list() if 'neurites' not in col]

# concat
# cellmorph_descriptors = pd.concat([cellmorph_descriptors, circular_means[circular_means_cols]], axis = 1)
    

# %% load AIS data

# load
axon_data = pd.read_excel(join(cell_morph_descrip_dir, 'axon_data.xlsx'), index_col='cell_ID')

# define columns to add
axon_cols = ['distance to soma']
axon_cols_ext = ['length(µm)', 'distance to soma', 'length + distance', 'source']

# rename columns with neurite_type
axon_data.rename(columns = {'distance to soma'   : 'AIS distance to soma'},
                inplace = True)

# concat
cellmorph_descriptors = pd.concat([cellmorph_descriptors, axon_data['AIS distance to soma']], axis = 1)


# %% load spine categorization

# load
spines_df = pd.read_excel(join(cell_morph_descrip_dir, 'spines.xlsx'), index_col='cell_ID')

# %%

# remove cells that do not contain all analysed values
cellmorph_descriptors.drop(index = cell_IDs_toDrop, inplace = True)

# replace nan values in cell descriptor table
# cellmorph_descriptors.fillna(value = 0, inplace = True)

# drop cells without axons
cellmorph_descriptors = cellmorph_descriptors.dropna(axis=0)


# %% sort column names after their name
sorted_columns = []

for neurite_type in ['dendrites', 'axons']: ### !!!!axons
    
    for col_parameter in ['width', 'height', 'n_primaries', 'n_terminals', 'bifurcation_ratio', 'total_cable_length', 'critical_radius', 'enclosing_radius', 'max_intersections']: #, 'circmean']:
    
        sorted_columns.append(f'{neurite_type}-{col_parameter}')

# append last column
sorted_columns.append('AIS distance to soma')

# apply sorting
cellmorph_descriptors = cellmorph_descriptors[sorted_columns]

# test: drop axons-n_primaries
cellmorph_descriptors.drop(columns = ['axons-n_primaries', 'axons-n_terminals'], inplace = True) #, 'axons-width', 'axons-height']

# %% normalise cell descriptors

# min-max normalize cellmorph matrix
cellmorph_descriptors_minmax = (cellmorph_descriptors - cellmorph_descriptors.min()) / (cellmorph_descriptors.max() - cellmorph_descriptors.min())

# z-score cellmorph matrix
cellmorph_descriptors_zscored = (cellmorph_descriptors - cellmorph_descriptors.mean()) / cellmorph_descriptors.std()


# from scipy.stats import circmean, circstd

# # exception in zscoring for circular mean
# for neurite_type in ['dendrites']:  ### !!!!axons
    
#     # get list of circular means
#     circ_means_pertype = cellmorph_descriptors[f'{neurite_type}-circmean']
    
#     # calc circular mean and std of means
#     population_circ_mean = circmean(circ_means_pertype)
#     population_circ_std = circstd(circ_means_pertype)
    
#     cellmorph_descriptors_zscored[f'{neurite_type}-circmean'] = [(cell_circ_mean - population_circ_mean) / population_circ_std for cell_circ_mean in circ_means_pertype]



# %% initialize plotting

from functions.functions_plotting import save_figures, get_colors, get_figure_size
from functions.functions_useful import round_up_to_base, round_down_to_base
from cellmorphology.cellmorph_parameters import cell_coordinates_field_of_view


# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% distributions of all parameters

# distributions_data = cellmorph_descriptors_zscored

def create_data_distribution_figure(distributions_data = cellmorph_descriptors_zscored):
    
    fig_dist, ax_dist = plt.subplots(nrows = 1,
                                     ncols = 1,
                                     layout = 'constrained',
                                     figsize = get_figure_size(width = 150, height = 125),
                                     dpi = 600)
    
    # melt dataframe to plot
    distributions_data_melted = distributions_data.melt(var_name = 'parameter')
    
    
    # create list of parameters
    parameters = distributions_data.columns.to_list()
    
    # create list of colors for parameters
    p_cmap = plt.get_cmap('viridis', 20)
    
    # plot violins
        
    violins = sbn.violinplot(data = distributions_data_melted,
                             x = 'parameter',
                             y = 'value',
                             ax = ax_dist,
                             linewidth = 1,
                             inner = 'quart',
                             scale='width',
                             hue = True, hue_order=[True, False], split = True)
    
    
    for p_idx, param in enumerate(parameters):
        # set line color of quarts
        for l_idx in range(3):
            violins.lines[p_idx * 3 + l_idx].set_color(p_cmap(p_idx))
        
        # set edge color of violin
        violins.collections[p_idx].set_edgecolor(p_cmap(p_idx))
    
        # set facecolor of violin
        violins.collections[p_idx].set_facecolor('None')
    
    # plot swarm
    
    swarms = sbn.swarmplot(data = distributions_data_melted,
                           x = 'parameter',
                           y = 'value', 
                           ax = ax_dist,
                           s = 1,
                           color=colors_dict['primecolor'])
    
    # plot error bar
    for p_idx, param in enumerate(parameters):
        ax_dist.errorbar(x = p_idx+0.3,
                         y = distributions_data[param].mean(),
                         yerr = distributions_data[param].std(),
                         fmt='_', 
                         markersize = 4,
                         markerfacecolor = 'none',
                         capsize = 1,
                         color=p_cmap(p_idx),
                         linewidth = 1,
                         label = '_nolegend_')
        
    # edit seaborn legend
    ax_dist.legend().set_visible(False)
    
    # edit axis
    # x
    xmin = 0
    xmax = len(parameters) - 1
    xpad = 0.6
    
    ax_dist.set_xlim([xmin - xpad, xmax + xpad])
    ax_dist.set_xticklabels(labels = parameters,
                            rotation = 90)
    ax_dist.spines['bottom'].set_bounds([xmin, xmax])
    
    ax_dist.set_xlabel('')
    
    # y
    ymin = round_down_to_base(distributions_data_melted.min().value, 2)
    ymax = round_up_to_base(distributions_data_melted.max().value, 2)
    ypad = 0.5
    ystep = 2
    ystepminor = 0.5
    
    ax_dist.set_ylim([ymin - ypad, ymax + ypad])
    ax_dist.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystep))
    ax_dist.set_yticks(ticks = np.arange(ymin, ymax, ystepminor), minor = True)
    ax_dist.spines['left'].set_bounds([ymin, ymax])
    
    # remove spines
    [ax_dist.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    return fig_dist, ax_dist
    
    
# # %% plot z-scored distribution

# fig_dist_zscored, ax_dist_zscored = create_data_distribution_figure(cellmorph_descriptors_zscored)    

# # set axis label
# ax_dist_zscored.set_ylabel('Z-scored parameter value [std]')

# # show plot
# plt.show()


# # %% plot min max distribution

# fig_dist_minmax, ax_dist_minmax = create_data_distribution_figure(cellmorph_descriptors_minmax)    

# # set axis label
# ax_dist_minmax.set_ylabel('Min-Max-Normalised parameter value')

# # edit axis
# ymin = -0.5
# ymax = 1.5
# ypad = 0.15
# ystep = 1
# ystepminor = 0.25
    
# ax_dist_minmax.set_ylim([ymin - ypad, ymax + ypad])
# ax_dist_minmax.set_yticks(ticks = np.arange(0, ymax+ ystepminor, ystep))
# ax_dist_minmax.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystepminor), minor = True)
# ax_dist_minmax.spines['left'].set_bounds([ymin, ymax])

# # show plot
# plt.show()
    

# # %% plot min max distribution

# fig_dist, ax_dist = create_data_distribution_figure(cellmorph_descriptors)    

# # set axis label
# ax_dist.set_ylabel('Parameter value [µm / # / ]')

# # edit axis
# ymin = 0
# ymax = 6000
# ypad = 100
# ystep = 2000
# ystepminor = 500
    
# ax_dist.set_ylim([ymin - ypad, ymax + ypad])
# ax_dist.set_yticks(ticks = np.arange(0, ymax+ ystepminor, ystep))
# ax_dist.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystepminor), minor = True)
# ax_dist.spines['left'].set_bounds([ymin, ymax])

# # show plot
# plt.show()    


# # %% create figure of all plots

# rewrite function
def plot_data_distribution(ax, distributions_data = cellmorph_descriptors_zscored):
        
    # melt dataframe to plot
    distributions_data_melted = distributions_data.melt(var_name = 'parameter')
    
    
    # create list of parameters
    parameters = distributions_data.columns.to_list()
    
    # create list of colors for parameters
    p_cmap = plt.get_cmap('viridis', 20)
    
    # plot violins
        
    violins = sbn.violinplot(data = distributions_data_melted,
                              x = 'parameter',
                              y = 'value',
                              ax = ax,
                              linewidth = 1,
                              inner = 'quart',
                              scale='width',
                              hue = True, hue_order=[True, False], split = True)
    
    
    for p_idx, param in enumerate(parameters):
        # set line color of quarts
        for l_idx in range(3):
            violins.lines[p_idx * 3 + l_idx].set_color(p_cmap(p_idx))
        
        # set edge color of violin
        violins.collections[p_idx].set_edgecolor(p_cmap(p_idx))
    
        # set facecolor of violin
        violins.collections[p_idx].set_facecolor('None')
    
    # plot swarm
    
    swarms = sbn.swarmplot(data = distributions_data_melted,
                            x = 'parameter',
                            y = 'value', 
                            ax = ax,
                            s = 1,
                            color=colors_dict['primecolor'])
    
    # plot error bar
    for p_idx, param in enumerate(parameters):
        ax.errorbar(x = p_idx+0.3,
                    y = distributions_data[param].mean(),
                    yerr = distributions_data[param].std(),
                    fmt='_', 
                    markersize = 4,
                    markerfacecolor = 'none',
                    capsize = 1,
                    color=p_cmap(p_idx),
                    linewidth = 1,
                    label = '_nolegend_')
        
    # edit seaborn legend
    ax.legend().set_visible(False)
    


# initialize figure
fig_dist_all, axs_dist_all = plt.subplots(nrows = 3,
                                          ncols = 1,
                                          layout = 'constrained',
                                          figsize = get_figure_size(width = 150, height = 250),
                                          dpi = 600)

for data_idx, data_df in enumerate([cellmorph_descriptors, cellmorph_descriptors_minmax, cellmorph_descriptors_zscored]):
    
    # set axis
    ax = axs_dist_all[data_idx]
    
    # plot
    plot_data_distribution(ax, distributions_data = data_df)
    
    # edit axis
    # x
    xmin = 0
    xmax = len(cellmorph_descriptors.columns) - 1
    xpad = 0.6
    

    # y
    if data_idx == 0:
        ymin = 0
        ymax = 6000
        ypad = 100
        ystep = 2000
        ystepminor = 500
        xticklabels = []
        ylabel = 'Parameter value [µm / # / ]'
        axistitle = 'Original parameters'
        
    elif data_idx == 1:
        ymin = -0.5
        ymax = 1.5
        ypad = 0.15
        ystep = 1
        ystepminor = 0.25
        xticklabels = []
        ylabel = 'Parameter value'
        axistitle = 'Min-Max-Normalized parameters'
            
    elif data_idx == 2:
        ymin = -4
        ymax = 6
        ypad = 0.25
        ystep = 2
        ystepminor = 0.25
        xticklabels = cellmorph_descriptors.columns
        ylabel = 'Parameter value [std]'
        axistitle = 'Z-scored parameters'
        
        
    # set axis title
    ax.set_title(axistitle,
                  fontsize=12, 
                  loc='left',
                  x = 0)
        
    # apply axis changes 
    ax.set_xlim([xmin - xpad, xmax + xpad])
    ax.set_xticklabels(labels = xticklabels,
                        rotation = 90)
    ax.spines['bottom'].set_bounds([xmin, xmax])
    
    ax.set_xlabel('')

    ax.set_ylim([ymin - ypad, ymax + ypad])
    ax.set_yticks(ticks = np.arange(round_up_to_base(ymin, 1), ymax+ ystepminor, ystep))
    ax.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystepminor), minor = True)
    ax.spines['left'].set_bounds([ymin, ymax])
    ax.set_ylabel(ylabel)
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
   
# align labels
fig_dist_all.align_labels()
    
# show figure
plt.show()

# # save figure
# save_figures(fig_dist_all, 'hierarchical_clustering-parameter_distributions-all', 
#               save_dir = cellmorph_clustering_fig_dir,
#               darkmode_bool= darkmode_bool,
#               figure_format= 'png')


# %% correlation analysis


# # initialize figure
# fig_corr_heat, axs_corr_heat = plt.subplots(nrows = 2,
#                                             ncols = 1,
#                                             layout = 'constrained',
#                                             figsize = get_figure_size(width = 160, height = 225),
#                                             dpi = 600,
#                                             sharey = True,
#                                             sharex = True)

# # set string for color map
# cmap_str = 'seismic'

# # get correlation matrix of dataframe
# cellmorph_descriptors_corr = cellmorph_descriptors.corr()

# # plot heatmap of correlation values
# sbn.heatmap(data= cellmorph_descriptors_corr,
#             cmap = cmap_str,
#             vmin = -1,
#             vmax = 1,
#             annot = False,
#             cbar_kws={'label': 'Correlation coefficient', 'ticks' : [-1, 0, 1]},
#             ax = axs_corr_heat[0])

# # set axis title
# axs_corr_heat[0].set_title('A: Pearson correlation coefficient',
#                            fontsize=12, 
#                            loc='left',
#                            x = -0.58)


# # plot heatmap with values above threshold only
# corr_threshold = 0.8

# # replace values in dataframe below threshold with nan
# cellmorph_descriptors_corr_thresh = cellmorph_descriptors_corr[(cellmorph_descriptors_corr > corr_threshold) | (cellmorph_descriptors_corr < -corr_threshold)]


# # define boundaries for color map with threshold
# color_boundaries = [-1, -corr_threshold, corr_threshold, 1]
# n_color = len(color_boundaries) -1

# norm = mtl.colors.BoundaryNorm(boundaries = [-1, -corr_threshold, corr_threshold, 1], ncolors=n_color)

# cmap = plt.get_cmap(cmap_str, n_color)



# # plot heatmap of correlation values
# sbn.heatmap(data= cellmorph_descriptors_corr,
#             cmap = cmap,
#             norm = norm,
#             vmin = -1,
#             vmax = 1,
#             annot = False,
#             cbar_kws={'label': 'Correlation coefficient', 'ticks' : [-1, -0.8, 0.8, 1]},
#             ax = axs_corr_heat[1],
#             xticklabels = 1)

# axs_corr_heat[1].set_title(r'B: Pearson correlation coefficient $\pm$' + str(corr_threshold),
#                            fontsize=12, 
#                            loc='left',
#                            x = -0.58)

# # show figure
# plt.show()

# # save figure
# save_figures(fig_corr_heat, 'hierarchical_clustering-correlation_matrices', 
#              save_dir = cellmorph_clustering_fig_dir,
#              darkmode_bool= darkmode_bool,
#              figure_format= 'png')


# %%

# sbn.clustermap(data = cellmorph_descriptors_zscored,
#                 method = 'ward',
#                 col_cluster=False,
#                 cbar_pos=(0.02, 0.2, 0.05, 0.18))

# plt.show()





# %% heatmap pre clustering

# # initialise figure
# fig_heat, ax_heat = plt.subplots(nrows = 1,
#                                   ncols = 1,
#                                   layout = 'constrained',
#                                   figsize = get_figure_size(width = 150, height = 150),
#                                   dpi = 600)

# #set figure title
# fig_heat.suptitle('Cell morphology parameters pre-clustering', fontsize = 12)

# # set min and max of heatmap
# heatmin = -2
# heatmax = 2

# # plot heatmap
# sbn.heatmap(cellmorph_descriptors_zscored,
#             vmin = heatmin,
#             vmax = heatmax,
#             square = False, 
#             ax = ax_heat, 
#             cmap="flare_r", 
#             yticklabels=False,
#             linewidth = 0,
#             cbar_kws={'label': 'Z-scored parameter value [std]', 'ticks' : np.arange(heatmin, heatmax+1, 1)}) 

    
# # show plot
# plt.show()


# %% hierarchical clustering

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import fcluster

# set scaling method
scaling = 'zscored'

# set dict for clustering parameters
clustering_data_dict = {'minmax' : {'df' : cellmorph_descriptors_minmax, 'heatmin' : 0, 'heatmax'  : 1, 'cbar_label' : 'Min-max normalized', 'c_threshold' : 3},
                        'zscored': {'df' : cellmorph_descriptors_zscored, 'heatmin' : -2, 'heatmax'  : 2, 'cbar_label' : 'Z-scored parameter value [std]', 'c_threshold' : 15}}

# set variables
df_tocluster, heatmin, heatmax, cbar_label, c_threshold = clustering_data_dict[scaling].values()

# calc distance between clusters
ward_clustering_linkage = linkage(df_tocluster, method="ward", metric="euclidean")

# get last few clusters that have been merged by the linkage functions
last_clusters = ward_clustering_linkage[-20:, 2]

# reverse list 
last_clusters_rev = last_clusters[::-1]

# set a list of indices
last_clusters_idc = np.arange(1, len(last_clusters) +1)


# calculate the acceleration of lost distance between clusters
# calculated as the 2nd derivative
acceleration = np.diff(last_clusters, 2)   

# reverse list
acceleration_rev = acceleration[::-1]


### elbow plot ###

fig_elbow, ax_elbow = plt.subplots(nrows = 1,
                                    ncols = 1,
                                    layout = 'constrained',
                                    figsize = get_figure_size(width = 100, height = 100),
                                    dpi = 600)

# set title
fig_elbow.suptitle(f'hierarchical clustering\n{scaling} scaling - elbow plot',
                    fontsize = 12)

# plot cluster distances
ax_elbow.plot(last_clusters_idc, last_clusters_rev,
              lw = 1,
              c = colors_dict['primecolor'],
              label = 'Cluster distance')

# plot acceleration
ax_elbow.plot(last_clusters_idc[:-2] + 1, acceleration_rev,
              lw = 1,
              c = colors_dict['color3'],
              label = 'Acceleration of cluster\ndistance growth')

# set legend
ax_elbow.legend(prop={'size': 9})


# edit axis
# x
ax_elbow.set_xlabel('Number of clusters [#]')

xmin = 1
xmax = last_clusters_idc[-1]
xpad = 0.5

ax_elbow.set_xlim([xmin - xpad, xmax + xpad])
ax_elbow.set_xticks(ticks = np.arange(2, len(last_clusters)+1, 2))
ax_elbow.set_xticks(ticks = last_clusters_idc, minor = True)
ax_elbow.spines['bottom'].set_bounds([xmin, xmax])

# remove spines
[ax_elbow.spines[spine].set_visible(False) for spine in ['top', 'right']]

# show plot
plt.show()

# save figure
save_figures(fig_elbow, f'hierarchical_clustering-elbow_plot-{scaling}', 
              save_dir = cellmorph_clustering_fig_dir,
              darkmode_bool= darkmode_bool,
              figure_format= 'png')




# %% create auxillary dataframe for additional parameters

# get cell_IDs from cellmorph_descriptors
cell_IDs_clustering = cellmorph_descriptors_zscored.index.to_list()

# create dataframe for auxillary heatmap
aux_cellmorph_descriptors = pd.DataFrame(index = cell_IDs_clustering)

# get auxillary data
aux_cellmorph_descriptors['Region'] = MetaData.loc[cell_IDs_clustering, 'Region']
aux_cellmorph_descriptors['spinyness'] = spines_df.loc[cell_IDs_clustering, 'Unnamed: 2']
aux_cellmorph_descriptors['dendrites-circ_mean'] = circular_means['dendrites-circmean']
aux_cellmorph_descriptors['axons-circ_mean'] = circular_means['axons-circmean']
aux_cellmorph_descriptors['AIS root'] = axon_data['source']

# transform region to values
aux_cellmorph_descriptors['Region'] = aux_cellmorph_descriptors['Region'].map({'MeA' : 0, 'BAOT/MeA' : 0.5, 'BAOT' : 1})
aux_cellmorph_descriptors['AIS root'] = aux_cellmorph_descriptors['AIS root'].map({'dendritic' : 0, 'somatic' : 1})

# define dict for cmaps in heatmaps
from cellmorphology.cellmorph_colors import spines_color_dict, neurite_color_dict

heatmap_dict = {'Region'              : {'cmap' : [region_colors[region] for region in ['MeA', 'BAOT/MeA', 'BAOT']]},
                'spinyness'           : {'cmap' : [spines_color_dict['both'][spinyness] for spinyness in ['low', 'moderate', 'high']]},
                'dendrites-circ_mean' : {'cmap' : plt.get_cmap('twilight'), 'vmin' : 0, 'vmax' : np.pi * 2},
                'axons-circ_mean'     : {'cmap' : plt.get_cmap('twilight'), 'vmin' : 0, 'vmax' : np.pi * 2},
                'AIS root'            : {'cmap' : [neurite_color_dict[region]['axons'] for region in ['MeA', 'BAOT']]}}

cbar_dict = {'Region'              : {'ticks' : [0.166, 0.5, 0.833],                  'labels' : ['MeA', 'BAOT/MeA', 'BAOT'], 'range' : [0, 1]},
             'spinyness'           : {'ticks' : [0.333, 1.0, 1.666],                  'labels' : ['low', 'moderate', 'high'], 'range' : [0, 2]},
             'dendrites-circ_mean' : {'ticks' : np.arange(0, np.pi *2 +0.1, np.pi/2), 'labels' : ['p', 'd', 'a', 'v', 'p'],   'range' : [0, np.pi * 2]},
             'axons-circ_mean'     : {'ticks' : np.arange(0, np.pi *2 +0.1, np.pi/2), 'labels' : ['p', 'd', 'a', 'v', 'p'],   'range' : [0, np.pi * 2]},
             'AIS root'            : {'ticks' : [0.25, 0.75],                         'labels' : ['dendritic', 'somatic'],    'range' : [0, 1]}} 

axs_dict = {'Region'              : {'heatmap' :3, 'cbar' : 10},
            'spinyness'           : {'heatmap' :4, 'cbar' : 11}, 
            'dendrites-circ_mean' : {'heatmap' :5, 'cbar' : np.nan},
            'axons-circ_mean'     : {'heatmap' :6, 'cbar' : 12},
            'AIS root'            : {'heatmap' :7, 'cbar' : 13}}


# heatmap test

from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.colors import LinearSegmentedColormap

# calc distance between clusters
ward_clustering_linkage = linkage(cellmorph_descriptors_zscored, method="ward", metric="euclidean")

# get last few clusters that have been merged by the linkage functions
last_clusters = ward_clustering_linkage[-20:, 2]

# reverse list 
last_clusters_rev = last_clusters[::-1]

# set number of clusters
n_clusters = 7

# as halfway point between n_clusters-1 and n_clusters
c_threshold = last_clusters_rev[n_clusters-1] + (last_clusters_rev[n_clusters-2] - last_clusters_rev[n_clusters-1]) / 2


# calculate width ratio
n_cols_heatmap = cellmorph_descriptors_zscored.shape[1]
n_additional_cols = 10
heatmap_width_r = n_cols_heatmap / (n_cols_heatmap + n_additional_cols)
single_col_width_r = 1 / (n_cols_heatmap + n_additional_cols)




# initialise figure
fig_clustering, axs_clustering = plt.subplots(nrows = 1,
                                              ncols = 14,
                                              layout = 'constrained',
                                              figsize = get_figure_size(width = 160, height = 205),
                                              width_ratios=[0.2, heatmap_width_r, 0.015] + [single_col_width_r] * (5) + [0.025] + [single_col_width_r/2] + [single_col_width_r] * (4),
                                              dpi = 600)

# adjust layout of constrained setup
fig_clustering.set_constrained_layout_pads(h_pad=4./72., hspace=0./72.,
                                           w_pad=0.01, wspace=0)

### dendrogram ###
# set axis for dendrogram
ax_dendro = axs_clustering[0]

# plot dendrogram
dendrogram_plot = dendrogram(Z = ward_clustering_linkage, 
                             labels = cellmorph_descriptors_zscored.index.to_list(), 
                             ax = ax_dendro,
                             orientation = 'left',
                             color_threshold = c_threshold)

# plot cluster threshold
ax_dendro.axvline(x = c_threshold, 
                  color = 'k',
                  lw = 1,
                  ls = 'dashed')

plt.rcParams['lines.linewidth'] = 1 

# edit axis
ax_dendro.set_xlabel('Distance')
ax_dendro.set_xticks(ticks = np.arange(0, 25, 10))
ax_dendro.set_xticks(ticks = np.arange(0, 25+1, 5), minor = True)

# remove spines
[ax_dendro.spines[spine].set_visible(False) for spine in ['top', 'right', 'left']]

# get cell IDs of leaves
leave_cell_IDs = dendrogram_plot['ivl']

# invert order of leave cell IDs
leave_cell_IDs = list(reversed(leave_cell_IDs))

# resort dataframe indices
cellmorph_descriptors_zscored_clustered = cellmorph_descriptors_zscored.reindex(leave_cell_IDs)
aux_cellmorph_descriptors_clustered = aux_cellmorph_descriptors.reindex(leave_cell_IDs)


### heatmap ###
# plot heatmap
sbn.heatmap(cellmorph_descriptors_zscored_clustered,
            vmin = heatmin,
            vmax = heatmax,
            square = False,
            xticklabels= 1,
            ax = axs_clustering[1], 
            cmap="coolwarm", 
            yticklabels=False,
            linewidth = 0,
            cbar = False) 

# remove ylabel
axs_clustering[1].set_ylabel('')


### heatmap colorbar ###
ax_cbar = axs_clustering[9]
cmin = -2
cmax = 2
crange = (cmax - cmin) * 0.3

# create normalize object
norm = mtl.colors.Normalize(-2, 2)

# create mappable cmap object
cmap = mtl.cm.ScalarMappable(norm=norm, cmap='coolwarm')
   
# plot colorbar in axis
cbar = fig_clustering.colorbar(mappable = cmap, 
                               cax = ax_cbar, 
                               label = '', 
                               orientation='vertical',
                               drawedges = False)

# calc new limits
cmin_new = cmin - crange
cmax_new = cmax + crange

# apply changes
# axis
ax_cbar.set_ylim([cmin_new, cmax_new])
ax_cbar.spines['left'].set_bounds([cmin_new, cmax_new])

# colorbar
cbar.set_ticks(np.arange(cmin, cmax+0.1, 1))
cbar.outline.set_visible(False)
cbar.set_ticklabels(np.arange(cmin, cmax+0.1, 1, dtype = int))

ax_cbar.set_ylabel('Z-scored parameters [std]')
ax_cbar.yaxis.set_label_position('left')


# get list of auxillary parameters
aux_parameters = aux_cellmorph_descriptors.columns.to_list()

for a_idx, aux_parameter in enumerate(aux_parameters):
    
    # set heat map axis
    ax_heatmap = axs_clustering[axs_dict[aux_parameter]['heatmap']]
    
    # plot heatmap
    heatmap = sbn.heatmap(data = aux_cellmorph_descriptors_clustered.loc[:, [aux_parameter]],
                          ax = ax_heatmap,
                          square= False,
                          xticklabels=1,
                          yticklabels=False,
                          linewidth = 0,
                          cbar = False,
                          **heatmap_dict[aux_parameter])
    
    # rotate xticklabels
    [label.set_rotation(90) for label in ax_heatmap.get_xticklabels()] 
    
    if type(axs_dict[aux_parameter]['cbar']) == int:
        # set colorbar axis
        ax_cbar = axs_clustering[axs_dict[aux_parameter]['cbar']]
        
        # create color map from list of colors
        if 'circ_mean' not in aux_parameter: 
            n_colors = len(heatmap_dict[aux_parameter]['cmap'])
            cmap = LinearSegmentedColormap.from_list(aux_parameter, heatmap_dict[aux_parameter]['cmap'], N=n_colors)
        
            # # min max normalize time for color-code
            norm_min = aux_cellmorph_descriptors.loc[:, [aux_parameter]].min()
            norm_max = aux_cellmorph_descriptors.loc[:, [aux_parameter]].max()
     
        else:
            cmap = heatmap_dict[aux_parameter]['cmap']
            
            # min max normalize time for color-code
            norm_min = 0
            norm_max = np.pi * 2
            
        # create normalize object
        norm = mtl.colors.Normalize(norm_min, norm_max)
        
        # create mappable cmap object
        cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap)
           
        # plot colorbar in axis
        cbar = fig_clustering.colorbar(mappable = cmap, 
                                       cax = ax_cbar, 
                                       label = '', 
                                       orientation='vertical',
                                       drawedges = False)
        
        
        
        # colorbar annotation
        
        # get min and max of y axis to shrink subplot
        cmin = cbar_dict[aux_parameter]['range'][0]
        cmax = cbar_dict[aux_parameter]['range'][-1]
        crange = (cmax - cmin) * 0.75
        
        # calc new limits
        cmin_new = cmin - crange
        cmax_new = cmax + crange
        
        # apply changes
        ax_cbar.set_ylim([cmin_new, cmax_new])
        
        for ctick, clabel in zip(cbar_dict[aux_parameter]['ticks'], cbar_dict[aux_parameter]['labels']):
            
            ax_cbar.text(x = 0.55,
                         y = ctick,
                         s = clabel,
                         fontsize = 7,
                         rotation = 90,
                         va = 'center',
                         ha = 'center')
        
        cbar.set_ticks([])
        cbar.outline.set_visible(False)

# empty axis
for ax_empty in [axs_clustering[2], axs_clustering[8]]:
    [ax_empty.spines[spine].set_visible(False) for spine in ['top', 'right', 'left', 'bottom']]
    ax_empty.set_yticks([])
    ax_empty.set_xticks([])

# set title for subplots
ax_titles = ['A', 'B', 'C', 'D']

for ax_idx, t_idx in zip([0, 1, 3, 9], [0, 1, 2, 3]):
    axs_clustering[ax_idx].set_title(ax_titles[t_idx],
                                     fontsize = 12,
                                     loc = 'left')

# show plot
plt.show()

# save figure
clustering_fig_dir = join(cell_morph_plots_dir, 'figure-hierarchical_clustering')
save_figures(fig_clustering, 'figure-hierarchical_clustering-dendrogram_heatmap_aux-z_scoring-cells_wAxon',
             save_dir = clustering_fig_dir,
             darkmode_bool= darkmode_bool,
             figure_format= 'both')


# save dataset
cellmorph_descriptors_zscored_clustered.to_excel(join(clustering_fig_dir, 'cellmorph_descriptors-cells_wAxons-z_scored-clustered.xlsx'), index_label = 'cell_ID')
aux_cellmorph_descriptors_clustered.to_excel(join(clustering_fig_dir, 'cellmorph_auxillary_descriptors-cells_wAxons-z_scored-clustered.xlsx'), index_label = 'cell_ID')

