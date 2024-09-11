# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 16:56:04 2024

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
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 12})




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


# ### elbow plot ###

# fig_elbow, ax_elbow = plt.subplots(nrows = 1,
#                                     ncols = 1,
#                                     layout = 'constrained',
#                                     figsize = get_figure_size(width = 100, height = 100),
#                                     dpi = 600)

# # set title
# fig_elbow.suptitle(f'hierarchical clustering\n{scaling} scaling - elbow plot',
#                     fontsize = 12)

# # plot cluster distances
# ax_elbow.plot(last_clusters_idc, last_clusters_rev,
#               lw = 1,
#               c = colors_dict['primecolor'],
#               label = 'Cluster distance')

# # plot acceleration
# ax_elbow.plot(last_clusters_idc[:-2] + 1, acceleration_rev,
#               lw = 1,
#               c = colors_dict['color3'],
#               label = 'Acceleration of cluster\ndistance growth')

# # set legend
# ax_elbow.legend(prop={'size': 9})


# # edit axis
# # x
# ax_elbow.set_xlabel('Number of clusters [#]')

# xmin = 1
# xmax = last_clusters_idc[-1]
# xpad = 0.5

# ax_elbow.set_xlim([xmin - xpad, xmax + xpad])
# ax_elbow.set_xticks(ticks = np.arange(2, len(last_clusters)+1, 2))
# ax_elbow.set_xticks(ticks = last_clusters_idc, minor = True)
# ax_elbow.spines['bottom'].set_bounds([xmin, xmax])

# # remove spines
# [ax_elbow.spines[spine].set_visible(False) for spine in ['top', 'right']]

# # show plot
# plt.show()

# # save figure
# # save_figures(fig_elbow, f'hierarchical_clustering-elbow_plot-{scaling}', 
# #               save_dir = cellmorph_clustering_fig_dir,
# #               darkmode_bool= darkmode_bool,
# #               figure_format= 'png')




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

cbar_dict = {'Region'              : {'ticks' : [0.166, 0.5, 0.833],                  'labels' : ['MeA', '', 'BAOT'], 'range' : [0, 1]},
             'spinyness'           : {'ticks' : [0.333, 1.0, 1.666],                  'labels' : ['low', 'moderate', 'high'], 'range' : [0, 2]},
             'dendrites-circ_mean' : {'ticks' : np.arange(0, np.pi *2 +0.1, np.pi/2), 'labels' : ['p', 'd', 'a', 'v', 'p'],   'range' : [0, np.pi * 2]},
             'axons-circ_mean'     : {'ticks' : np.arange(0, np.pi *2 +0.1, np.pi/2), 'labels' : ['p', 'd', 'a', 'v', 'p'],   'range' : [0, np.pi * 2]},
             'AIS root'            : {'ticks' : [0.25, 0.75],                         'labels' : ['dendritic', 'somatic'],    'range' : [0, 1]}} 

axs_dict = {'Region'              : {'heatmap' :4,  'cbar' : 5},
            'spinyness'           : {'heatmap' :6,  'cbar' : 7}, 
            'dendrites-circ_mean' : {'heatmap' :8,  'cbar' : 9},
            'axons-circ_mean'     : {'heatmap' :10, 'cbar' : 11},
            'AIS root'            : {'heatmap' :12, 'cbar' : 13}}


# heatmap test

from scipy.cluster import hierarchy
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
n_additional_cols = 5
heatmap_width_r = n_cols_heatmap / (n_cols_heatmap + n_additional_cols)
single_col_width_r = 1 / (n_cols_heatmap + n_additional_cols)



# initialise figure
fig_clustering, axs_clustering = plt.subplots(nrows = 7,
                                              ncols = 2,
                                              layout = 'constrained',
                                              figsize = get_figure_size(height = 136.983, width = 277.25),
                                              width_ratios = [0.8, 0.2],
                                              height_ratios = [0.15, single_col_width_r*n_cols_heatmap, single_col_width_r, single_col_width_r, single_col_width_r, single_col_width_r, single_col_width_r],
                                              dpi = 600)

# adjust layout of constrained setup
fig_clustering.set_constrained_layout_pads(w_pad=0.01, hspace=0.02, h_pad = 0.02)

# flatten axes array
axs_clustering = axs_clustering.flatten()

# ### dendrogram ###
# set axis for dendrogram
ax_dendro = axs_clustering[0]

# set color palette for dendrogram
#                                   0.1        0.2        0.3        0.4        0.5        0.6        0.7
# hierarchy.set_link_color_palette(['#E6E6E6', '#CCCCCC', '#B3B3B3', '#999999', '#808080', '#666666', '#4D4D4D'])

hierarchy.set_link_color_palette(['#E6E6E6', '#CCCCCC', '#B3B3B3', '#999999', '#808080', '#666666'])

# plot dendrogram
dendrogram_plot = dendrogram(Z = ward_clustering_linkage, 
                             labels = cellmorph_descriptors_zscored.index.to_list(), 
                             ax = ax_dendro,
                             orientation = 'top',
                             color_threshold = c_threshold,
                             above_threshold_color='w')

# plot cluster threshold
ax_dendro.axhline(y = c_threshold, 
                  color = 'w',
                  lw = 0.75,
                  ls = 'dashed')

plt.rcParams['lines.linewidth'] = 0.5 

# edit axis
ax_dendro.set_ylabel('Distance')
ax_dendro.set_yticks(ticks = np.arange(0, 25, 10))
ax_dendro.set_yticks(ticks = np.arange(0, 25+1, 5), minor = True)

# remove spines
[ax_dendro.spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom']]

# get cell IDs of leaves
leave_cell_IDs = dendrogram_plot['ivl']

# invert order of leave cell IDs
# leave_cell_IDs = list(reversed(leave_cell_IDs))

# resort dataframe indices
cellmorph_descriptors_zscored_clustered = cellmorph_descriptors_zscored.reindex(leave_cell_IDs)
aux_cellmorph_descriptors_clustered = aux_cellmorph_descriptors.reindex(leave_cell_IDs)


# # # heatmap # # #

# colormap for heatmap
c_map_str = 'icefire'

# custom colormaps
# from matplotlib.colors import LinearSegmentedColormap
# c_map_str = sbn.diverging_palette(h_neg=60, h_pos=300, s=75, l=75, center="dark", as_cmap=True)
# c_map_str = LinearSegmentedColormap.from_list('diverging_k', ['#FF00FF', '#000000', '#FFFF00'])


# plot heatmap
sbn.heatmap(cellmorph_descriptors_zscored_clustered.T,
            vmin = heatmin,
            vmax = heatmax,
            square = False,
            xticklabels= False,
            ax = axs_clustering[2], 
            cmap=c_map_str,
            center = 0,
            yticklabels=1,
            linewidth = 0,
            cbar = False) 

# remove ylabel
axs_clustering[2].set_xlabel('')
# axs_clustering[1].set_xticks([])

# axs_clustering[1].xaxis.set_label_position('top') 

### heatmap colorbar ###
ax_cbar = axs_clustering[3]


# create normalize object
norm = mtl.colors.Normalize(-2, 2)

# create mappable cmap object
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=c_map_str)
   
# plot colorbar in axis
cbar = fig_clustering.colorbar(mappable = cmap, 
                                cax = ax_cbar, 
                                label = '', 
                                orientation='vertical',
                                drawedges = False)

# calc new limits
cmin = -2
cmax = 2
crange = (cmax - cmin) * 0.05
cmin_new = cmin - crange
cmax_new = cmax + crange

# apply changes
# y axis
ax_cbar.set_ylim([cmin_new, cmax_new])
ax_cbar.spines['left'].set_bounds([cmin_new, cmax_new])

# x axis
ax_cbar.set_xlim([-3.5, 6.5])
ax_cbar.spines.right.set_position(('data', 1))

# colorbar
cbar.set_ticks(np.arange(cmin, cmax+0.1, 1))
cbar.outline.set_visible(False)
cbar.set_ticklabels(np.arange(cmin, cmax+0.1, 1, dtype = int))

ax_cbar.set_ylabel('Z-scored parameters [std]')
# ax_cbar.yaxis.set_label_position('left')
# ax_cbar.yaxis.set_ticks_position('left')

# get list of auxillary parameters
aux_parameters = aux_cellmorph_descriptors.columns.to_list()

for a_idx, aux_parameter in enumerate(aux_parameters):
    
    # set heat map axis
    ax_heatmap = axs_clustering[axs_dict[aux_parameter]['heatmap']]
    
    # plot heatmap
    heatmap = sbn.heatmap(data = aux_cellmorph_descriptors_clustered.loc[:, [aux_parameter]].T,
                          ax = ax_heatmap,
                          square= False,
                          xticklabels=False,
                          yticklabels=False,
                          linewidth = 0,
                          cbar = False,
                          **heatmap_dict[aux_parameter])
    
    # rotate xticklabels
    ax_heatmap.set_yticks(ticks = [0.5], labels = [aux_parameter]) 
    
    # ylim
    ax_heatmap.set_ylim([0, 1])
    
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
                                        orientation='horizontal',
                                        drawedges = False)
        
        
        
        # colorbar annotation
        
        # get min and max of y axis to shrink subplot
        cmin = cbar_dict[aux_parameter]['range'][0]
        cmax = cbar_dict[aux_parameter]['range'][-1]
        crange = (cmax - cmin) * 0
        
        # calc new limits
        cmin_new = cmin - crange
        cmax_new = cmax + crange
        
        # apply changes
        ax_cbar.set_ylim([0, 1])
        
        for ctick, clabel in zip(cbar_dict[aux_parameter]['ticks'], cbar_dict[aux_parameter]['labels']):
            
            ax_cbar.text(x = ctick,
                          y = 0.5,
                          s = clabel,
                          fontsize = 9,
                          rotation = 0,
                          va = 'center',
                          ha = 'center')
        
        cbar.set_ticks([])
        cbar.outline.set_visible(False)

# delete axis
fig_clustering.delaxes(axs_clustering[1])

# show plot
plt.show()

# save figure
clustering_fig_dir = "C:/Users/nesseler/Desktop/Poster_iBehave"
save_figures(fig_clustering, 'figure-hierarchical_clustering-dendrogram_heatmap_aux-z_scoring-cells_wAxon',
              save_dir = clustering_fig_dir,
              darkmode_bool= darkmode_bool,
              figure_format= 'both')


# # # save dataset
# # cellmorph_descriptors_zscored_clustered.to_excel(join(clustering_fig_dir, 'cellmorph_descriptors-cells_wAxons-z_scored-clustered.xlsx'), index_label = 'cell_ID')
# # aux_cellmorph_descriptors_clustered.to_excel(join(clustering_fig_dir, 'cellmorph_auxillary_descriptors-cells_wAxons-z_scored-clustered.xlsx'), index_label = 'cell_ID')

