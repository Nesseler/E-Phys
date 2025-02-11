#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 16:01:21 2024

@author: moritznesseler
"""

from functions.initialize_packages import *

# load celldescriptors
from parameters.directories_win import clustering_dir
celldescriptors = pd.read_excel(join(clustering_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')


# get cell_IDs
cell_IDs = celldescriptors.index.to_list()

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)


# %% z-score celldescriptors

# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()


# %% initialize plotting

from functions.initialize_plotting import *


# %% z-score celldescriptors

# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()


# %% hierarchical clustering

from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import fcluster

# set min and max for heatmap
heatmin = -2.5
heatmax = 2.5

# calc distance between clusters
ward_clustering_linkage = linkage(celldescriptors_zscored, method="ward", metric="euclidean")

# get last few clusters that have been merged by the linkage functions
last_clusters = ward_clustering_linkage[-20:, 2]

# reverse list 
last_clusters_rev = last_clusters[::-1]

# set a list of indices
last_clusters_idc = np.arange(1, len(last_clusters) +1)

# calculate the acceleration of lost distance between clusters
# calculated as the 2nd derivative
acceleration = np.diff(last_clusters_rev, 2)   

# calc change of distance (for complettness)
# calculated as the 1nd derivative
change = np.diff(last_clusters_rev, 1)   


# %% elbow plot

fig_elbow, ax_elbow = plt.subplots(nrows = 1,
                                   ncols = 1,
                                   layout = 'constrained',
                                   figsize = get_figure_size(height = 120, width = 120),
                                   dpi = 600)

# set axis title
ax_elbow.set_title('hierarchical clustering z-transformed - elbow plot',
                    fontsize = 12)

# plot cluster distances
ax_elbow.plot(last_clusters_idc, last_clusters_rev,
              lw = 1,
              c = colors_dict['primecolor'],
              label = 'Cluster distance')

# plot change
ax_elbow.plot(last_clusters_idc[:-1] + 1, change,
              lw = 1,
              ls = 'dashed',
              c = colors_dict['color3'],
              alpha = 0.5,
              label = 'Change of cluster distance')

# plot acceleration
ax_elbow.plot(last_clusters_idc[:-2] + 1, acceleration,
              lw = 1,
              c = colors_dict['color3'],
              label = 'Acceleration of cluster\ndistance growth')

# edit axis
# y
ydict = {'ax_min' : -10,
         'ax_max' : 25,
         'pad' : 1,
         'step' : 5,
         'stepminor' : 1,
         'label' : 'Distance between clusters [pA]'}

# edit axis
apply_axis_settings(ax_elbow, axis = 'y', **ydict)

# x
xdict = {'ax_min' : 0,
         'ax_max' : 20,
         'pad' : 0.5,
         'step' : 5,
         'stepminor' : 1,
         'label' : 'Number of clusters [#]'}

# edit axis
apply_axis_settings(ax_elbow, axis = 'x', **xdict)

# edit legend
ax_elbow.legend(prop={'size': 9})
ax_elbow.get_legend().get_frame().set_linewidth(0.0)

# remove spines
[ax_elbow.spines[spine].set_visible(False) for spine in ['top', 'right']]

# show plot
plt.show()

# set directory for figure
fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_elbow, 'figure-hierarchical_clustering-elbow_plot', 
             save_dir = fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')


# %% auxillary data

# create dataframe for auxillary heatmap
aux_celldescriptors = pd.DataFrame(index = cell_IDs)

# get auxillary data
aux_celldescriptors['Region'] = MetaData.loc[cell_IDs, 'Region']

# transform region to values
aux_celldescriptors['Region'] = aux_celldescriptors['Region'].map({'MeA' : 0, 'BAOT/MeA' : 0.5, 'BAOT' : 1})

# set number of clusters
n_clusters = 6

# as halfway point between n_clusters-1 and n_clusters
c_threshold = last_clusters_rev[n_clusters-1] + (last_clusters_rev[n_clusters-2] - last_clusters_rev[n_clusters-1]) / 2

# set cmap
cmap_str = 'icefire'

# calculate width ratio
n_cols_heatmap = celldescriptors_zscored.shape[1]
n_additional_cols = 1
heatmap_width_r = n_cols_heatmap / (n_cols_heatmap + n_additional_cols)
single_col_width_r = 1 / (n_cols_heatmap + n_additional_cols)


# %% set dicts for plotting

heatmap_dict = {'Region'              : {'cmap' : [region_colors[region] for region in ['MeA', 'BAOT/MeA', 'BAOT']]}}

cbar_dict =    {'Region'              : {'ticks'  : [1/6, 3/6, 5/6],
                                         'labels' : ['MeA', 'BAOT/MeA', 'BAOT'], 'range' : [0, 1]}} 

axs_dict =     {'Region'              : {'heatmap' :2, 'cbar' : 4}}


# %% combined heatmap figure

# set linewidth
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams.update()

# initialise figure
fig_clustering, axs_clustering = plt.subplots(nrows = 1,
                                              ncols = 5,
                                              layout = 'constrained',
                                              figsize = get_figure_size(),
                                              width_ratios=[0.2, heatmap_width_r] + [single_col_width_r] + [0.025] + [0.025],
                                              dpi = 600)

# adjust layout of constrained setup
fig_clustering.set_constrained_layout_pads(w_pad=0.01, wspace=0.05)

### dendrogram ###
# set axis for dendrogram
ax_dendro = axs_clustering[0]

# set color palette for dendrogram
hierarchy.set_link_color_palette(None)

# plot dendrogram
dendrogram_plot = dendrogram(Z = ward_clustering_linkage, 
                              labels = celldescriptors_zscored.index.to_list(), 
                              ax = ax_dendro,
                              orientation = 'left',
                              color_threshold = c_threshold)

# plot cluster threshold
ax_dendro.axvline(x = c_threshold, 
                  color = colors_dict['primecolor'],
                  lw = 1,
                  ls = 'dashed')

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
celldescriptors_zscored_clustered = celldescriptors_zscored.reindex(leave_cell_IDs)
aux_celldescriptors_clustered = aux_celldescriptors.reindex(leave_cell_IDs)


# # # heatmap # # #
ax_heat = axs_clustering[1]

# plot heatmap
sbn.heatmap(celldescriptors_zscored_clustered,
            vmin = heatmin,
            vmax = heatmax,
            center = 0,
            square = False,
            xticklabels = 1,
            ax = ax_heat, 
            cmap = cmap_str, 
            yticklabels = False,
            linewidth = 0,
            cbar = False) 

# remove ylabel
ax_heat.set_ylabel('')


# # # heatmap colorbar # # #
ax_cbar = axs_clustering[3]
cmin = heatmin
cmax = heatmax
crange = (cmax - cmin) * 0.3

# create normalize object
norm = mtl.colors.Normalize(heatmin, heatmax)

# create mappable cmap object
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
   
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


# # # region auxillary parameter # # #

# set heat map axis
ax_heatmap = axs_clustering[axs_dict['Region']['heatmap']]

# plot heatmap
heatmap = sbn.heatmap(data = aux_celldescriptors_clustered.loc[:, ['Region']],
                      ax = ax_heatmap,
                      square= False,
                      xticklabels=1,
                      yticklabels=False,
                      linewidth = 0,
                      cbar = False,
                      **heatmap_dict['Region'])

# rotate xticklabels
[label.set_rotation(90) for label in ax_heatmap.get_xticklabels()] 

# set colorbar axis
ax_cbar = axs_clustering[axs_dict['Region']['cbar']]

n_colors = len(heatmap_dict['Region']['cmap'])
cmap = LinearSegmentedColormap.from_list('Region', heatmap_dict['Region']['cmap'], N=n_colors)

# # min max normalize time for color-code
norm_min = aux_celldescriptors_clustered.loc[:, ['Region']].min()
norm_max = aux_celldescriptors_clustered.loc[:, ['Region']].max()

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
cmin = 0
cmax = 1
crange = (cmax - cmin) * 0.3

# calc new limits
cmin_new = cmin - crange
cmax_new = cmax + crange

# apply changes
ax_cbar.set_ylim([cmin_new, cmax_new])

for ctick, clabel in zip(cbar_dict['Region']['ticks'], cbar_dict['Region']['labels']):

    ax_cbar.text(x = 0.55,
                  y = ctick,
                  s = clabel,
                  fontsize = 7,
                  rotation = 90,
                  va = 'center',
                  ha = 'center')

cbar.set_ticks([])
cbar.outline.set_visible(False)


# show plot
plt.show()

# set directory for figure
fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_clustering, 'figure-hierarchical_clustering-dendro_heat_aux', 
              save_dir = fig_dir,
              darkmode_bool = darkmode_bool,
              figure_format = 'png')


# %% ##################


# set linewidth
plt.rcParams['lines.linewidth'] = 0.0
plt.rcParams.update()

# initialise figure
fig_clustering, axs_clustering = plt.subplots(nrows = 1,
                                              ncols = 1,
                                              layout = 'constrained',
                                              figsize = get_figure_size(height = 160.5, width = 1000),
                                              # height_ratios=[0.2, heatmap_width_r],
                                              dpi = 600)

# adjust layout of constrained setup
fig_clustering.set_constrained_layout_pads(w_pad=0.01, wspace=0.05)

### dendrogram ###
# set axis for dendrogram
# ax_dendro = axs_clustering[0]

# set color palette for dendrogram
hierarchy.set_link_color_palette(None)

# # plot dendrogram
# dendrogram_plot = dendrogram(Z = ward_clustering_linkage, 
#                               labels = celldescriptors_zscored.index.to_list(), 
#                               ax = ax_dendro,
#                               orientation = 'top',
#                               color_threshold = c_threshold)

# # plot cluster threshold
# ax_dendro.axhline(y = c_threshold, 
#                   color = colors_dict['primecolor'],
#                   lw = 1,
#                   ls = 'dashed')

# # edit axis
# ax_dendro.set_ylabel('Distance')
# ax_dendro.set_yticks(ticks = np.arange(0, 25, 10))
# ax_dendro.set_yticks(ticks = np.arange(0, 25+1, 5), minor = True)

# # remove spines
# [ax_dendro.spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom']]

# get cell IDs of leaves
leave_cell_IDs = dendrogram_plot['ivl']

# invert order of leave cell IDs
leave_cell_IDs = list(leave_cell_IDs)

# resort dataframe indices
celldescriptors_zscored_clustered = celldescriptors_zscored.reindex(leave_cell_IDs)
aux_celldescriptors_clustered = aux_celldescriptors.reindex(leave_cell_IDs)


# # # heatmap # # #
ax_heat = axs_clustering

# plot heatmap
sbn.heatmap(celldescriptors_zscored_clustered.T,
            vmin = heatmin,
            vmax = heatmax,
            center = 0,
            square = False,
            xticklabels = False,
            ax = ax_heat, 
            cmap = cmap_str, 
            # xticklabels = False,
            linewidth = 0,
            cbar = False) 

# remove ylabel
ax_heat.set_ylabel('')


# # # # heatmap colorbar # # #
# ax_cbar = axs_clustering[3]
# cmin = heatmin
# cmax = heatmax
# crange = (cmax - cmin) * 0.3

# # create normalize object
# norm = mtl.colors.Normalize(heatmin, heatmax)

# # create mappable cmap object
# cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)
   
# # plot colorbar in axis
# cbar = fig_clustering.colorbar(mappable = cmap, 
#                                 cax = ax_cbar, 
#                                 label = '', 
#                                 orientation='vertical',
#                                 drawedges = False)

# # calc new limits
# cmin_new = cmin - crange
# cmax_new = cmax + crange

# # apply changes
# # axis
# ax_cbar.set_ylim([cmin_new, cmax_new])
# ax_cbar.spines['left'].set_bounds([cmin_new, cmax_new])

# # colorbar
# cbar.set_ticks(np.arange(cmin, cmax+0.1, 1))
# cbar.outline.set_visible(False)
# cbar.set_ticklabels(np.arange(cmin, cmax+0.1, 1, dtype = int))

# ax_cbar.set_ylabel('Z-scored parameters [std]')
# ax_cbar.yaxis.set_label_position('left')


# # # # region auxillary parameter # # #

# # set heat map axis
# ax_heatmap = axs_clustering[axs_dict['Region']['heatmap']]

# # plot heatmap
# heatmap = sbn.heatmap(data = aux_celldescriptors_clustered.loc[:, ['Region']],
#                       ax = ax_heatmap,
#                       square= False,
#                       xticklabels=1,
#                       yticklabels=False,
#                       linewidth = 0,
#                       cbar = False,
#                       **heatmap_dict['Region'])

# # rotate xticklabels
# [label.set_rotation(90) for label in ax_heatmap.get_xticklabels()] 

# # set colorbar axis
# ax_cbar = axs_clustering[axs_dict['Region']['cbar']]

# n_colors = len(heatmap_dict['Region']['cmap'])
# cmap = LinearSegmentedColormap.from_list('Region', heatmap_dict['Region']['cmap'], N=n_colors)

# # # min max normalize time for color-code
# norm_min = aux_celldescriptors_clustered.loc[:, ['Region']].min()
# norm_max = aux_celldescriptors_clustered.loc[:, ['Region']].max()

# # create normalize object
# norm = mtl.colors.Normalize(norm_min, norm_max)

# # create mappable cmap object
# cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap)
   
# # plot colorbar in axis
# cbar = fig_clustering.colorbar(mappable = cmap, 
#                                 cax = ax_cbar, 
#                                 label = '', 
#                                 orientation='vertical',
#                                 drawedges = False)

# # colorbar annotation

# # get min and max of y axis to shrink subplot
# cmin = 0
# cmax = 1
# crange = (cmax - cmin) * 0.3

# # calc new limits
# cmin_new = cmin - crange
# cmax_new = cmax + crange

# # apply changes
# ax_cbar.set_ylim([cmin_new, cmax_new])

# for ctick, clabel in zip(cbar_dict['Region']['ticks'], cbar_dict['Region']['labels']):

#     ax_cbar.text(x = 0.55,
#                   y = ctick,
#                   s = clabel,
#                   fontsize = 7,
#                   rotation = 90,
#                   va = 'center',
#                   ha = 'center')

# cbar.set_ticks([])
# cbar.outline.set_visible(False)


# show plot
plt.show()

# set directory for figure
fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_clustering, 'figure-hierarchical_clustering-onlyflipped_heatmap', 
              save_dir = fig_dir,
              darkmode_bool = darkmode_bool,
              figure_format = 'both')



############################

# %% get cell_IDs and cluster

# cluster the original parameter
celldescriptors_clustered = celldescriptors.reindex(leave_cell_IDs)

# assign cluster id to cells
celldescriptors_clustered['hierarchical_cluster'] = [int(clus.replace('C', '')) for clus in dendrogram_plot['leaves_color_list'][::-1]]

# save dataframe
celldescriptors_clustered.to_excel(join(clustering_dir, 'ePhys_celldescriptors-clustered.xlsx'), 
                                   index_label = 'cell_ID')

