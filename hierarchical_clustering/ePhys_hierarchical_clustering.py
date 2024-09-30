#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 16:01:21 2024

@author: moritznesseler
"""

import pandas as pd
import scipy as sc
import numpy as np
from os.path import join

from parameters.directories_win import table_file, hierarchical_dir, cell_morph_descrip_dir



# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# load celldescriptors
celldescriptors = pd.read_excel(join(hierarchical_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')

# load cellmorpho clustering
cellmorph_descriptors = pd.read_excel(join(cell_morph_descrip_dir, 'cellmorpho_celldescriptors_n_aux_n_cluster.xlsx'), index_col = 'cell_ID')

# drop columns
from hierarchical_clustering.ePhys_hierarchical_parameters import parameters_toDrop
celldescriptors.drop(columns = parameters_toDrop, inplace = True)

# %% drop BAOT/MeA (unclassified) cells

cell_IDs = celldescriptors.index
cell_IDs_n_regions = MetaData.loc[cell_IDs, 'Region']
cell_IDs = cell_IDs_n_regions[cell_IDs_n_regions != 'BAOT/MeA']

celldescriptors = celldescriptors.loc[cell_IDs.index, :]

# %% initialize plotting

import matplotlib as mtl
import matplotlib.pyplot as plt
import seaborn as sbn

from functions.functions_plotting import save_figures, get_colors, get_figure_size, apply_axis_settings

from matplotlib.colors import LinearSegmentedColormap

# set colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% z-score celldescriptors

# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()

# sbn.clustermap(data = celldescriptors_zscored,
#                method = 'ward',
#                col_cluster= False,
#                center = 0,
#                cmap = 'icefire',
#                vmin = -2,
#                vmax = 2,
#                tree_kws= {'lw' : 2})

# %% hierarchical clustering

from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import fcluster

# set min and max for heatmap
heatmin = -2
heatmax = 2

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


### elbow plot ###
fig_elbow, ax_elbow = plt.subplots(nrows = 1,
                                   ncols = 1,
                                   layout = 'constrained',
                                   figsize = get_figure_size(width = 165.5),
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
fig_dir = join(hierarchical_dir, 'temp_figs')

# save figure
save_figures(fig_elbow, 'figure-hierarchical_clustering-elbow_plot', 
             save_dir = fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')


# %% 

# get cell_IDs from cellmorph_descriptors
cell_IDs_clustering = celldescriptors_zscored.index.to_list()

# create dataframe for auxillary heatmap
aux_celldescriptors = pd.DataFrame(index = cell_IDs_clustering)

# get auxillary data
aux_celldescriptors['Region'] = MetaData.loc[cell_IDs_clustering, 'Region']

# transform region to values
aux_celldescriptors['Region'] = aux_celldescriptors['Region'].map({'MeA' : 0, 'BAOT/MeA' : 0.5, 'BAOT' : 1})


# set number of clusters
n_clusters = 5

# as halfway point between n_clusters-1 and n_clusters
c_threshold = last_clusters_rev[n_clusters-1] + (last_clusters_rev[n_clusters-2] - last_clusters_rev[n_clusters-1]) / 2

# set cmap
cmap_str = 'icefire'

# calculate width ratio
n_cols_heatmap = celldescriptors_zscored.shape[1]
n_additional_cols = 2
heatmap_width_r = n_cols_heatmap / (n_cols_heatmap + n_additional_cols)
single_col_width_r = 1 / (n_cols_heatmap + n_additional_cols)

# %% set dicts for plotting

heatmap_dict = {'Region'              : {'cmap' : [region_colors[region] for region in ['MeA', 'BAOT']]},
                'cellmorph_category'  : {'cmap' : ['#E6E6E6', '#CCCCCC', '#B3B3B3', '#999999', '#808080', '#666666']}}

cbar_dict = {'Region'              : {'ticks' : [0.25, 0.75],           'labels' : ['MeA', 'BAOT'], 'range' : [0, 1]},
             'cellmorph_category'  : {'ticks' : np.arange(0, 8, 1) + 0.5,      'labels' : ['None', '1', '2', '3', '4', '5', '6'], 'range' : [0, 6]}} 

axs_dict = {'Region'              : {'heatmap' :2, 'cbar' : 5},
            'cellmorph_category'  : {'heatmap' :3, 'cbar' : 6}}


# %% combined heatmap figure

# set linewidth
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams.update()

# initialise figure
fig_clustering, axs_clustering = plt.subplots(nrows = 1,
                                              ncols = 7,
                                              layout = 'constrained',
                                              figsize = get_figure_size(),
                                              width_ratios=[0.2, heatmap_width_r] + [single_col_width_r] * (n_additional_cols) + [0.025] * (n_additional_cols+1),
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


### heatmap colorbar ###
ax_cbar = axs_clustering[4]
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


# get list of auxillary parameters
aux_parameters = aux_celldescriptors_clustered.columns.to_list()

for a_idx, aux_parameter in enumerate(aux_parameters):
    
    # set heat map axis
    ax_heatmap = axs_clustering[axs_dict[aux_parameter]['heatmap']]
    
    # plot heatmap
    heatmap = sbn.heatmap(data = aux_celldescriptors_clustered.loc[:, [aux_parameter]],
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
            norm_min = aux_celldescriptors_clustered.loc[:, [aux_parameter]].min()
            norm_max = aux_celldescriptors_clustered.loc[:, [aux_parameter]].max()
     
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




# show plot
plt.show()

