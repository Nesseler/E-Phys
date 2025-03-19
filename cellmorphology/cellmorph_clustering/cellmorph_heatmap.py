# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:07:51 2025

@author: nesseler
"""


from functions.initialize_packages import *


# %% load dataframes to describe cells

from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_metrics_dir, cellmorph_analysis_dir

cellmetrics = pd.read_excel(join(cellmorph_metrics_dir, 'cellmorph_cellmetrics' + '.xlsx'), index_col = 'cell_ID')
cellmetrics_aux = pd.read_excel(join(cellmorph_metrics_dir, 'cellmorph_cellmetrics_auxillary' + '.xlsx'), index_col = 'cell_ID')

cell_IDs = cellmetrics.index.to_list()


# %% transform

# z-score cellmorph matrix
cellmetrics_zscored = (cellmetrics - cellmetrics.mean()) / cellmetrics.std()


# %% initialize plotting

# init plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *


# %% hierarchical clustering

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import fcluster

# set variables
n_clusters = 7
heatmin = -2
heatmax = 2
cbar_label = 'Z-scored parameter value [std]'
c_threshold = 15

# calc distance between clusters
ward_clustering_linkage = linkage(cellmetrics_zscored, method="ward", metric="euclidean")

# get last few clusters that have been merged by the linkage functions
last_clusters = ward_clustering_linkage[-20:, 2]

# reverse list 
last_clusters_rev = last_clusters[::-1]



# %% remap categorical values

# transform region to values
cellmetrics_aux['Region'] = cellmetrics_aux['Region'].map({'MeA' : 0, 'BAOT' : 1})
cellmetrics_aux['AIS root'] = cellmetrics_aux['AIS root'].map({'dendritic' : 0, 'somatic' : 1})
cellmetrics_aux['spinyness'] = cellmetrics_aux['spinyness'].map({'low' : 0, 'moderate' : 1, 'high' : 2})


# %% heatmap

heatmap_dict = {'Region'                  : {'cmap' : [region_colors[region] for region in ['MeA', 'BAOT']]},
                'spinyness'               : {'cmap' : [spines_color_dict['both'][spinyness] for spinyness in ['low', 'moderate', 'high']]},
                'dendrites-circular_mean' : {'cmap' : plt.get_cmap('twilight'), 'vmin' : 0, 'vmax' : np.pi * 2},
                'axons-circular_mean'     : {'cmap' : plt.get_cmap('twilight'), 'vmin' : 0, 'vmax' : np.pi * 2},
                'AIS root'                : {'cmap' : [neurite_color_dict[region]['axons'] for region in ['MeA', 'BAOT']]}}

cbar_dict = {'Region'                  : {'ticks' : [0.25, 0.75],                         'labels' : ['MeA', 'BAOT'],             'range' : [0, 1]},
             'spinyness'               : {'ticks' : [0.333, 1.0, 1.666],                  'labels' : ['low', 'moderate', 'high'], 'range' : [0, 2]},
             'dendrites-circular_mean' : {'ticks' : np.arange(0, np.pi *2 +0.1, np.pi/2), 'labels' : ['p', 'd', 'a', 'v', 'p'],   'range' : [0, np.pi * 2]},
             'axons-circular_mean'     : {'ticks' : np.arange(0, np.pi *2 +0.1, np.pi/2), 'labels' : ['p', 'd', 'a', 'v', 'p'],   'range' : [0, np.pi * 2]},
             'AIS root'                : {'ticks' : [0.25, 0.75],                         'labels' : ['dendritic', 'somatic'],    'range' : [0, 1]}} 

axs_dict = {'Region'                  : {'heatmap' :3, 'cbar' : 10},
            'spinyness'               : {'heatmap' :4, 'cbar' : 11}, 
            'dendrites-circular_mean' : {'heatmap' :5, 'cbar' : np.nan},
            'axons-circular_mean'     : {'heatmap' :6, 'cbar' : 12},
            'AIS root'                : {'heatmap' :7, 'cbar' : 13}}


# as halfway point between n_clusters-1 and n_clusters
c_threshold = last_clusters_rev[n_clusters-1] + (last_clusters_rev[n_clusters-2] - last_clusters_rev[n_clusters-1]) / 2

# calculate width ratio
n_cols_heatmap = cellmetrics_zscored.shape[1]
n_additional_cols = 10
heatmap_width_r = n_cols_heatmap / (n_cols_heatmap + n_additional_cols)
single_col_width_r = 1 / (n_cols_heatmap + n_additional_cols)




# initialise figure
fig_clustering, axs_clustering = plt.subplots(nrows = 1,
                                              ncols = 14,
                                              layout = 'constrained',
                                              figsize = get_figure_size(width = 260.334, height = 295.47),
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
                             labels = cellmetrics_zscored.index.to_list(), 
                             ax = ax_dendro,
                             orientation = 'left',
                             color_threshold = c_threshold)

# plot cluster threshold
ax_dendro.axvline(x = c_threshold, 
                  color = colors_dict['primecolor'],
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
cellmetrics_zscored_clustered = cellmetrics_zscored.reindex(leave_cell_IDs)
cellmetrics_aux_clustered = cellmetrics_aux.reindex(leave_cell_IDs)


### heatmap ###
# plot heatmap
sbn.heatmap(cellmetrics_zscored_clustered,
            vmin = heatmin,
            vmax = heatmax,
            square = False,
            xticklabels= 1,
            ax = axs_clustering[1], 
            cmap="coolwarm", 
            yticklabels=False,
            linewidth = 0,
            cbar = False) 


# rewrite variables for labels


# remove ylabel
axs_clustering[1].set_ylabel('')
axs_clustering[1].set_xticks(ticks = np.arange(0, cellmetrics_zscored_clustered.shape[1], 1) + 0.5,
                             labels = [c.replace('dendrites-', '').replace('axons-', '') for c in cellmetrics_zscored_clustered.columns.to_list()])

# secondary axis
sec = axs_clustering[1].secondary_xaxis(location=-0.25)
apply_axis_settings(sec, axis = 'x',  **{'ax_min' : 0.5, 'ax_max' : 8.5,  'pad' : None, 'step' : 8, 'stepminor' : 8, 'label' : ''})
sec.set_xticks(ticks = [0.5, 4.5, 8.5], labels = ['', 'Dendrites', ''])

sec2 = axs_clustering[1].secondary_xaxis(location=-0.25)
apply_axis_settings(sec2, axis = 'x', **{'ax_min' : 9.5, 'ax_max' : 16.5, 'pad' : None, 'step' : 8, 'stepminor' : 8, 'label' : ''})
sec2.set_xticks(ticks = [9.5, 13, 16.5], labels = ['', 'Axons', ''])


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
aux_parameters = cellmetrics_aux_clustered.columns.to_list()

for a_idx, aux_parameter in enumerate(aux_parameters):
    
    # set heat map axis
    ax_heatmap = axs_clustering[axs_dict[aux_parameter]['heatmap']]
    
    # plot heatmap
    heatmap = sbn.heatmap(data = cellmetrics_aux_clustered.loc[:, [aux_parameter]],
                          ax = ax_heatmap,
                          square= False,
                          xticklabels=1,
                          yticklabels=False,
                          linewidth = 0,
                          cbar = False,
                          **heatmap_dict[aux_parameter])
    
    # remove ylabel
    ax_heatmap.set_ylabel('')
    
    # rotate xticklabels
    [label.set_rotation(90) for label in ax_heatmap.get_xticklabels()] 
    
    if type(axs_dict[aux_parameter]['cbar']) == int:
        # set colorbar axis
        ax_cbar = axs_clustering[axs_dict[aux_parameter]['cbar']]
        
        # create color map from list of colors
        if 'circular_mean' not in aux_parameter: 
            n_colors = len(heatmap_dict[aux_parameter]['cmap'])
            cmap = LinearSegmentedColormap.from_list(aux_parameter, heatmap_dict[aux_parameter]['cmap'], N=n_colors)
        
            # # min max normalize time for color-code
            norm_min = cellmetrics_aux_clustered.loc[:, [aux_parameter]].min()
            norm_max = cellmetrics_aux_clustered.loc[:, [aux_parameter]].max()
     
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
        # ax_cbar.set_ylabel('')
        
        for ctick, clabel in zip(cbar_dict[aux_parameter]['ticks'], cbar_dict[aux_parameter]['labels']):
            
            ax_cbar.text(x = 0.55,
                         y = ctick,
                         s = clabel,
                         fontsize = 12,
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
                                     fontsize = 14,
                                     loc = 'left')

# show plot
plt.show()

# save figure
save_figures(fig_clustering, f'hierarchical_clustering-heatmap', 
             save_dir = join(cellmorph_analysis_dir, 'plots-hierarchical_clustering'),
             darkmode_bool = darkmode_bool,
             figure_format = 'both')