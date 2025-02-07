# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:05:10 2025

@author: nesseler
"""

from functions.initialize_packages import *

# load celldescriptors
from parameters.directories_win import clustering_dir
celldescriptors_clustered = pd.read_excel(join(clustering_dir, 'ePhys_celldescriptors-clustered.xlsx'), index_col = 'cell_ID')

# get cell_IDs
cell_IDs = celldescriptors_clustered.index.to_list()

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# add region to dataframe
celldescriptors_clustered['Region'] = MetaData.loc[cell_IDs, 'Region']


# %% initialize plotting

from functions.initialize_plotting import *


# %% violinplot parameter

n_cluster = celldescriptors_clustered['hierarchical_cluster'].max()

for parameter in tqdm(celldescriptors_clustered.drop(columns=['Region', 'hierarchical_cluster']).columns.to_list()):
    
    # initialize figure
    fig, ax = plt.subplots(nrows = 1,
                           ncols = 1,
                           layout = 'constrained',
                           figsize = get_figure_size(width = 200, height = 100),
                           dpi = 600)
    
    # set plotting dicts
    swarm_dict = {'size' : 4,
                  'linewidth' : 0.5,
                  'edgecolor' : colors_dict['seccolor']}
    
    # iterate through cluster
    for cluster in range(n_cluster+1):
        
        # get cell_IDs for cluster
        cluster_cell_IDs = celldescriptors_clustered[celldescriptors_clustered['hierarchical_cluster'] == cluster].index.to_list()
    
        # get data and x position
        data = celldescriptors_clustered.loc[cluster_cell_IDs, parameter]
        x = [cluster] * len(celldescriptors_clustered.loc[cluster_cell_IDs, parameter])
        
        # get assigned region for each cell in cluster
        cluster_regions = MetaData.loc[cluster_cell_IDs, 'Region'].to_list()
        
        # assign colors to dots
        colors = [region_colors[c_region] for c_region in cluster_regions]
    
        # plot data points as swarm
        sbn.swarmplot(y = data,
                      x = x,
                      ax = ax,
                      c = colors,
                      **swarm_dict)
        
        # width relative to total number of cell
        v_width = 0.25
        
        if len(data) > 2:
            # plot half violin
            plot_half_violin(data = data.to_list(), 
                             ax = ax,
                             v_direction = 1,
                             v_resolution = 0.01,
                             v_kde_cutoff = 0.1,
                             v_position = cluster,
                             v_offset = 0.1,
                             v_width = v_width,
                             v_baseline = False,
                             v_color = colors_dict['primecolor'],
                             v_zorder = 0,
                             v_lw = 0.75)
            
            e_position = cluster + (0.15 * 1)
        
            # plot errorbar
            ax.errorbar(x = e_position,
                        y = data.mean(),
                        yerr = data.std(),
                        fmt='_', 
                        markersize = 4,
                        markerfacecolor = 'none',
                        capsize = 2,
                        color = colors_dict['primecolor'],
                        linewidth = 0.75,
                        label = '_nolegend_',
                        zorder = 2)
            
            # plot median
            ax.scatter(x = e_position,
                       y = data.median(),
                       marker='D', 
                       s = 4,
                       color = colors_dict['primecolor'],
                       linewidth = 1,
                       label = '_nolegend_',
                       zorder = 3)
            
    # edit axis
    apply_axis_settings(ax, axis = 'y', **axis_dict[parameter])
    
    # x
    xdict = {'ax_min' : 0,
             'ax_max' : n_cluster,
             'pad' : 0.5,
             'step' : 1,
             'stepminor' : 1,
             'label' : 'Cluster'}
    
    apply_axis_settings(ax, axis = 'x', **xdict)
            
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # align labels
    fig.align_labels()
    
    # create saving path and save
    save_figures(fig, f'vplot-cluster_single_parameters-{parameter}', join(clustering_dir, 'temp_figs'), darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()