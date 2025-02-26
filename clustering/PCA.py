# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 17:05:06 2025

@author: nesseler

ressource
https://www.geo.fu-berlin.de/en/v/soga-py/Advanced-statistics/Multivariate-Approaches/Principal-Component-Analysis/PCA-the-basics/Choose-Principal-Components/index.html

"""

from functions.initialize_packages import *

# specific package
from sklearn.decomposition import PCA

# load celldescriptors
from parameters.directories_win import clustering_dir
celldescriptors = pd.read_excel(join(clustering_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')
celldescriptors_clustered = pd.read_excel(join(clustering_dir, 'ePhys_celldescriptors-clustered.xlsx'), index_col = 'cell_ID')


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


# %% principal component analysis

celldescriptors_PCA = PCA().fit(celldescriptors_zscored)

celldescriptors_PCA_components = celldescriptors_PCA.transform(celldescriptors_zscored)

celldescriptors_PCA_explained_variance = celldescriptors_PCA.explained_variance_

# # get the eigenvectors and eigenvalues
#     # index: celldescriptors
#     # columns: principal components
celldescriptors_PCA_eigen = pd.DataFrame(celldescriptors_PCA.components_.T,
                                         columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])],
                                         index = celldescriptors.columns)


prinicpal_components = pd.DataFrame(celldescriptors_PCA_components,
                                    index = cell_IDs,
                                    columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])])


# colors = ['#FFEC9DFF', '#FAC881FF', '#F4A464FF', '#E87444FF', '#D9402AFF',
#         '#BF2729FF', '#912534FF', '#64243EFF', '#3D1B28FF', '#161212FF']

colors = [(0.5529411764705883, 0.8274509803921568, 0.7803921568627451),
         (0.996078431372549, 1.0, 0.7019607843137254),
         (0.7490196078431373, 0.7333333333333333, 0.8509803921568627),
         (0.9803921568627451, 0.5058823529411764, 0.4549019607843137),
         (0.5058823529411764, 0.6941176470588235, 0.8235294117647058),
         (0.9921568627450981, 0.7058823529411765, 0.3843137254901961),
         (0.7019607843137254, 0.8705882352941177, 0.4117647058823529),
         (0.7372549019607844, 0.5098039215686274, 0.7411764705882353),
         (0.8, 0.9215686274509803, 0.7686274509803922),
         (1.0, 0.9294117647058824, 0.43529411764705883)]

c_colors = [colors[i] for i in celldescriptors_clustered.loc[cell_IDs, 'hierarchical_cluster'].to_list()]


# get number of clusters
n_clusters = celldescriptors_clustered.loc[cell_IDs, 'hierarchical_cluster'].nunique()

# %% pca and vectors

fig_PCA2, axs_PCA2 = plt.subplots(nrows = 1,
                                 ncols = 2,
                                 layout = 'constrained',
                                 figsize = get_figure_size(),
                                 dpi = 600)

# plot zero lines
for ax in axs_PCA2:
    
    axline_dict = {'linewidth' : 1, 
                   'color' : 'w',
                   'zorder' : 0,
                   'linestyle' : 'dashed',
                   'dashes' : [4, 7],
                   'alpha' : 0.75}
    
    ax.axline((0, 1), (0, 2), **axline_dict)
    ax.axline((1, 0), (2, 0), **axline_dict)

# plot each cluster separatly for coloring
for cluster_idx in range(n_clusters):
    
    # get cell_IDs of cells in cluster
    cluster_cellIDs = celldescriptors_clustered[celldescriptors_clustered['hierarchical_cluster'] == cluster_idx].index.to_list()

    # plot cluster
    axs_PCA2[0].scatter(x = prinicpal_components.loc[cluster_cellIDs, 'PC1'],
                        y = prinicpal_components.loc[cluster_cellIDs, 'PC2'],
                        c = [colors[cluster_idx]] * len(cluster_cellIDs),
                        s = 60,
                        label = cluster_idx)

# plot cell_IDs
for cell_ID in prinicpal_components.index.to_list():
    axs_PCA2[0].text(x = prinicpal_components.at[cell_ID, 'PC1'],
                     y = prinicpal_components.at[cell_ID, 'PC2'],
                     s = cell_ID,
                     color = "k",
                     ha = "center",
                     va = "center",
                     fontsize = 3)

# legend
h, l = axs_PCA2[0].get_legend_handles_labels() 
axs_PCA2[0].legend(h, l, title = 'cluster id', 
                   frameon = False, ncol = 2, loc = 'upper left')

# plot the variables as vectors
for i in range(celldescriptors_PCA.components_.shape[1]):
    axs_PCA2[1].arrow(x = 0,
                     y = 0,
                     dx = celldescriptors_PCA.components_[0, i]*8.5,
                     dy = celldescriptors_PCA.components_[1, i]*8.5,
                     head_width=0.05,
                     head_length=0.05,
                     linewidth=0.75,
                     color="w",
                     )
    
# Plot annotations
for parameter in celldescriptors_PCA_eigen.index.to_list():
    axs_PCA2[1].text(x = (celldescriptors_PCA_eigen.loc[parameter, 'PC1'])*10,
                     y = (celldescriptors_PCA_eigen.loc[parameter, 'PC2'])*10,
                     s = parameter,
                     color="w",
                     fontsize = 9,
                     ha = 'center')

PC_dict = {'ax_min' : -6,
           'ax_max' : 8.5,
           'pad' : None,
           'step' : 2,
           'stepminor' : 0.5,
           'label' : None}

# edit axis
[apply_axis_settings(ax, axis = 'x', **PC_dict) for ax in axs_PCA2]
[apply_axis_settings(ax, axis = 'y', **PC_dict) for ax in axs_PCA2]

for ax in axs_PCA2:
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs_PCA2 for spine in ['top', 'right']]

# set directory for figure
fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_PCA2, 'figure-PCA-components_plot', 
             save_dir = fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'both')

# show plot
plt.show()


# %% regions

fig_PCA_region, axs_PCA_region = plt.subplots(nrows = 1,
                                              ncols = 2,
                                              layout = 'constrained',
                                              figsize = get_figure_size(),
                                              dpi = 600)

# plot zero lines
for ax in axs_PCA_region:
    
    axline_dict = {'linewidth' : 1, 
                   'color' : 'w',
                   'zorder' : 0,
                   'linestyle' : 'dashed',
                   'dashes' : [4, 7],
                   'alpha' : 0.75}
    
    ax.axline((0, 1), (0, 2), **axline_dict)
    ax.axline((1, 0), (2, 0), **axline_dict)

# plot cells with region
for region in ['BAOT', 'BAOT/MeA', 'MeA']:
    
    # get cell_IDs per region
    region_cellIDs = MetaData[MetaData['Region'] == region].index.to_list()

    # plot region cells
    axs_PCA_region[0].scatter(x = prinicpal_components.loc[region_cellIDs, 'PC1'],
                              y = prinicpal_components.loc[region_cellIDs, 'PC2'],
                              c = [region_colors[region]] * len(region_cellIDs),
                              s = 60,
                              label = cluster_idx)

# plot cell_IDs
for cell_ID in prinicpal_components.index.to_list():
    axs_PCA_region[0].text(x = prinicpal_components.at[cell_ID, 'PC1'],
                     y = prinicpal_components.at[cell_ID, 'PC2'],
                     s = cell_ID,
                     color = "k",
                     ha = "center",
                     va = "center",
                     fontsize = 3)

# legend
h, l = axs_PCA_region[0].get_legend_handles_labels() 
axs_PCA_region[0].legend(h, ['BAOT', 'BAOT/MeA', 'MeA'], 
                         title = 'Region', 
                         frameon = False, 
                         ncol = 1, 
                         loc = 'upper left')


# plot the variables as vectors
for i in range(celldescriptors_PCA.components_.shape[1]):
    axs_PCA_region[1].arrow(x = 0,
                     y = 0,
                     dx = celldescriptors_PCA.components_[0, i]*8.5,
                     dy = celldescriptors_PCA.components_[1, i]*8.5,
                     head_width=0.05,
                     head_length=0.05,
                     linewidth=0.75,
                     color="w",
                     )
    

# Plot annotations
for parameter in celldescriptors_PCA_eigen.index.to_list():
    axs_PCA_region[1].text(x = (celldescriptors_PCA_eigen.loc[parameter, 'PC1'])*10,
                     y = (celldescriptors_PCA_eigen.loc[parameter, 'PC2'])*10,
                     s = parameter,
                     color="w",
                     fontsize = 9,
                     ha = 'center'
                     )


PC_dict = {'ax_min' : -6,
           'ax_max' : 8.5,
           'pad' : None,
           'step' : 2,
           'stepminor' : 0.5,
           'label' : None}

# edit axis
[apply_axis_settings(ax, axis = 'x', **PC_dict) for ax in axs_PCA_region]
[apply_axis_settings(ax, axis = 'y', **PC_dict) for ax in axs_PCA_region]

for ax in axs_PCA_region:
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs_PCA_region for spine in ['top', 'right']]

# set directory for figure
fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_PCA_region, 'figure-PCA-components_plot-region', 
             save_dir = fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'both')

# show plot
plt.show()


# %% single parameter

# set parameter
parameter = 'r_input'
parameter_i = celldescriptors.columns.to_list().index(parameter)

# initialize figure
fig_PCA2, axs_PCA2 = plt.subplots(nrows = 1,
                                  ncols = 2,
                                  layout = 'constrained',
                                  figsize = get_figure_size(),
                                  dpi = 600)

# plot zero lines
for ax in axs_PCA2:
    
    axline_dict = {'linewidth' : 1, 
                   'color' : 'w',
                   'zorder' : 0,
                   'linestyle' : 'dashed',
                   'dashes' : [4, 7],
                   'alpha' : 0.75}
    
    ax.axline((0, 1), (0, 2), **axline_dict)
    ax.axline((1, 0), (2, 0), **axline_dict)

# set first axis
ax = axs_PCA2[0]

# plot scatter
sbn.scatterplot(data = prinicpal_components,
                x = 'PC1',
                y = 'PC2',
                hue = celldescriptors[parameter],
                ax = ax,
                linewidth = 0,
                s = 60)

sbn.move_legend(ax, "upper left",
                frameon = False,
                fontsize = 9,
                title = parameter,
                ncol = 1,
                title_fontsize = 9,
                columnspacing = 0.6,
                handletextpad = 0.0)

ax = axs_PCA2[1]

# plot the variables as vectors
ax.arrow(x = 0,
         y = 0,
         dx = celldescriptors_PCA.components_[0, parameter_i]*8.5,
         dy = celldescriptors_PCA.components_[1, parameter_i]*8.5,
         head_width=0.05,
         head_length=0.05,
         linewidth=0.75,
         color="w",
         )
  
# Plot annotations
ax.text(x = (celldescriptors_PCA_eigen.loc[parameter, 'PC1'])*10,
        y = (celldescriptors_PCA_eigen.loc[parameter, 'PC2'])*10,
        s = parameter,
        color="w",
        fontsize = 12,
        ha = 'center'
        )


# edit axis
for ax in axs_PCA2:
    apply_axis_settings(ax, axis = 'x', **PC_dict)
    apply_axis_settings(ax, axis = 'y', **PC_dict)
    
    # set axis labels
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# set directory for figure
fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_PCA2, f'figure-PCA-{parameter}', 
              save_dir = fig_dir,
              darkmode_bool = darkmode_bool,
              figure_format = 'both')

# show plot
plt.show()

# %%


# # %% single parameter color code

# parameter = 'r_input'

# for parameter in celldescriptors_PCA_eigen.index.to_list():

#     parameter_i = celldescriptors_PCA_eigen.index.to_list().index(parameter)
    
#     fig_PCA3, axs_PCA3 = plt.subplots(nrows = 1,
#                                      ncols = 2,
#                                      layout = 'constrained',
#                                      figsize = get_figure_size(height = 120, width = 240),
#                                      dpi = 600)
    
#     sbn.scatterplot(data = prinicpal_components,
#                     x = 'PC1',
#                     y = 'PC2',
#                     # hue = celldescriptors_clustered.loc[cell_IDs, 'hierarchical_cluster'].to_list(),
#                     # palette = colors[::1],
#                     hue = celldescriptors[parameter],
#                     # hue = MetaData.loc[cell_IDs, 'Region'],
#                     # palette = region_colors,
#                     ax = axs_PCA3[0])
    
#     sbn.move_legend(axs_PCA3[0], "upper left")
    
    
#     # plot cell_IDs
#     for cell_ID in prinicpal_components.index.to_list():
#         axs_PCA3[0].text(x = prinicpal_components.at[cell_ID, 'PC1'],
#                          y = prinicpal_components.at[cell_ID, 'PC2'],
#                          s = cell_ID,
#                          color = "k",
#                          ha = "center",
#                          va = "center",
#                          fontsize = 2)
    
    
    
#     # plot the variables as vectors
#     axs_PCA3[1].arrow(x = 0,
#                       y = 0,
#                       dx = celldescriptors_PCA.components_[0, parameter_i]*9,
#                       dy = celldescriptors_PCA.components_[1, parameter_i]*9,
#                       head_width=0.05,
#                       head_length=0.05,
#                       linewidth=0.75,
#                       color="w",
#                       )
        
    
#     # Plot annotations
#     axs_PCA3[1].text(x = (celldescriptors_PCA_eigen.loc[parameter, 'PC1'])*10,
#                      y = (celldescriptors_PCA_eigen.loc[parameter, 'PC2'])*10,
#                      s = parameter,
#                      color="red",
#                      fontsize = 6,
#                      )
    
    
#     PC_dict = {'ax_min' : -6,
#                'ax_max' : 8.5,
#                'pad' : None,
#                'step' : 2,
#                'stepminor' : 0.5,
#                'label' : None}
    
#     # edit axis
#     [apply_axis_settings(ax, axis = 'x', **PC_dict) for ax in axs_PCA3]
#     [apply_axis_settings(ax, axis = 'y', **PC_dict) for ax in axs_PCA3]
    
#     # remove spines
#     [ax.spines[spine].set_visible(False) for ax in axs_PCA3 for spine in ['top', 'right']]
    
#     # set directory for figure
#     fig_dir = join(clustering_dir, 'temp_figs')
    
#     # save figure
#     save_figures(fig_PCA3, f'figure-PCA-components_plot-{parameter}', 
#                  save_dir = fig_dir,
#                  darkmode_bool = darkmode_bool,
#                  figure_format = 'png')
    
#     # show plot
#     plt.show()




# %% scree plot

celldescriptors_PCA_dict = pd.DataFrame(columns = [f'PC{i+1}' for i in range(celldescriptors_zscored.shape[1])],
                                        index = celldescriptors_zscored.columns)

# get eigenvalues and calcuate the explained variance ratio 
celldescriptors_PCA_dict.loc["eigenvalue", :] = celldescriptors_PCA.explained_variance_
celldescriptors_PCA_dict.loc["explained_variance", :] = celldescriptors_PCA.explained_variance_ / celldescriptors_zscored.var().sum()

# calc cumulative sum of explained variance ratio
celldescriptors_PCA_dict.loc["cumsum_explained_variance", :] = celldescriptors_PCA_dict.loc["explained_variance", :].cumsum()



fig_scree, axs_scree = plt.subplots(nrows = 1,
                                   ncols = 2,
                                   layout = 'constrained',
                                   figsize = get_figure_size(height = 100, width = 200),
                                   dpi = 600)

axs_scree[0].plot(celldescriptors_PCA_dict.loc['explained_variance', :],
                  marker = 'x',
                  lw = 0.75)


axs_scree[1].plot(celldescriptors_PCA_dict.loc["cumsum_explained_variance", :],
                  marker = 'x',
                  lw = 0.75)


# edit axis
xscree_dict = {'ax_min' : 0,
               'ax_max' : celldescriptors_PCA_eigen.shape[1]-1,
               'pad' : None,
               'step' : 1,
               'stepminor' : 1,
               'label' : None,
               'ticklabels' : celldescriptors_PCA_eigen.columns.to_list(),
               'rotation' : 90}

[apply_axis_settings(ax, axis = 'x', **xscree_dict) for ax in axs_scree]


yscree_dict = {'ax_min' : 0,
               'ax_max' : 1,
               'pad' : None,
               'step' : 0.2,
               'stepminor' : 0.05,
               'label' : 'explained variance'}

# apply_axis_settings(ax_scree, axis = 'y', **yscree_dict)

axs_scree[0].set_ylabel('explained variance')
axs_scree[1].set_ylabel('cumultative explained variance')

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs_scree for spine in ['top', 'right']]

# set directory for figure
fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_scree, 'figure-PCA-scree_plot', 
             save_dir = fig_dir,
             darkmode_bool = darkmode_bool,
             figure_format = 'png')

# show plot
plt.show()


# %% 3d plot

# initialize figure
fig = plt.figure(layout = 'constrained',
                  dpi = 100)
ax = fig.add_subplot(projection='3d')

# plot each cluster separatly for coloring
for cluster_idx in range(n_clusters):
    
    # get cell_IDs of cells in cluster
    cluster_cellIDs = celldescriptors_clustered[celldescriptors_clustered['hierarchical_cluster'] == cluster_idx].index.to_list()

    # plot cluster
    ax.scatter(xs = prinicpal_components.loc[cluster_cellIDs, 'PC1'],
                ys = prinicpal_components.loc[cluster_cellIDs, 'PC2'],
                zs = prinicpal_components.loc[cluster_cellIDs, 'PC3'],
                c = [colors[cluster_idx]] * len(cluster_cellIDs))

PC_dict = {'ax_min' : -6,
            'ax_max' : 8.5,
            'pad' : None,
            'step' : 2,
            'stepminor' : 0.5,
            'label' : None}

# edit axis
apply_axis_settings(ax, axis = 'x', **PC_dict)
apply_axis_settings(ax, axis = 'y', **PC_dict)
apply_axis_settings(ax, axis = 'z', **PC_dict)

# set titles
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')

# set background colors
for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
    axis.set_pane_color((1.0, 1.0, 1.0), 0.0) 
    axis._axinfo['grid']['linewidth'] = 0.0
    axis._axinfo['grid']['color'] = "#d1d1d1"
    axis._axinfo['tick']['inward_factor'] = 0.2
    axis._axinfo['tick']['outward_factor'] = 0.0

# remove grid
ax.grid(False)

# set initial viewing angle
ax.view_init(elev = 90, azim = -90, roll = 0)

# display figure
plt.show()



# %% PCA 3d as subplots
fig_3d_subp, axs_3ds = plt.subplots(nrows = 2,
                                    ncols = 4,
                                    layout = 'constrained',
                                    figsize = get_figure_size(height = 120, width = 240),
                                    dpi = 600)

# flatten axis array
axs_3ds = axs_3ds.flatten()

PC_dict = {'PC1': {'ax_min' : -6,
                   'ax_max' : 9,
                   'pad' : None,
                   'step' : 2,
                   'stepminor' : 1,
                   'label' : None},
           'PC2': {'ax_min' : -6,
                   'ax_max' : 7,
                   'pad' : None,
                   'step' : 2,
                   'stepminor' : 1,
                   'label' : None},
           'PC3': {'ax_min' : -4,
                   'ax_max' : 6,
                   'pad' : None,
                   'step' : 2,
                   'stepminor' : 1,
                   'label' : None}}


for ax_i, PCs in zip([0, 1, 4], [('PC1', 'PC2'), ('PC3', 'PC2'), ('PC1', 'PC3')]):
    
    # plot each cluster separatly for coloring
    for cluster_idx in range(n_clusters):
        
        # set axis
        ax = axs_3ds[ax_i]
        
        # get cell_IDs of cells in cluster
        cluster_cellIDs = celldescriptors_clustered[celldescriptors_clustered['hierarchical_cluster'] == cluster_idx].index.to_list()

        # if cluster_idx in [2, 3]:

        # plot cluster
        ax.scatter(x = prinicpal_components.loc[cluster_cellIDs, PCs[0]],
                   y = prinicpal_components.loc[cluster_cellIDs, PCs[1]],
                   c = [colors[cluster_idx]] * len(cluster_cellIDs),
                   s = 5,
                   label = cluster_idx)
    
        # set axis labels
        ax.set_xlabel(PCs[0])
        ax.set_ylabel(PCs[1])
        
        # set axis
        apply_axis_settings(ax = ax, axis = 'x', **PC_dict[PCs[0]])
        apply_axis_settings(ax = ax, axis = 'y', **PC_dict[PCs[1]])
 
        
for ax_i, PCs in zip([2, 3, 6], [('PC1', 'PC2'), ('PC3', 'PC2'), ('PC1', 'PC3')]):
    
    # plot each cluster separatly for coloring
    for cluster_idx in range(n_clusters):
        
        # set axis
        ax = axs_3ds[ax_i]
        
        # plot the variables as vectors
        for i in range(celldescriptors_PCA.components_.shape[1]):
            ax.arrow(x = 0,
                     y = 0,
                     dx = celldescriptors_PCA.components_[int(PCs[0][-1])-1, i]*4.5,
                     dy = celldescriptors_PCA.components_[int(PCs[1][-1])-1, i]*4.5,
                     head_width=0.05,
                     head_length=0.05,
                     linewidth=0.75,
                     color="w",
                     )
   
        
        # Plot annotations
        for parameter in celldescriptors_PCA_eigen.index.to_list():
            ax.text(x = (celldescriptors_PCA_eigen.loc[parameter, PCs[0]])*5,
                    y = (celldescriptors_PCA_eigen.loc[parameter, PCs[1]])*5,
                    s = parameter,
                    color="red",
                    fontsize = 4,
                    )
   
        # set axis labels
        ax.set_xlabel(PCs[0])
        ax.set_ylabel(PCs[1])
        
        # set axis
        apply_axis_settings(ax = ax, axis = 'x', **PC_dict[PCs[0]])
        apply_axis_settings(ax = ax, axis = 'y', **PC_dict[PCs[1]])       
       
        
# legend
h, l = axs_3ds[4].get_legend_handles_labels() 
axs_3ds[5].legend(h, l, title = 'cluster id', 
                  frameon = False, ncol = 2, loc = 'center left')

# remove spines, ticks, and ticklabels
remove_spines_n_ticks([axs_3ds[5]], axis = 'x')
remove_spines_n_ticks([axs_3ds[5]], axis = 'y')
axs_3ds[5].set_xticks([])
axs_3ds[5].set_yticks([])

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs_3ds for spine in ['top', 'right']]

# remove unused subplots
fig_3d_subp.delaxes(axs_3ds[7])

# show plot
plt.show()


# %% loadings plot


eigenvectors = celldescriptors_PCA.components_

celldescriptors_PCA_eigenvectors = pd.DataFrame(eigenvectors,
                                                columns = celldescriptors.columns,
                                                index = [f'PC{i+1}' for i in range(celldescriptors.shape[1])])

loadings = eigenvectors.T * np.sqrt(celldescriptors_PCA_explained_variance)

celldescriptors_PCA_loadings = pd.DataFrame(loadings,
                                            columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])],
                                            index = celldescriptors.columns)


plt.barh(y = celldescriptors.columns,
         width = celldescriptors_PCA_loadings['PC1'].to_list(),
         left = 0,
         zorder = 1)

plt.barh(y = celldescriptors.columns,
         width = celldescriptors_PCA_loadings['PC2'].to_list(),
         left = 2,
         zorder = 1)


axline_dict = {'linewidth' : 1.5, 
               'color' : 'w',
               'zorder' : 0,
               'linestyle' : 'solid',
               'alpha' : 0.75}

plt.axline((0, 1), (0, 2), **axline_dict)
plt.axline((2, 1), (2, 2), **axline_dict)

# remove spines
[plt.gca().spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom', 'left']]

plt.xticks(ticks = [0, 2],
           labels = ['PC1', 'PC2'])

# invert y axis
plt.gca().invert_yaxis()

plt.xlim([-1, 3])
