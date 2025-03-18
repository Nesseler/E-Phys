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


# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()


# %% principal component analysis

# set up PCA
celldescriptors_PCA = PCA().fit(celldescriptors_zscored)

# perform transfrom
celldescriptors_PCA_components = celldescriptors_PCA.transform(celldescriptors_zscored)

# write to dataframe
prinicpal_components = pd.DataFrame(celldescriptors_PCA_components,
                                    index = cell_IDs,
                                    columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])])

# get eigenvectors
eigenvectors = celldescriptors_PCA.components_

celldescriptors_PCA_eigenvectors = pd.DataFrame(eigenvectors,
                                                columns = celldescriptors.columns,
                                                index = [f'PC{i+1}' for i in range(celldescriptors.shape[1])])

# get the eigenvalues
celldescriptors_PCA_eigenvalues = celldescriptors_PCA.explained_variance_

# get components loading
loadings = eigenvectors.T * np.sqrt(celldescriptors_PCA_eigenvalues)

celldescriptors_PCA_loadings = pd.DataFrame(loadings,
                                            columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])],
                                            index = celldescriptors.columns)

# combine to one dataframe 
celldescriptors_PCA_metrics = pd.DataFrame(index = [f'PC{i+1}' for i in range(celldescriptors_zscored.shape[1])])

# get eigenvalues and calcuate the explained variance ratio 
celldescriptors_PCA_metrics.loc[:, "eigenvalue"] = celldescriptors_PCA_eigenvalues
celldescriptors_PCA_metrics.loc[:, "explained_variance"] = celldescriptors_PCA_eigenvalues / celldescriptors_zscored.var().sum()

# calc cumulative sum of explained variance ratio
celldescriptors_PCA_metrics.loc[:, "cumsum_explained_variance"] = celldescriptors_PCA_metrics.loc[:, "explained_variance"].cumsum()




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




# %% initialize plotting

from functions.initialize_plotting import *

# %% pca and vectors

fig_PCA2, axs_PCA2 = plt.subplots(nrows = 1,
                                 ncols = 2,
                                 layout = 'constrained',
                                 figsize = get_figure_size(width = 220, height = 120),
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
                        s = 30,
                        label = cluster_idx)

# # plot cell_IDs
# for cell_ID in prinicpal_components.index.to_list():
#     axs_PCA2[0].text(x = prinicpal_components.at[cell_ID, 'PC1'],
#                      y = prinicpal_components.at[cell_ID, 'PC2'],
#                      s = cell_ID,
#                      color = "k",
#                      ha = "center",
#                      va = "center",
#                      fontsize = 3)

# legend
h, l = axs_PCA2[0].get_legend_handles_labels() 
axs_PCA2[0].legend(h, l, 
                   title = 'Cluster ID',
                   title_fontsize = 12,
                   frameon = False, 
                   ncol = 2, 
                   loc = 'upper left',
                   fontsize = 12,
                   handletextpad = 0.05,
                   columnspacing = 0.4)

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
for parameter in celldescriptors_PCA_eigenvectors.columns.to_list():
    axs_PCA2[1].text(x = (celldescriptors_PCA_eigenvectors.loc['PC1', parameter])*10,
                     y = (celldescriptors_PCA_eigenvectors.loc['PC2', parameter])*10,
                     s = parameter,
                     color="w",
                     fontsize = 9,
                     ha = 'center',
                     wrap = True)

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
    ax.set_aspect(1)
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
                                              figsize = get_figure_size(width = 220, height = 120),
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
                              s = 30,
                              label = cluster_idx)

# plot cell_IDs
# for cell_ID in prinicpal_components.index.to_list():
#     axs_PCA_region[0].text(x = prinicpal_components.at[cell_ID, 'PC1'],
#                      y = prinicpal_components.at[cell_ID, 'PC2'],
#                      s = cell_ID,
#                      color = "k",
#                      ha = "center",
#                      va = "center",
#                      fontsize = 3)

# legend
h, l = axs_PCA_region[0].get_legend_handles_labels() 
axs_PCA_region[0].legend(h, ['BAOT', 'BAOT/MeA', 'MeA'], 
                         title = 'Region', 
                         title_fontsize = 12,
                         frameon = False, 
                         ncol = 1, 
                         loc = 'upper left',
                         fontsize = 12,
                         handletextpad = 0.05,
                         columnspacing = 0.4)


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
for parameter in celldescriptors_PCA_eigenvectors.columns.to_list():
    axs_PCA_region[1].text(x = (celldescriptors_PCA_eigenvectors.loc['PC1', parameter])*10,
                           y = (celldescriptors_PCA_eigenvectors.loc['PC2', parameter])*10,
                           s = parameter,
                           color="w",
                           fontsize = 9,
                           ha = 'center',
                           wrap = True)


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
    ax.set_aspect(1)
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
parameter = 'rheobase_rel'
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
ax.text(x = (celldescriptors_PCA_eigenvectors.loc['PC1', parameter])*10,
        y = (celldescriptors_PCA_eigenvectors.loc['PC2', parameter])*10,
        s = parameter,
        color="w",
        fontsize = 9,
        ha = 'center',
        wrap = True)


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

# %% all single parameter


# # %% single parameter color code

# parameter = 'r_input'

# for parameter in celldescriptors_PCA_eigenvectors.columns.to_list():

#     parameter_i = celldescriptors_PCA_eigenvectors.columns.to_list().index(parameter)
    
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
        
    
# # Plot annotations
#     axs_PCA3[1].text(x = (celldescriptors_PCA_eigenvectors.loc['PC1', parameter])*10,
#                      y = (celldescriptors_PCA_eigenvectors.loc['PC2', parameter])*10,
#                      s = parameter,
#                      color="w",
#                      fontsize = 9,
#                      ha = 'center',
#                      wrap = True)
    
    
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

# set font size
mtl.rcParams.update({'font.size': 9, 'font.family' : 'Arial'})

fig_scree, axs_scree = plt.subplots(nrows = 1,
                                    ncols = 1,
                                    layout = 'constrained',
                                    figsize = get_figure_size(height = 100, width = 120),
                                    dpi = 300)

axs_scree.plot(celldescriptors_PCA_metrics.loc[:, 'explained_variance'],
                  marker = 'x',
                  markersize = 3,
                  markeredgewidth = 0.5,
                  lw = 0.75,
                  alpha = 0.9,
                  color = colors_dict['primecolor'])

ax_cumsum = axs_scree.twinx()

ax_cumsum.plot(celldescriptors_PCA_metrics.loc[:, "cumsum_explained_variance"],
                  marker = 'x',
                  markersize = 3,
                  markeredgewidth = 0.5,
                  lw = 0.75,
                  alpha = 0.9,
                  color = 'r')

ax_cumsum.spines['right'].set_color('r')
ax_cumsum.tick_params(axis='y', color = 'r', which = 'both')
ax_cumsum.yaxis.label.set_color('w')

# edit axis
xscree_dict = {'ax_min' : 0,
               'ax_max' : celldescriptors_PCA_eigenvectors.shape[0]-1,
               'pad' : None,
               'step' : 1,
               'stepminor' : 1,
               'label' : None,
               'ticklabels' : celldescriptors_PCA_eigenvectors.index.to_list(),
               'rotation' : 90}

apply_axis_settings(axs_scree, axis = 'x', **xscree_dict)
apply_axis_settings(ax_cumsum, axis = 'x', **xscree_dict)

yscree_dict = {'ax_min' : 0,
                'ax_max' : 0.4,
                'pad' : None,
                'step' : 0.1,
                'stepminor' : 0.05,
                'label' : 'explained variance'}

apply_axis_settings(axs_scree, axis = 'y', **yscree_dict)
axs_scree.set_ylabel('explained variance')

ycumsum_dict = {'ax_min' : 0.2,
                'ax_max' : 1,
                'pad' : None,
                'step' : 0.1,
                'stepminor' : 0.05,
                'label' : ''}


apply_axis_settings(ax_cumsum, axis = 'y', **ycumsum_dict)
ax_cumsum.spines['right'].set_bounds([ycumsum_dict['ax_min'], ycumsum_dict['ax_max']])
ax_cumsum.set_ylabel('cumultative explained variance')

# remove spines
[axs_scree.spines[spines].set_visible(False) for spines in ['top', 'right']]
[ax_cumsum.spines[spines].set_visible(False) for spines in ['top', 'left']]

# set directory for figure
fig_dir = join(clustering_dir, 'temp_figs')

# save figure
save_figures(fig_scree, 'figure-PCA-scree_plot', 
              save_dir = fig_dir,
              darkmode_bool = darkmode_bool,
              figure_format = 'both')

# show plot
plt.show()

# set font size
mtl.rcParams.update({'font.size': 12, 'font.family' : 'Arial'})

# %% 3d plot

# %matplotlib qt5

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


# %matplotlib inline

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
    for cluster_idx in [4, 6, 7]: #range(n_clusters):
        
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
        for parameter in celldescriptors_PCA_eigenvectors.columns.to_list():
            ax.text(x = (celldescriptors_PCA_eigenvectors.loc[PCs[0], parameter])*10,
                    y = (celldescriptors_PCA_eigenvectors.loc[PCs[1], parameter])*10,
                    s = parameter,
                    color="w",
                    fontsize = 6,
                    ha = 'center',
                    wrap = True)
        
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

# set font size
mtl.rcParams.update({'font.size': 9, 'font.family' : 'Arial'})

axline_dict = {'linewidth' : 1.5, 
               'color' : 'w',
               'zorder' : 0,
               'linestyle' : 'solid',
               'alpha' : 0.75}


fig, axs= plt.subplots(nrows = 1,
                       ncols = 3,
                       layout = 'constrained',
                       figsize = get_figure_size(height = 100, width = 150),
                       dpi = 300,
                       sharey = True)

for ax, PC, color in zip(axs, ['PC1', 'PC2', 'PC3'], ['mediumturquoise', 'sandybrown', 'mediumorchid']):

    # plot PC loadings as bar graph
    ax.barh(y = celldescriptors_PCA_loadings.index,
            width = celldescriptors_PCA_loadings.loc[:, PC].to_list(),
            left = 0,
            zorder = 1,
            color = color)

    # plot zero line
    ax.axline((0, 1), (0, 2), **axline_dict)
    
    # set x axis
    ax.set_xticks(ticks = [0], labels = [PC])
    ax.set_xticks(ticks = np.arange(-1, 1+0.1, 0.5), labels = [], minor = True)
    ax.set_xlim([-1.15, 1.15])
    
    # set y axis
    ax.set_ylim([0-0.5, 20-0.5])
    ax.invert_yaxis()

    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom', 'left']]

    # set grid
    ax.grid(axis = 'x', which = 'both', lw = 0.5, alpha = 0.5)
    ax.set_axisbelow(True)
    
    
# remove y ticks
remove_spines_n_ticks(axs[1:], axis = 'y')

# display figure
plt.show()


# %% sorted loadings

fig, axs= plt.subplots(nrows = 3,
                       ncols = 1,
                       layout = 'constrained',
                       figsize = get_figure_size(height = 230, width = 100),
                       dpi = 300,
                       sharey = False)

for ax, PC, color in zip(axs, ['PC1', 'PC2', 'PC3'], ['mediumturquoise', 'sandybrown', 'mediumorchid']):

    # sort     
    PC_loadings = celldescriptors_PCA_loadings.sort_values(by = PC, key = abs, ascending=False)
    labels = PC_loadings.index

    # plot PC loadings as bar graph
    ax.barh(y = labels,
            width = PC_loadings.loc[:, PC].to_list(),
            left = 0,
            zorder = 1,
            color = color)

    # plot zero line
    ax.axline((0, 1), (0, 2), **axline_dict)
    
    # set x axis
    ax.set_xticks(ticks = [0], labels = [PC])
    ax.set_xticks(ticks = np.arange(-1, 1+0.1, 0.5), labels = [], minor = True)
    ax.set_xlim([-1.15, 1.15])
    
    # set y axis
    ax.set_ylim([0-0.5, 20-0.5])
    ax.invert_yaxis()

    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom', 'left']]

    # set grid
    ax.grid(axis = 'x', which = 'both', lw = 0.5, alpha = 0.5)
    ax.set_axisbelow(True)
    
    
# remove y ticks
remove_spines_n_ticks(axs[1:], axis = 'y')

# display figure
plt.show()



# %% length of eigenvectors

combined_first_twoPCs_loadings_sorted = celldescriptors_PCA_loadings.loc[:, ['PC1', 'PC2']].abs().sum(axis = 1).sort_values(ascending = False)



# vector_length = np.sqrt((celldescriptors_PCA_loadings**2).loc[:, ['PC1', 'PC2']].sum(axis = 1))
vector_length = np.sqrt((celldescriptors_PCA_eigenvectors.T**2).loc[:, ['PC1', 'PC2']].sum(axis = 1))
vector_length_sorted = vector_length.sort_values(ascending = False)


fig, ax = plt.subplots(nrows = 1,
                       ncols = 1,
                       layout = 'constrained',
                       figsize = get_figure_size(height = 100, width = 100),
                       dpi = 300,
                       sharey = False)


# plot PC loadings as bar graph
ax.barh(y = vector_length_sorted.index,
        width = vector_length_sorted.to_list(),
        left = 0,
        zorder = 1,
        color = 'cornflowerblue')

# plot zero line
ax.axline((0, 1), (0, 2), **axline_dict)

# set x axis
ax.set_xticks(ticks = [0.5], labels = [r'$\sqrt{PC1^2 + PC2^2}$'])
ax.set_xticks(ticks = np.arange(-1, 1+0.1, 0.5), labels = [], minor = True)
ax.set_xlim([-0.2, 1.2])

# set y axis
ax.set_ylim([0-0.5, 20-0.5])
ax.invert_yaxis()

# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right', 'bottom', 'left']]

# set grid
ax.grid(axis = 'x', which = 'both', lw = 0.5, alpha = 0.5)
ax.set_axisbelow(True)



# %% PCA region joint histogram

bin_width = 1
bin_edges = np.arange(-6, 11, bin_width)

PC_histogram = pd.DataFrame(index = bin_edges)

for PC in ['PC1', 'PC2']:
    for region in ['BAOT/MeA', 'MeA', 'BAOT']:
        
        # get cell_IDs per region
        region_cellIDs = MetaData[MetaData['Region'] == region].index.to_list()
        
        # get histogram
        counts, bins = np.histogram(prinicpal_components.loc[region_cellIDs, PC], 
                                    bins = bin_edges)
        
        # write to dataframe
        PC_histogram.loc[bins[:-1], f'{region}-{PC}'] = counts



# set font size
mtl.rcParams.update({'font.size': 12, 'font.family' : 'Arial'})

fig, axs = plt.subplots(nrows = 2,
                        ncols = 2,
                        layout = 'constrained',
                        figsize = get_figure_size(width = 160.5),
                        dpi = 300,
                        height_ratios = [0.8, 0.2],
                        width_ratios = [0.2, 0.8],
                        sharex = 'col',
                        sharey = 'row')

# flatten axis
axs = axs.flatten()

# set PCA axis
ax = axs[1]

# axis lines
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
    ax.scatter(x = prinicpal_components.loc[cluster_cellIDs, 'PC1'],
               y = prinicpal_components.loc[cluster_cellIDs, 'PC2'],
               c = [colors[cluster_idx]] * len(cluster_cellIDs),
               s = 60,
               label = cluster_idx)

# plot cell_IDs
for cell_ID in prinicpal_components.index.to_list():
    ax.text(x = prinicpal_components.at[cell_ID, 'PC1'],
                     y = prinicpal_components.at[cell_ID, 'PC2'],
                     s = cell_ID,
                     color = "k",
                     ha = "center",
                     va = "center",
                     fontsize = 3)

# legend
h, l = ax.get_legend_handles_labels() 
ax.legend(h, l, 
          title = 'cluster id', 
          frameon = False, 
          ncol = 2, 
          loc = 'upper left',
          fontsize = 9)


PC_dict = {'ax_min' : -6,
           'ax_max' : 8.5,
           'pad' : None,
           'step' : 2,
           'stepminor' : 0.5,
           'label' : None}

# edit axis
apply_axis_settings(ax, axis = 'x', **PC_dict)
apply_axis_settings(ax, axis = 'y', **PC_dict)

# axis labels
axs[3].set_xlabel('PC 1')
axs[0].set_ylabel('PC 2')

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]


# set axis to histogram
ax = axs[3]

# set array for stacking
bottom = np.zeros_like(bin_edges)

# plot stacked bar plot
for ri, region in enumerate(['BAOT/MeA', 'MeA', 'BAOT']):
    
    # get counts
    region_counts = PC_histogram[f'{region}-PC1']
    
    p = ax.bar(PC_histogram.index, 
               region_counts,
               align = 'edge', 
               width = bin_width, 
               label = region, 
               bottom = bottom,
                edgecolor = 'None',
                facecolor = region_colors[region],
                alpha = 0.5,
                zorder = 3-ri,
                lw = 1.5)
    
    p = ax.bar(PC_histogram.index, 
               region_counts,
               align = 'edge', 
               width = bin_width, 
               label = region, 
               bottom = bottom,
                edgecolor = region_colors[region],
                facecolor = 'None',
                zorder = 3-ri,
                lw = 1.5)
    
    # add to stacking array
    bottom += region_counts


# set axis to histogram
ax = axs[0]

# set array for stacking
bottom = np.zeros_like(bin_edges)

# plot stacked bar plot
for ri, region in enumerate(['BAOT/MeA', 'MeA', 'BAOT']):
    
    # get counts
    region_counts = PC_histogram[f'{region}-PC2']
    
    p = ax.barh(PC_histogram.index, 
                region_counts, 
                height = bin_width,
                align = 'edge',
                label = region,
                left = bottom,
                edgecolor = 'None',
                facecolor = region_colors[region],
                alpha = 0.5,
                zorder = 3-ri,
                lw = 1.5)
    
    p = ax.barh(PC_histogram.index, 
                region_counts, 
                height = bin_width,
                align = 'edge',
                label = region,
                left = bottom,
                edgecolor = region_colors[region],
                facecolor = 'None',
                zorder = 3-ri,
                lw = 1.5)
    
    # add to stacking array
    bottom += region_counts

# delete unused axis
fig.delaxes(axs[2])

# # set directory for figure
# fig_dir = join(clustering_dir, 'temp_figs')

# # save figure
# save_figures(fig_PCA_region, 'figure-PCA-components_plot-region', 
#              save_dir = fig_dir,
#              darkmode_bool = darkmode_bool,
#              figure_format = 'both')

# show plot
plt.show()

