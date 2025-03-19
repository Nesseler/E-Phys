# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:00:13 2025

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

# calc distance between clusters
ward_clustering_linkage = linkage(cellmetrics_zscored, method="ward", metric="euclidean")

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
fig_elbow.suptitle(f'hierarchical clustering - elbow plot',
                    fontsize = 12)

# plot cluster distances
ax_elbow.plot(last_clusters_idc, last_clusters_rev,
              lw = 1,
              c = colors_dict['primecolor'],
              label = 'Cluster distance')

# plot acceleration
ax_elbow.plot(last_clusters_idc[:-2] + 1, acceleration_rev,
              lw = 1,
              c = 'r',
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
save_figures(fig_elbow, f'hierarchical_clustering-elbow_plot', 
             save_dir = join(cellmorph_analysis_dir, 'plots-hierarchical_clustering'),
             darkmode_bool = darkmode_bool,
             figure_format = 'png')
