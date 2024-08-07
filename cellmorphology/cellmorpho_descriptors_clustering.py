# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:22:45 2024

@author: nesseler
"""

from cellmorphology.cellmorph_packages import mtl, plt, sbn, pd, np, join

from parameters.directories_win import table_file, cell_morph_descrip_dir

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# get cell_IDs:
cell_IDs = MetaData[MetaData['reconstructed'] == 1].index.to_list()

# set list of neurite types
neurite_types = ['neurites', 'dendrites', 'axons']


# %%

# create dataframe that contains all parameters
cellmorph_descriptors = pd.DataFrame(index = cell_IDs)
cellmorph_descriptors.index.name = 'cell_ID'


# %% load height & width

# load
height_width = pd.read_excel(join(cell_morph_descrip_dir, 'height_width.xlsx'), index_col = 'cell_IDs')

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


# %% remove cells that do not contain all analysed values

cellmorph_descriptors.drop(index = ['E-126', 'E-158'], inplace = True)

# %% replace nan values in cell descriptor table

cellmorph_descriptors.fillna(value = 0, inplace = True)

# %% normalise cell descriptors

# min-max normalize cellmorph matrix
cellmorph_descriptors_minmax = (cellmorph_descriptors - cellmorph_descriptors.min()) / (cellmorph_descriptors.max() - cellmorph_descriptors.min())

# z-score cellmorph matrix
cellmorph_descriptors_zscored = (cellmorph_descriptors - cellmorph_descriptors.mean()) / cellmorph_descriptors.std()


# %% sort dataframe

# cellmorph_descriptors_zscored.sort_values(['total_cable_length-axons'], inplace = True)


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
    
    
# %% plot z-scored distribution

fig_dist_zscored, ax_dist_zscored = create_data_distribution_figure(cellmorph_descriptors_zscored)    

# set axis label
ax_dist_zscored.set_ylabel('Z-scored parameter value [std]')

# show plot
plt.show()


# %% plot min max distribution

fig_dist_minmax, ax_dist_minmax = create_data_distribution_figure(cellmorph_descriptors_minmax)    

# set axis label
ax_dist_minmax.set_ylabel('Min-Max-Normalised parameter value')

# edit axis
ymin = -0.5
ymax = 1.5
ypad = 0.15
ystep = 1
ystepminor = 0.25
    
ax_dist_minmax.set_ylim([ymin - ypad, ymax + ypad])
ax_dist_minmax.set_yticks(ticks = np.arange(0, ymax+ ystepminor, ystep))
ax_dist_minmax.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystepminor), minor = True)
ax_dist_minmax.spines['left'].set_bounds([ymin, ymax])

# show plot
plt.show()
    

# %% plot min max distribution

fig_dist, ax_dist = create_data_distribution_figure(cellmorph_descriptors)    

# set axis label
ax_dist.set_ylabel('Parameter value [µm / # / ]')

# edit axis
ymin = 0
ymax = 6000
ypad = 100
ystep = 2000
ystepminor = 500
    
ax_dist.set_ylim([ymin - ypad, ymax + ypad])
ax_dist.set_yticks(ticks = np.arange(0, ymax+ ystepminor, ystep))
ax_dist.set_yticks(ticks = np.arange(ymin, ymax+ ystepminor, ystepminor), minor = True)
ax_dist.spines['left'].set_bounds([ymin, ymax])

# show plot
plt.show()    


# %% create figure of all plots

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
    
# show plot
plt.show


# %% correlation analysis


# initialize figure
fig_corr_heat, axs_corr_heat = plt.subplots(nrows = 2,
                                            ncols = 1,
                                            layout = 'constrained',
                                            figsize = get_figure_size(width = 150, height = 225),
                                            dpi = 600,
                                            sharey = True,
                                            sharex = True)

# set string for color map
cmap_str = 'seismic'

# get correlation matrix of dataframe
cellmorph_descriptors_corr = cellmorph_descriptors.corr()

# plot heatmap of correlation values
sbn.heatmap(data= cellmorph_descriptors_corr,
            cmap = cmap_str,
            vmin = -1,
            vmax = 1,
            annot = False,
            cbar_kws={'label': 'Correlation coefficient', 'ticks' : [-1, 0, 1]},
            ax = axs_corr_heat[0])

# set axis title
axs_corr_heat[0].set_title('A: Pearson correlation coefficient',
                           fontsize=12, 
                           loc='left',
                           x = -0.62)


# plot heatmap with values above threshold only
corr_threshold = 0.8

# replace values in dataframe below threshold with nan
cellmorph_descriptors_corr_thresh = cellmorph_descriptors_corr[(cellmorph_descriptors_corr > corr_threshold) | (cellmorph_descriptors_corr < -corr_threshold)]

# plot heatmap of correlation values
sbn.heatmap(data= cellmorph_descriptors_corr_thresh,
            cmap = cmap_str,
            vmin = -1,
            vmax = 1,
            annot = False,
            cbar_kws={'label': 'Correlation coefficient', 'ticks' : [-1, -0.8, 0.8, 1]},
            ax = axs_corr_heat[1],
            xticklabels = 1)

axs_corr_heat[1].set_title(r'B: Pearson correlation coefficient $\pm$' + str(corr_threshold),
                           fontsize=12, 
                           loc='left',
                           x = -0.62)

plt.show()


# %%
# sbn.pairplot(cellmorph_descriptors, kind="reg")
# plt.show()






# %% heatmap


fig_heat, ax_heat = plt.subplots(nrows = 1,
                                  ncols = 1,
                                  layout = 'constrained',
                                  figsize = get_figure_size(width = 150, height = 125),
                                  dpi = 600)


sbn.heatmap(cellmorph_descriptors_zscored,
            vmin = -2,
            vmax = 2,
            square = False, 
            ax = ax_heat, 
            cmap="flare_r", 
            yticklabels=False,
            linewidth = 0) 

    
# show plot
plt.show()