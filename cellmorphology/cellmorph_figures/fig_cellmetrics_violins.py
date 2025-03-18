# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 11:19:51 2025

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# init plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *

# get MetaData
from functions.functions_import import get_MetaData
MetaData = get_MetaData()

# define regions
regions = ['MeA', 'BAOT']

# set neurite types
neurite_types = ['dendrites', 'axons']


# %% load data

from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_analysis_dir, cellmorph_metrics_dir, cellmorph_figures_dir

# set dataframe
cellmetrics = pd.DataFrame()

# set filename
metrics = ['total_cable_length', 'n_primary', 'n_terminal', 'bifurcation_ratio']

# iterate through files
for filename in metrics:
    loaded_file = pd.read_excel(join(cellmorph_metrics_dir, filename + '.xlsx'),
                                index_col = 'cell_ID')
    
    # concat
    cellmetrics = pd.concat([cellmetrics, loaded_file], axis = 1)
    
# get cell_IDs
cell_IDs = cellmetrics.index.to_list()[:-20]

# limit dataframe
cellmetrics = cellmetrics.loc[cell_IDs, :]


# %% cell_IDs dicts

# get all cell_IDs
cell_IDs = cellmetrics.index.to_list()

# get cell_IDs
cell_IDs_dict = {'all': cell_IDs,
                 'MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in cell_IDs],
                 'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in cell_IDs]}


# %% 

import matplotlib.transforms as transforms

# initialise figure
fig, axs = plt.subplots(nrows=1,
                        ncols=4,
                        figsize=get_figure_size(width=260.334, height=90),
                        layout='constrained',
                        dpi=600)

# flatten axes array
axs = axs.flatten()

for m_idx, metric in enumerate(metrics):
    
    # set axis
    ax = axs[m_idx]
    
    for n_idx, ntype in enumerate(['dendrites', 'axons']):
        
        for (r_idx, region), v_direction in zip(enumerate(['MeA', 'BAOT']), [-1, 1]):
            
            # get region cell_Id
            region_cellIDs = cell_IDs_dict[region]
            
            # get data
            metrics_toplot = cellmetrics.loc[region_cellIDs, f'{metric}-{ntype}'].dropna()
            
            # calc mean, etc.
            metrics_mean = metrics_toplot.mean()
            metrics_median = metrics_toplot.median()
            metrics_std = metrics_toplot.std()
            
            # set swarm x
            swarm_x = n_idx*2 + r_idx
            
            # plot swarmplot
            sbn.swarmplot(x = [swarm_x] * len(metrics_toplot),
                          y = metrics_toplot,
                          color = colors_dict['primecolor'],
                          ax = ax,
                          size = 2,
                          zorder = 1)
            
            # calc violin position
            x = swarm_x - (v_direction*0.3)
            
            # plot half violin
            plot_half_violin(data = metrics_toplot, 
                              ax = ax,
                              v_position = x,
                              v_direction = v_direction,
                              v_offset = 0,
                              v_lw = 1.5,
                              v_color = neurite_color_dict[region][ntype],
                              v_zorder = 2,
                              v_width = 0.8,
                              v_abs_cutoff = [0, np.nan])
            
            # errorbar
            ax.errorbar(x = x,
                        y = metrics_mean,
                        yerr = metrics_std,
                        fmt='_', 
                        markersize = 6,
                        markerfacecolor = 'none',
                        capsize = 2,
                        color = neurite_color_dict[region][ntype],
                        linewidth = 1,
                        label = '_nolegend_',
                        zorder = 3)
            
            # plot median
            ax.scatter(x = x,
                       y = metrics_median,
                       marker='D', 
                       s = 5,
                       color = neurite_color_dict[region][ntype],
                       linewidth = 1,
                       label = '_nolegend_',
                       zorder = 4)
            
    # edit main plot axes
    # y
    ydict = {'total_cable_length': {'ax_min' : 0, 'ax_max' : 5500, 'pad' : None, 'step' : 2000,'stepminor' : 1000, 'label' : 'Total cable length [µm]'        , 'pad_factor' : 0.05},
             'n_primary'         : {'ax_min' : 0, 'ax_max' : 10,   'pad' : None, 'step' : 2,   'stepminor' : 1,    'label' : 'Number of\nprimary branches [#]' , 'pad_factor' : 0.05},
             'n_terminal'        : {'ax_min' : 0, 'ax_max' : 40,   'pad' : None, 'step' : 10,  'stepminor' : 5,    'label' : 'Number of\nterminal branches [#]', 'pad_factor' : 0.05},
             'bifurcation_ratio' : {'ax_min' : 0, 'ax_max' : 17,   'pad' : None, 'step' : 5,   'stepminor' : 1,    'label' : 'Bifurcation ratio'               , 'pad_factor' : 0.05}}
    
    apply_axis_settings(ax, axis = 'y', **ydict[metric])
    
    # x
    xdict = {'ax_min' : 0,
             'ax_max' : 3,
             'pad' : 0.8,
             'step' : 1,
             'stepminor' : 1,
             'label' : ''}
    
    apply_axis_settings(ax, axis = 'x', **xdict)
    

    ax.set_xticks(ticks = [0, 1, 2, 3], 
                  labels = regions + regions,
                  rotation = 90)
    
    
    sec = ax.secondary_xaxis(location=-0.32)
    sec.set_xticks(ticks = [0.5, 2.5], labels = [n.capitalize() for n in neurite_types])
    
    xdict2 = {'ax_min' : 0.5,
              'ax_max' : 2.5,
              'pad' : 1.8,
              'step' : 2,
              'stepminor' : 2,
              'label' : ''}
    
    apply_axis_settings(sec, axis = 'x', **xdict2)
    
    # remove top and right spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    

            
# align labels
fig.align_labels()

# show figure
plt.show()

# save figure
save_figures(fig, 
              figure_name = 'cell_metrics', 
              save_dir = join(cellmorph_figures_dir, 'cell_metrics'),
              darkmode_bool=darkmode_bool, 
              figure_format='both')