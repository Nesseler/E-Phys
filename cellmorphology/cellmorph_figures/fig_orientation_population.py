# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 16:34:00 2025

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
regions = ['BAOT', 'MeA']

# set neurite types
neurite_types = ['dendrites', 'axons']


# %% load data

from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_metrics_dir

population_orientation = pd.read_excel(join(cellmorph_metrics_dir, 'population_orientation.xlsx'), 
                                       index_col = 'ntype-region-measure')


# %% polar plot specific

from cellmorphology.cellmorph_functions.cellmorph_polarplot_functions import n_bins, hist_bins, binsize, orientation_labels


# %% figure

# inititalise figure
fig, axs = plt.subplots(nrows=2,
                        ncols=4,
                        figsize=get_figure_size(width=260.334, height=145.54),
                        layout='constrained',
                        dpi=600,
                        )

# flatten axes array
axs = axs.flatten()

# set titles
axis_titles = {0: 'BAOT',
               2: 'MeA',
               4: 'Dendrites',
               5: 'Axons',
               6: 'Dendrites',
               7: 'Axons'}

# stacked and normed to all
for r_idx, region in enumerate(regions):
    
    # get plot idx
    plot_idx = r_idx *2
    
    # polar histogram
    change_projection(fig, axs, axs[plot_idx], projection = 'polar')

    # set axis
    ax = axs[plot_idx]
    
    # subplot title
    ax.set_title(axis_titles[plot_idx], 
                 fontsize = 14, 
                 loc='left')
    
    # define bottom for stack histogram
    bottom = [0]*n_bins
    
    for n_idx, ntype in enumerate(neurite_types):

        # get occurances
        hist_angles_occu = population_orientation.loc[f'{ntype}-{region}-norm_toall']

        # plot histogram as barplot
        ax.bar(hist_bins, hist_angles_occu,
               bottom = bottom,
               width = binsize, 
               align = 'edge',
               edgecolor = 'none',
               color = neurite_color_dict[region][ntype],
               label = ntype.capitalize())
        
        # add to bottom list for next step
        bottom = np.add(bottom, hist_angles_occu)
        
    # y
    ax.set_ylim([0, 0.25])
    ax.set_yticks(ticks = [0.1, 0.2, 0.25], labels = ['', '',str(int(0.25*100)) + ' %'])
    
    # legend

    # get legend handles
    h, l = ax.get_legend_handles_labels()

    # legend
    axs[plot_idx+1].legend(h, l, 
                           title = 'Neurite type',
                           title_fontsize = 14,
                           frameon = False, 
                           ncol = 1, 
                           loc = 'center',
                           fontsize = 12)
    
    # remove spines, ticks, and ticklabels
    remove_spines_n_ticks([axs[plot_idx+1]], axis = 'x')
    remove_spines_n_ticks([axs[plot_idx+1]], axis = 'y')
    [axs[plot_idx+1].spines[spine].set_visible(False) for spine in ['top', 'right']]
    axs[plot_idx+1].set_xticks([])
    axs[plot_idx+1].set_yticks([])



# normed to type
for r_idx, region in enumerate(regions):
    
    for n_idx, ntype in enumerate(neurite_types):
    
        # get plot idx
        plot_idx = n_idx + (r_idx*2) + 4
        
        # polar histogram
        change_projection(fig, axs, axs[plot_idx], projection = 'polar')
    
        # set axis
        ax = axs[plot_idx]
        
        # subplot title
        ax.set_title(axis_titles[plot_idx], 
                     fontsize = 14, 
                     loc='left')
        
        # get occurances
        hist_angles_occu = population_orientation.loc[f'{ntype}-{region}-norm_totype']

        # plot histogram as barplot
        ax.bar(hist_bins, hist_angles_occu,
               width = binsize, 
               align = 'edge',
               edgecolor = 'none',
               color = neurite_color_dict[region][ntype])
               # label = branch_label)
            
        # y axis
        ymax = np.ceil(max(hist_angles_occu)/0.1)*0.1

        # y
        ax.set_ylim([0, ymax])
        ax.set_yticks(ticks = np.arange(0, ymax +0.05, 0.1), labels = (['']*int((ymax/0.1))) + [str(int(ymax*100)) + ' %'])
        

for ax in axs[[0,2,4,5,6,7]]:
    # x axis
    ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
    ax.set_xticklabels(['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv'])

    # set grid
    ax.grid(True, alpha = 0.5)
    
    # set grid behind plot
    ax.set_axisbelow(True)
    
# align labels
fig.align_labels()

# show figure
plt.show()

# save figure
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_polarpop_dir
save_figures(fig, 
              figure_name = 'population_orientation', 
              save_dir = cellmorph_polarpop_dir,
              darkmode_bool=darkmode_bool, 
              figure_format='both')