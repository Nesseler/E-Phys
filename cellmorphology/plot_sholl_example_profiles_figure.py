# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 16:43:20 2024

@author: nesseler
"""

# general packages
from cellmorphology.cellmorph_colors import neurite_color_dict
from functions.functions_plotting import get_colors, save_figures, set_font_sizes, get_figure_size
from cellmorphology.cellmorph_packages import plt, mtl, sbn, pd, join, np

# script specific directories / parameters / functions
from parameters.directories_win import cell_morph_descrip_dir, table_file, cell_morph_plots_dir

from cellmorphology.cellmorph_parameters import sholl_step_size, sholl_profile_range
from cellmorphology.cellmorph_colors import neurite_color_dict

# initialise lists
neurite_types = ['neurites', 'dendrites', 'axons']
regions = ['MeA', 'BAOT']

# example cell_IDs
cell_IDs = {'MeA' : ['E-140', 'E-147'],
            'BAOT' : ['E-076', 'E-089']}

# get all cell_IDs from dict
all_cell_IDs = [cell_ID for region in regions for cell_ID in cell_IDs[region]]

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# %% load data

# initialise dict & dataframe for sholl profiles
sholl_profiles = {'neurites' : pd.DataFrame(columns = all_cell_IDs, index = sholl_profile_range),
                  'dendrites' : pd.DataFrame(columns = all_cell_IDs, index = sholl_profile_range),
                  'axons' : pd.DataFrame(columns = all_cell_IDs, index = sholl_profile_range)}

# loop through neurite types
for neurite_type in neurite_types:
    
    # set index name
    sholl_profiles[neurite_type].index.name = 'Radius'
    
    # define path to all profiles file
    sholl_profiles_path = join(cell_morph_descrip_dir, 'sholl_profiles_' + neurite_type + '.xlsx')
    
    # load all sholl profiles
    sholl_profiles[neurite_type] = pd.read_excel(sholl_profiles_path, index_col='Radius')
    
    
# %% initialize plotting

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% create figure

fig_ex, axs_ex = plt.subplots(nrows = 4,
                              ncols = 1,
                              figsize = get_figure_size(100, 225),
                              dpi = 600,
                              layout = 'constrained')

# adjust padding and space between subplots
fig_ex.set_constrained_layout_pads(w_pad=0.1, h_pad= 0.0, hspace=0.10, wspace=0.2)

# set list of alphabetical labels
alpha_labels = ['$\mathregular{A_{ii}}$', '$\mathregular{B_{ii}}$', '$\mathregular{C_{ii}}$', '$\mathregular{D_{ii}}$']

# loop through example cell_IDs
for c_idx, cell_ID in enumerate(all_cell_IDs):

    # set region of cell_ID
    region = MetaData.at[cell_ID, 'Region']
    
    # loop through neurite_types
    for neurite_type in neurite_types:
        
        # set dataframe to plot from
        neurite_sholl = sholl_profiles[neurite_type]
        
        # set axis
        ax = axs_ex[all_cell_IDs.index(cell_ID)]
        
        # subplot title
        ax.set_title(f'{alpha_labels[c_idx]}: {region} {cell_ID}', fontsize=12, loc='left')
        
        if cell_ID in neurite_sholl.columns.to_list():
            
            # plot sholl profile
            ax.plot(sholl_profile_range, neurite_sholl[cell_ID],
                    c = neurite_color_dict[region][neurite_type],
                    label = neurite_type)
        
        # set legend
        ax.legend(prop={'size': 9})
        
        # edit axes
        # x
        xmin = 0
        xmax = 450
        xstep = 100
        xstepminor = 25
        ax.set_xlim(xmin-10, xmax-10)
        ax.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep))
        ax.set_xticks(ticks=np.arange(xmin, xmax + xstepminor, xstepminor), minor=True)
        ax.spines['bottom'].set_bounds([xmin, xmax])
        
        # y
        ymin = 0
        ymax = 30
        ystep = 10
        ystepminor = 2
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
        ax.set_yticks(ticks=np.arange(ymin, ymax + ystepminor, ystepminor), minor=True)
        
        # remove spines
        [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# axis labels
[ax.set_ylabel('Number of\nintersections [#]', fontsize = 9) for ax in axs_ex]
[ax.set_xlabel('Radius  [µm]', fontsize = 9) for ax in axs_ex]


# align all axis labels
fig_ex.align_labels()

# show figure
plt.show()

# save figure
save_figures(fig_ex, 'sholl_profiles_cell_examples_figure', join(cell_morph_plots_dir, 'sholl_plots'),
             darkmode_bool=darkmode_bool, figure_format='both')


