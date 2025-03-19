# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:17:07 2024

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

# define sholl radius step size
sholl_step_size = 1 # µm


# %% load data

from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_analysis_dir, cellmorph_metrics_dir, cellmorph_shollfigs_dir

# get mean profiles
avg_sholl_profiles = pd.read_excel(join(cellmorph_analysis_dir, 'sholl-combined_profiles', 'sholl_profiles-avg.xlsx'),
                                   index_col = 'Radius')

# get sholl metrics
sholl_metrics = pd.read_excel(join(cellmorph_metrics_dir, 'sholl_metrics.xlsx'),
                              index_col = 'cell_ID')


# %% cell_IDs dicts

# get all cell_IDs
cell_IDs = sholl_metrics.index.to_list()

# get cell_IDs
cell_IDs_dict = {'all': cell_IDs,
                 'MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in cell_IDs],
                 'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in cell_IDs],
                 'BAOT/MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list() if cell_ID in cell_IDs]}

# get cell_IDs for cells with axon
cell_IDs_w_axon_dict = {'all': sholl_metrics['critical_radius-axons'].dropna().index.to_list(),
                        'MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in sholl_metrics['critical_radius-axons'].dropna().index.to_list()],
                        'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in sholl_metrics['critical_radius-axons'].dropna().index.to_list()],
                        'BAOT/MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT/MeA'].index.to_list() if cell_ID in sholl_metrics['critical_radius-axons'].dropna().index.to_list()]}


# %% sholl profiles figure 

# initialise figure
fig_sholl, axs_sholl = plt.subplot_mosaic('BBCC;EEFF',#'BD;BD;BE;CE;CF;CF',
                                          figsize=get_figure_size(height=125.103, width=260.334),
                                          layout='tight',
                                          height_ratios = [2.5,1],
                                          dpi=600)

# flatten nd array of axes
# axs_sholl = axs_sholl.flatten()
axs_keys = {'BAOT': 'B', 'MeA': 'C'}

# specify line
line_dict = {'lw': 1, 'alpha': 1}

# set labels for subplots
axs_titles = {'MeA' : 'B: MeA',
              'BAOT': 'A: BAOT',
              'neurites' : '$\mathregular{C_{i}}$: Neurites',
              'dendrites' : '$\mathregular{C_{i}}$: Dendrites',
              'axons' : '$\mathregular{C_{ii}}$: Axons'}

# plot all profiles together
for region in regions:
    for neurite_type in ['dendrites', 'axons']:

        # set axis
        ax = axs_sholl[axs_keys[region]]
        
        # calc upper and lower border
        y_bottom = avg_sholl_profiles[f'mean-sholl_profile-{neurite_type}-{region}'] - avg_sholl_profiles[f'std-sholl_profile-{neurite_type}-{region}']
        y_top = avg_sholl_profiles[f'mean-sholl_profile-{neurite_type}-{region}'] + avg_sholl_profiles[f'std-sholl_profile-{neurite_type}-{region}']

        # plot std as shade
        ax.fill_between(x = avg_sholl_profiles.index.to_list(),
                        y1 = y_bottom,
                        y2 = y_top,
                        color = neurite_color_dict[region][neurite_type],
                        alpha = 0.5,
                        linewidth = 0,
                        zorder = 1,
                        label = '_nolegend_')

        # plot mean
        ax.plot(avg_sholl_profiles[f'mean-sholl_profile-{neurite_type}-{region}'],
                color = neurite_color_dict[region][neurite_type],
                zorder = 3,
                label = neurite_type.title())

# edit sholl profile axes
for region in regions:

    # set axis
    ax = axs_sholl[axs_keys[region]]

    # subplot title
    ax.set_title(axs_titles[region], fontsize=14, loc='left')

    # legend
    ax.legend(fontsize = 12,
              frameon = False,
              title = 'Neurite types',
              title_fontsize = 14)

    # y
    ymin = 0
    ymax = 20
    ystep = 5
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# comparisons of neurite types between regions
axs_keys_neurites = {'dendrites': 'E', 'axons': 'F'}

# add color for 'all' to region_colors
region_colors['all'] = 'gray'

for neurite_type in neurite_types:

    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]

    for region in ['BAOT', 'MeA']:

        ax.set_title(axs_titles[neurite_type], fontsize=14, loc='left')

        ax.plot(avg_sholl_profiles[f'mean-sholl_profile-{neurite_type}-{region}'],
                color = neurite_color_dict[region][neurite_type],
                zorder = 3,
                label  =region)

    # legend
    ax.legend(fontsize = 12,
              frameon = False,
              title = 'Region',
              title_fontsize = 14)


# configure axis of neurites and dendrites
for neurite_type in ['dendrites']:

    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]

    # y
    ymin = 0
    ymax = 15
    ystep = 15
    ystepminor = 5
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


for neurite_type in ['axons']:

    # set axis
    ax = axs_sholl[axs_keys_neurites[neurite_type]]

    # y
    ymin = 0
    ymax = 4
    ystep = 4
    ystepminor = 1
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystep, ystep))
    ax.set_yticks(ticks=np.arange(
        ymin, ymax + ystepminor, ystepminor), minor=True)

    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# set all x axis
for key in ['B', 'C', 'E', 'F']:
    # set axis
    ax = axs_sholl[key]

    # x
    xmin = 0
    xmax = 400
    xstep = 100
    xstepminor = 25
    ax.set_xlim(xmin-10, xmax-10)
    ax.set_xticks(ticks=np.arange(xmin, xmax + 1, xstep))
    ax.set_xticks(ticks=np.arange(
        xmin, xmax + ystepminor, xstepminor), minor=True)
    ax.spines['bottom'].set_bounds([xmin, xmax])


# axis labels
[axs_sholl[key].set_ylabel('Number of\nintersections [#]') for key in ['B', 'E']]
[axs_sholl[key].set_xlabel('Radius  [µm]') for key in ['B', 'C', 'E', 'F']]

# align all axis labels
fig_sholl.align_labels()

# show figure
plt.show()

# save figure
save_figures(fig_sholl, 
              figure_name = 'sholl_profiles_figure', 
              save_dir = cellmorph_shollfigs_dir,
              darkmode_bool=darkmode_bool, 
              figure_format='both')

