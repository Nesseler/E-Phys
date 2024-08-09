# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 16:51:28 2024

@author: nesseler
"""

from cellmorphology.cellmorph_packages import plt, mtl, sbn, pd, join, np

# script specific directories / parameters / functions
from parameters.directories_win import cell_morph_descrip_dir, table_file, cell_morph_plots_dir



# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# load histogram data
polar_plot_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_occurrances.xlsx'), index_col = 'orientation_rad')
polar_plot_dendrites_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_dendrites_occurrances.xlsx'), index_col = 'orientation_rad')
polar_plot_axons_occurrances = pd.read_excel(join(cell_morph_descrip_dir, 'polar_plot_axons_occurrances.xlsx'), index_col = 'orientation_rad')

# organise histogram data into dict
polar_occu_dict = {'neurites' : polar_plot_occurrances, 
                   'dendrites' : polar_plot_dendrites_occurrances,
                   'axons' : polar_plot_axons_occurrances}


orientation_labels = ['p', '', 'd', '', 'a', '', 'v', '']

bins_angles = polar_plot_occurrances.index.to_list()

# get binsize
resul_binsize = np.diff(bins_angles)[-1]

# define directory
cell_morph_plots_polar_pop_dir = join(cell_morph_plots_dir, 'polar_plots_population')

# initialise lists
neurite_types = ['neurites', 'dendrites', 'axons']
regions = ['MeA', 'BAOT']


# %% drop cells with no axon from axon occurances

# get cell_IDs of cells with axon
cell_IDs_w_axon = (polar_plot_axons_occurrances.sum() > 0)
cell_IDs_w_axon = cell_IDs_w_axon[cell_IDs_w_axon == True].index.to_list()

# redefine axon occurrances
polar_plot_axons_occurrances = polar_plot_axons_occurrances[cell_IDs_w_axon]


# %% get cell_IDs for regions

# get cell_IDs
cell_IDs = polar_plot_occurrances.columns.to_list()

region_cell_IDs = {'all' : cell_IDs,
                   'MeA' : MetaData.loc[cell_IDs, 'Region'][MetaData.loc[cell_IDs, 'Region'] == 'MeA'].index.to_list(),
                   'BAOT': MetaData.loc[cell_IDs, 'Region'][MetaData.loc[cell_IDs, 'Region'] == 'BAOT'].index.to_list()}

region_cell_IDs_w_axon = {'all' : cell_IDs_w_axon,
                          'MeA' : MetaData.loc[cell_IDs_w_axon, 'Region'][MetaData.loc[cell_IDs_w_axon, 'Region'] == 'MeA'].index.to_list(),
                          'BAOT': MetaData.loc[cell_IDs_w_axon, 'Region'][MetaData.loc[cell_IDs_w_axon, 'Region'] == 'BAOT'].index.to_list()}


# %% normalise polar occurences

# initialize dataframe for saveing normalised occurances
polar_occu_toNeurites = pd.DataFrame(index = bins_angles)
polar_occu_toType = pd.DataFrame(index = bins_angles)

for r_idx, region in enumerate(regions):
    
    for n_idx, neurite_type in enumerate(neurite_types):
        
        # get current cell_IDs
        if neurite_type != 'axons':
            cur_region_cell_IDs = region_cell_IDs[region]
        elif neurite_type == 'axons':
            cur_region_cell_IDs = region_cell_IDs_w_axon[region]
        
        # set polar occurrances dataframe
        polar_occu = polar_occu_dict[neurite_type][cur_region_cell_IDs]
        
        # set polar occurances of only neurites
        polar_occu_neurites = polar_occu_dict['neurites'][cur_region_cell_IDs]
        
        # get sum and divide all occurrances by sum
        polar_occu_normed_per_neurites = polar_occu.div(polar_occu_neurites.sum(), axis = 1)
        
        # get mean per bin
        polar_occu_normed_per_neurites_mean_per_bin = polar_occu_normed_per_neurites.mean(axis = 1)
    
        # write to Dataframe
        polar_occu_toNeurites[f'{neurite_type}-{region}'] = polar_occu_normed_per_neurites_mean_per_bin
        
        
        # get sum and divide all occurrances by sum
        polar_occu_normed_per_type = polar_occu.div(polar_occu.sum(), axis = 1)
        
        # get mean per bin
        polar_occu_normed_per_type_mean_per_bin = polar_occu_normed_per_type.mean(axis = 1)
    
        # write to Dataframe
        polar_occu_toType[f'{neurite_type}-{region}'] = polar_occu_normed_per_type_mean_per_bin
    

# %% load number of terminal points

# load excel sheet
n_terminal_points_df = pd.read_excel(join(cell_morph_descrip_dir, 'n_terminal_points.xlsx'), index_col = 'cell_ID')

# melt df to plot
n_terminal_points = pd.melt(n_terminal_points_df, var_name= 'neurite_type', value_name= 'n_terminal_points', ignore_index=False)

# get Regions
n_terminal_points['Region'] = MetaData.loc[n_terminal_points_df.index, 'Region']

# calc mean and std of number of terminal points
n_terminal_points_means = pd.DataFrame(columns = ['neuritic_terminals', 'dendritic_terminals', 'axonic_terminals'],
                                        index = ['all', 'MeA', 'BAOT'])
n_terminal_points_stds = pd.DataFrame(columns = ['neuritic_terminals', 'dendritic_terminals', 'axonic_terminals'],
                                      index = ['all', 'MeA', 'BAOT'])

# get cell_IDs in regions
all_cell_IDs = n_terminal_points_df.index.to_list()
region_cell_IDs = {'all' : all_cell_IDs,
                    'MeA' : MetaData.loc[all_cell_IDs, 'Region'][MetaData.loc[all_cell_IDs, 'Region'] == 'MeA'].index.to_list(),
                    'BAOT': MetaData.loc[all_cell_IDs, 'Region'][MetaData.loc[all_cell_IDs, 'Region'] == 'BAOT'].index.to_list()}

# calc per region
for region in ['all', 'MeA', 'BAOT']:
    n_terminal_points_means.loc[region, :] = n_terminal_points_df.loc[region_cell_IDs[region], :].mean()
    n_terminal_points_stds.loc[region, :] = n_terminal_points_df.loc[region_cell_IDs[region], :].std()


# %% initialize plotting

from functions.functions_plotting import save_figures, get_colors, get_figure_size, set_font_sizes, change_projection
from cellmorphology.cellmorph_colors import neurite_color_dict
from inspect import cleandoc as indentstring # package for multiline string with indentation

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% set dicts for plot (region, neurite_types)

# set color dict for neurite_types and regions
neurites_regions_color_dict = {'all'  : {'dendrites' : 'grey',                'axons' : 'lightgrey'},
                                'MeA'  : {'dendrites' : region_colors['MeA'],  'axons' : colors_dict['MeA_lighter']},
                                'BAOT' : {'dendrites' : region_colors['BAOT'], 'axons' : colors_dict['BAOT_lighter']}}

# set dict for titles of subplots
alpha_labels_dict = {'MeA'  : {'neurites'   : '$\mathregular{A_{i}}$',
                               'dendrites'  : '$\mathregular{A_{ii}}$',
                               'axons'      : '$\mathregular{A_{iii}}$'},
                     'BAOT' : {'neurites'   : '$\mathregular{B_{i}}$',
                               'dendrites'  : '$\mathregular{B_{ii}}$',
                               'axons'      : '$\mathregular{B_{iii}}$'},
                     'n_terminals' : {'neurites'   : '$\mathregular{C_{i}}$',
                                      'dendrites'  : '$\mathregular{C_{ii}}$',
                                      'axons'      : '$\mathregular{C_{iii}}$'}}



# %% figure

# polar histogram 
fig_norm, axs_norm = plt.subplots(nrows = 3,
                                  ncols = 3,
                                  layout = 'constrained',
                                  figsize = get_figure_size(width = 150, height = 150),
                                  dpi = 600)

# flatten axes array
axs_norm = axs_norm.flatten()

### polar plots ###

# # change projection for polar plots
# for ax in axs_norm[:6]:
#     change_projection(fig_norm, axs_norm, ax, projection = 'polar')
    
# # set polar_plot idx
# polar_idx = 0

# # loop through regions (rows)
# for region in regions:
 
#     # loop through neurite types (cols)
#     for neurite_type in neurite_types:
 
#         # set axis
#         ax = axs_norm[polar_idx]
        
#         # iterate on polar index
#         polar_idx += 1
        
#         # set axis title
#         ax.set_title(alpha_labels_dict[region][neurite_type] + ': ' + region + ' ' + neurite_type,
#                      fontsize=12, 
#                      loc='left',
#                      x = -0.2)
    

        
#         # set indices for dataframe
#         polar_occu_cols = [f'{neurite_type_ls}-{region}' for neurite_type_ls in neurite_types]
        
#         # plot differently for neurites and dendrites / axons (stack bar plot for neurites)
#         if neurite_type == 'neurites':
            
#             # dendrites normed to all neurites (bottom)
#             ax.bar(bins_angles, polar_occu_toNeurites[f'dendrites-{region}'],
#                    width = resul_binsize, 
#                    align = 'edge',
#                    edgecolor = 'none',
#                    color = neurite_color_dict[region]['dendrites'],
#                    zorder = 2)
            
#             # axons normed to all neurites (top)
#             ax.bar(bins_angles, polar_occu_toNeurites[f'axons-{region}'],
#                    bottom = polar_occu_toNeurites[f'dendrites-{region}'],
#                    width = resul_binsize, 
#                    align = 'edge',
#                    edgecolor = 'none',
#                    color = neurite_color_dict[region]['axons'],
#                    zorder = 2)
            
#         else:
#             ax.bar(bins_angles, polar_occu_toType[f'{neurite_type}-{region}'],
#                    width = resul_binsize, 
#                    align = 'edge',
#                    edgecolor = 'none',
#                    color = neurite_color_dict[region][neurite_type],
#                    zorder = 2)  
           

# ### y axis

# # neurites
# for ax in axs_norm[:6:3]:
#     ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '25 %'], zorder = 3)
#     ax.set_rlabel_position(80)

# # dendrites
# for ax in axs_norm[1:6:3]:
#     ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '25 %'], zorder = 3)
#     ax.set_rlabel_position(80)

# # axons
# for ax in axs_norm[2:6:3]:
#     ax.set_yticks(ticks = np.arange(0, 0.40 + 0.01, 0.1), labels = ['', '', '20 %', '', '40 %'], va = 'top', zorder = 3)
#     ax.set_rlabel_position(-80)


# ### x axis
# for ax in axs_norm[:9]:
#     # x axis
#     ax.set_xticks(ticks = np.arange(0, np.pi*2, np.pi / 4), labels = orientation_labels, zorder = 3)

#     # grid
#     ax.grid(True, alpha = 0.5, color = 'gray', zorder = 0)
#     ax.set_axisbelow(True)     
        

### number of terminal points ###

# set neurite_type lookup dict for terminals
neurite_terminal_types = {'neurites' : 'neuritic_terminals',
                          'dendrites' : 'dendritic_terminals',
                          'axons' : 'axonic_terminals'}

for neurite_idx, neurite_type in enumerate(neurite_types):
    
    # set axis
    ax = axs_norm[neurite_idx + 6]  
    
    # set axis title
    ax.set_title(alpha_labels_dict['n_terminals'][neurite_type] + ': ' + neurite_type.title(),
                 fontsize=12, 
                 loc='left',
                 x = -0.2)
    
    # # set violin data
    # # get data for both regions
    # violin_data_bothregions = n_terminal_points[n_terminal_points['neurite_type'] == neurite_terminal_types[neurite_type]]
    # violin_data_bothregions.loc[:, 'Region'] = ['both_regions'] * violin_data_bothregions.shape[0]

    # get data for individual regions
    violin_data = n_terminal_points[n_terminal_points['neurite_type'] == neurite_terminal_types[neurite_type]]
    violin_data = violin_data[violin_data['Region'] != 'BAOT/MeA']
    
    # concat data
    # violin_data = pd.concat([violin_data, violin_data_bothregions], axis = 0)
    violin_data['X'] = [0] * violin_data.shape[0]
      
    # violin plots   
    violins = sbn.violinplot(data = violin_data,
                             x = 'X',
                             y = 'n_terminal_points',
                             hue = 'Region', 
                             # hue_order = ['MeA', 'BAOT'],
                             split = True,
                             gap = .1)
    
    # edit lines of quarts
    for r_idx, region in enumerate(['MeA', 'BAOT']):
        
        # for l_idx in np.arange(0, 3):
        
        #     all_l_idx = r_idx * 3 + l_idx
            
        #     # violins.lines[all_l_idx].set_color(neurite_color_dict[region][neurite_type])
        #     violins.lines[all_l_idx].set_color(colors_dict['primecolor'])
        #     violins.lines[all_l_idx].set_linestyle('solid')
        
        # set edge color of violin
        # violins.collections[r_idx].set_edgecolor(neurite_color_dict[region][neurite_type])
        violins.collections[r_idx].set_edgecolor('None')
        
        # set facecolor of violin
        # violins.collections[r_idx].set_facecolor('None')
        violins.collections[r_idx].set_facecolor(neurite_color_dict[region][neurite_type])
    
    
#     # swarmplots    
#     sbn.swarmplot(data = violin_data,
#                   x = 'Region', 
#                   y = 'n_terminal_points',
#                   ax = ax,
#                   s = 1.5,
#                   color=colors_dict['primecolor'],
#                   order = ['both_regions', 'MeA', 'BAOT'])
        
#     # errorbar
#     for r_idx, region in enumerate(['all', 'MeA', 'BAOT']):
#         ax.errorbar(x = r_idx+0.3,
#                     y = n_terminal_points_means.at[region, neurite_terminal_types[neurite_type]],
#                     yerr = n_terminal_points_stds.at[region, neurite_terminal_types[neurite_type]],
#                     fmt='_', 
#                     markersize = 6,
#                     markerfacecolor = 'none',
#                     capsize = 2,
#                     color=neurite_color_dict[region][neurite_type],
#                     linewidth = 1,
#                     label = '_nolegend_')
    
    
#     # edit axes
#     # x
#     ax.set_xlim(0-0.5, 2+0.5)
#     ax.set_xticks(ticks=np.arange(0, 2 + 0.5, 1),
#                   labels=['Both\nregions', 'MeA', 'BAOT'], rotation=60)
#     ax.spines['bottom'].set_bounds([0, 2])
#     ax.set_xlabel('')
    
#     # edit seaborn legend
#     ax.legend().set_visible(False)
    
    
#     # y
#     if neurite_type == 'neurites' or neurite_type == 'dendrites':
#         ymin = 0
#         ymax = 45
#         ystep = 10
#         ystepminor = 5
    
#     elif neurite_type == 'axons':
#         ymin = 0
#         ymax = 10
#         ystep = 5
#         ystepminor = 1
    
    
#     ypad = (ymax - ymin) * 0.05
#     ax.set_ylim(ymin-ypad, ymax+ypad)
#     ax.set_yticks(ticks=np.arange(ymin, ymax+1, ystep))
#     ax.set_yticks(ticks=np.arange(ymin, ymax + ystepminor, ystepminor), minor=True)
#     ax.spines['left'].set_bounds([ymin, ymax])
    
    
#     # remove spines
#     [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
#     # axis labels
#     ax.set_ylabel('Number of terminal\nbranches [#]')
#     ax.set_xlabel('Region')

# # align labels
# fig_norm.align_labels()

# show figure
plt.show()


