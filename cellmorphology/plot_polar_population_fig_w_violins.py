# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 19:01:10 2024

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

orientation_labels = ['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv']

bins_angles = polar_plot_occurrances.index.to_list()

# get binsize
resul_binsize = np.diff(bins_angles)[-1]

# define directory
cell_morph_plots_polar_pop_dir = join(cell_morph_plots_dir, 'polar_plots_population')


# %% min max normalised version

### normalize polar occurrences
# all
polar_plot_occurrances_normed = polar_plot_occurrances.div(polar_plot_occurrances.sum(), axis = 1)

# dendrites
polar_plot_dendrites_occurrances_normed = polar_plot_dendrites_occurrances.div(polar_plot_dendrites_occurrances.sum(), axis = 1)
pp_occu_dendrites_occu_normed_to_neurites_mean = polar_plot_dendrites_occurrances.div(polar_plot_occurrances.sum(), axis = 1).mean(axis = 1)

# axons
polar_plot_axons_occurrances_normed = polar_plot_axons_occurrances.div(polar_plot_axons_occurrances.sum(), axis = 1)
pp_occu_axons_occu_normed_to_neurites_mean = polar_plot_axons_occurrances.div(polar_plot_occurrances.sum(), axis = 1).mean(axis = 1)


### create average per bin
# all
polar_plot_occurrances_normed_mean = polar_plot_occurrances_normed.mean(axis = 1)
polar_plot_occurrances_normed_std = polar_plot_occurrances_normed.std(axis = 1)

# dendrites
polar_plot_dendrites_occurrances_normed_mean = polar_plot_dendrites_occurrances_normed.mean(axis = 1)
polar_plot_dendrites_occurrances_normed_std = polar_plot_dendrites_occurrances.div(polar_plot_occurrances.sum(), axis = 1).std(axis = 1)

# axons
polar_plot_axons_occurrances_normed_mean = polar_plot_axons_occurrances_normed.mean(axis = 1)
polar_plot_axons_occurrances_normed_std = polar_plot_axons_occurrances.div(polar_plot_occurrances.sum(), axis = 1).std(axis = 1)


### combine to one dataframe
ALL_pp_normed_totype_mean_occs = pd.DataFrame({'neurites' : polar_plot_occurrances_normed_mean.to_list(),
                                                'dendrites' : polar_plot_dendrites_occurrances_normed_mean.to_list(), 
                                                'axons' : polar_plot_axons_occurrances_normed_mean.to_list()}, 
                                               index = orientation_labels)

ALL_pp_normed_toneurites_mean_occs = pd.DataFrame({'neurites' : polar_plot_occurrances_normed_mean.to_list(),
                                                   'dendrites' : pp_occu_dendrites_occu_normed_to_neurites_mean.to_list(), 
                                                   'axons' : pp_occu_axons_occu_normed_to_neurites_mean.to_list()}, 
                                                  index = orientation_labels)

ALL_pp_normed_toneurites_std_occs = pd.DataFrame({'neurites' : polar_plot_occurrances_normed_std.to_list(),
                                                  'dendrites' : polar_plot_dendrites_occurrances_normed_std.to_list(), 
                                                  'axons' : polar_plot_axons_occurrances_normed_std.to_list()}, 
                                                 index = orientation_labels)



### create normalized polar occurrences per region


def region_polar_plot_occurrences(pp_all, pp_dendrites, pp_axons, region, neurite_types, orientation_labels):
    
# pp_all = polar_plot_occurrances
# pp_dendrites = polar_plot_dendrites_occurrances
# pp_axons = polar_plot_axons_occurrances
# region = 'BAOT'
# neurite_types = ['all', 'dendrites', 'axons']
# orientation_labels = orientation_labels

    # define output dataframe
    pp_normed_totype_mean_occs = pd.DataFrame(index = orientation_labels, columns = neurite_types)
    pp_normed_toneurites_mean_occs = pd.DataFrame(index = orientation_labels, columns = neurite_types)
    
    # get cell IDs in region
    region_cell_IDs = MetaData[MetaData['Region'] == region].index.to_list()
    region_cell_IDs = [c for c in region_cell_IDs if c in pp_all.columns]
    
    # loop through neurite types
    for occurrence_df, neurite_type in zip([pp_all, pp_dendrites, pp_axons], neurite_types):
        
        # limit dataframe to cells within region
        region_occ = occurrence_df[region_cell_IDs]
        region_neurite_occ = pp_all[region_cell_IDs]
        
        
        ## to type
        # get normed occurrences
        region_occ_normed = region_occ.div(region_occ.sum(), axis = 1)
        
        # create average per bin
        region_occ_normed_mean = region_occ_normed.mean(axis = 1).to_list()
    
        # add to dataframe
        pp_normed_totype_mean_occs[neurite_type] = region_occ_normed_mean
        
        
        ## to neurites
        # get normed occurrences
        region_occ_normed_mea_toneurites = region_occ.div(region_neurite_occ.sum(), axis = 1).mean(axis = 1).to_list()
    
        # add to dataframe
        pp_normed_toneurites_mean_occs[neurite_type] = region_occ_normed_mea_toneurites
        
        
    # return dataframe
    return pp_normed_totype_mean_occs, pp_normed_toneurites_mean_occs


# get normed mean occurrences for MeA
MeA_pp_normed_totype_mean_occs, MeA_pp_normed_toneurites_mean_occs = region_polar_plot_occurrences(pp_all = polar_plot_occurrances,
                                                                                                    pp_dendrites = polar_plot_dendrites_occurrances,
                                                                                                    pp_axons = polar_plot_axons_occurrances,
                                                                                                    region = 'MeA',
                                                                                                    neurite_types = ['neurites', 'dendrites', 'axons'],
                                                                                                    orientation_labels = orientation_labels)

# get normed mean occurrences for BAOT
BAOT_pp_normed_totype_mean_occs, BAOT_pp_normed_toneurites_mean_occs = region_polar_plot_occurrences(pp_all = polar_plot_occurrances,
                                                                                                     pp_dendrites = polar_plot_dendrites_occurrances,
                                                                                                     pp_axons = polar_plot_axons_occurrances,
                                                                                                     region = 'BAOT',
                                                                                                     neurite_types = ['neurites', 'dendrites', 'axons'],
                                                                                                     orientation_labels = orientation_labels)

# %% load number of terminal points

# load excel sheet
n_terminal_points_df = pd.read_excel(join(cell_morph_descrip_dir, 'n_terminal_points.xlsx'), index_col = 'cell_ID')

# melt df to plot
n_terminal_points = pd.melt(n_terminal_points_df, var_name= 'neurite_type', value_name= 'n_terminal_points', ignore_index=False)

# get Regions
n_terminal_points['Region'] = MetaData.loc[n_terminal_points_df.index, 'Region']

# remove number of terminal points of axonic
# n_terminal_points = n_terminal_points[n_terminal_points[n_terminal_points['neurite_type'] == 'axonic_terminals'] != 0]


# %% initialize plotting

from functions.functions_plotting import save_figures, get_colors, get_figure_size, set_font_sizes, change_projection
from cellmorphology.cellmorph_colors import neurite_color_dict

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})

# initialise lists
neurite_types = ['neurites', 'dendrites', 'axons']
regions = ['all', 'MeA', 'BAOT']


# %% figure

# polar histogram 
fig_norm, axs_norm = plt.subplots(nrows = 4,
                                  ncols = 3,
                                  layout = 'constrained',
                                  figsize = get_figure_size(width = 150, height = 200),
                                  dpi = 600)

# flatten axes array
axs_norm = axs_norm.flatten()

# ### polar plots ###

# # change projection for polar plots
# for ax in axs_norm[:9]:
#     change_projection(fig_norm, axs_norm, ax, projection = 'polar')
    
# # set dicts for plot (region, neurite_types)
# occs_dict = {'all' : {'neurites'  : ALL_pp_normed_toneurites_mean_occs,
#                       'dendrites' : ALL_pp_normed_totype_mean_occs['dendrites'],
#                       'axons'     : ALL_pp_normed_totype_mean_occs['axons']},
#              'MeA' : {'neurites'  : MeA_pp_normed_toneurites_mean_occs,
#                       'dendrites' : MeA_pp_normed_totype_mean_occs['dendrites'],
#                        'axons'    : MeA_pp_normed_totype_mean_occs['axons']},
#              'BAOT': {'neurites'  : BAOT_pp_normed_toneurites_mean_occs,
#                       'dendrites' : BAOT_pp_normed_totype_mean_occs['dendrites'],
#                        'axons'    : BAOT_pp_normed_totype_mean_occs['axons']}}


# # set color dict for neurite_types and regions
# neurites_regions_color_dict = {'all'  : {'dendrites' : 'grey',                'axons' : 'lightgrey'},
#                                'MeA'  : {'dendrites' : region_colors['MeA'],  'axons' : colors_dict['MeA_lighter']},
#                                'BAOT' : {'dendrites' : region_colors['BAOT'], 'axons' : colors_dict['BAOT_lighter']}}

# # set dict for titles of subplots
# alpha_labels_dict = {'all' : {'neurites'   : '$\mathregular{A_{i}}$',
#                               'dendrites'  : '$\mathregular{A_{ii}}$',
#                               'axons'      : '$\mathregular{A_{iii}}$',
#                               'n_terminals': '$\mathregular{A_{iv}}$',},
#                      'MeA' : {'neurites'   : '$\mathregular{B_{i}}$',
#                               'dendrites'  : '$\mathregular{B_{ii}}$',
#                               'axons'      : '$\mathregular{B_{iii}}$',
#                               'n_terminals': '$\mathregular{B_{iv}}$',},
#                      'BAOT': {'neurites'   : '$\mathregular{C_{i}}$',
#                               'dendrites'  : '$\mathregular{C_{ii}}$',
#                               'axons'      : '$\mathregular{C_{iii}}$',
#                               'n_terminals': '$\mathregular{C_{iv}}$'}}

# # set polar_plot idx
# polar_idx = 0


# # loop through neurite types (rows)
# for neurite_type in neurite_types:
    
#     # loop through regions (cols)
#     for region in regions:
        
#         # set axis
#         ax = axs_norm[polar_idx]
        
#         # set title
#         if region == 'all':
#             region_label = 'Both regions'
#         else:
#             region_label = region
#         ax.set_title(alpha_labels_dict[region][neurite_type] + ': ' + region_label + ' ' + neurite_type,
#                      fontsize=9, 
#                      loc='left',
#                      x = -0.2)
        
#         # iterate on polar index
#         polar_idx += 1
        
#         # set plotting data
#         occu = occs_dict[region][neurite_type]
        
#         # plot differently for neurites and dendrites / axons (stack bar plot for neurites)
#         if neurite_type == 'neurites':
            
#             # dendrites normed to all neurites (bottom)
#             ax.bar(bins_angles, occs_dict[region]['neurites']['dendrites'],
#                    width = resul_binsize, 
#                    align = 'edge',
#                    edgecolor = 'none',
#                    color = neurite_color_dict[region]['dendrites'])
            
#             # axons normed to all neurites (top)
#             ax.bar(bins_angles, occs_dict[region]['neurites']['axons'],
#                    bottom = occs_dict[region]['neurites']['dendrites'],
#                    width = resul_binsize, 
#                    align = 'edge',
#                    edgecolor = 'none',
#                    color = neurite_color_dict[region]['axons'])
            
#         else:
#             ax.bar(bins_angles, occu,
#                    width = resul_binsize, 
#                    align = 'edge',
#                    edgecolor = 'none',
#                    color = neurite_color_dict[region][neurite_type])  
           
# # adjust axis title position
# # plt.rcParams['axes.titlex'] = -0.1 
        
# # # column header
# # for ax, neurite_type in zip(axs_norm[::3], ['Neurites', 'Dendrites', 'Axons']):
# #     ax.annotate(neurite_type, xy=(0, 0.5), xytext=(-20, 0),
# #                 xycoords=ax.yaxis.label, textcoords='offset points',
# #                 ha='center', va='center', fontsize = 9, rotation = 90)

    
# # # row header
# # for ax, region_label in zip(axs_norm[:3], ['Both regions', 'MeA', 'BAOT']):
# #     ax.annotate(region_label, xy = (0.5, 1), xytext = (0, 30), 
# #                 xycoords='axes fraction', textcoords='offset points',
# #                 ha='center', va='baseline', fontsize = 9)


# ### y axis

# # neurites
# for ax in axs_norm[:3]:
#     ax.set_yticks(ticks = np.arange(0, 0.30 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '', '30 %'])
#     ax.set_rlabel_position(80)

# # dendrites
# for ax in axs_norm[3:6]:
#     ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '10 %', '', '20 %', ''])
#     ax.set_rlabel_position(80)

# # axons
# for ax in axs_norm[6:9]:
#     ax.set_yticks(ticks = np.arange(0, 0.35 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '', '30 %', ''], va = 'top')
#     ax.set_rlabel_position(-80)


# ### x axis
# for ax in axs_norm[:9]:
#     # x axis
#     ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
#     ax.set_xticklabels(orientation_labels)


#     # grid
#     ax.grid(True, alpha = 0.5, color = 'gray')
#     ax.set_axisbelow(True)     
        

### number of terminal points ###

for r_idx, region in enumerate(regions):
    
    if region == 'all':
        violin_data = n_terminal_points
    else:
        violin_data = n_terminal_points[n_terminal_points['Region'] == region]
    
    sbn.violinplot(data = violin_data,
                   x = 'neurite_type', 
                   y = 'n_terminal_points',
                   ax = axs_norm[r_idx + 9], 
                   scale='width')                       # will be density_norm in v0.13.0
    
    sbn.swarmplot(data = violin_data,
                   x = 'neurite_type', 
                   y = 'n_terminal_points',
                   ax = axs_norm[r_idx + 9],
                   s = 1)
                   # color=)



# align labels
fig_norm.align_labels()

# show figure
plt.show()




# %%

# flatten axis
axs_norm = axs_norm.flatten()

# set colors for dendrites and axons in both regions
both_regions_colors = ['grey', 'lightgrey']

## neurites (dendrites and axons stacked)
for plot_i, to_neurites_occu, color_pair in zip(np.arange(0, 9, 3), [ALL_pp_normed_toneurites_mean_occs, MeA_pp_normed_toneurites_mean_occs, BAOT_pp_normed_toneurites_mean_occs], [both_regions_colors, [region_colors['MeA'], colors_dict['MeA_lighter']], [region_colors['BAOT'], colors_dict['BAOT_lighter']]]):
    
    axs_norm[plot_i].bar(bins_angles, to_neurites_occu['dendrites'],
                          width = resul_binsize, 
                          align = 'edge',
                          edgecolor = 'none',
                          color = color_pair[0])
    
    axs_norm[plot_i].bar(bins_angles, to_neurites_occu['axons'],
                         bottom = to_neurites_occu['dendrites'],
                         width = resul_binsize, 
                         align = 'edge',
                         edgecolor = 'none',
                         color = color_pair[1])
    
# test standard deviation

# axs_norm[0].errorbar(x = np.add(bins_angles, ((2*np.pi)/(len(bins_angles))/2)),
#                      y = ALL_pp_normed_totype_mean_occs['all'],
#                      yerr = ALL_pp_normed_toneurites_std_occs['all'],
#                      fmt = '.',
#                      elinewidth = 0.5,
#                      color = colors_dict['primecolor'],
#                      markersize = 2,
#                      alpha = 0.5)





### regions
for plot_i, neurite_type in enumerate(['dendrites', 'axons']):
    ## both
    axs_norm[plot_i + 1].bar(bins_angles, ALL_pp_normed_totype_mean_occs[neurite_type],
                             width = resul_binsize, 
                             align = 'edge',
                             edgecolor = 'none',
                             color = neurites_regions_color_dict['all'][neurite_type])  
        
    ## MeA
    # plot histogram as barplot
    axs_norm[plot_i + 1 + 3].bar(bins_angles, MeA_pp_normed_totype_mean_occs[neurite_type],
                                 width = resul_binsize, 
                                 align = 'edge',
                                 edgecolor = 'none',
                                 color = neurites_regions_color_dict['MeA'][neurite_type])

    ## BAOT
    # plot histogram as barplot
    axs_norm[plot_i + 1 + 6].bar(bins_angles, BAOT_pp_normed_totype_mean_occs[neurite_type],
                                 width = resul_binsize, 
                                 align = 'edge',
                                 edgecolor = 'none',
                                 color = neurites_regions_color_dict['BAOT'][neurite_type])



# column header
for ax, neurite_type in zip(axs_norm[:3], ['Neurites', 'Dendrites', 'Axons']):
    ax.annotate(neurite_type, xy = (0.5, 1), xytext = (0, 30), 
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline', fontsize = 9)
    
# row header
for ax, plot_type in zip(axs_norm[::3], ['Both regions', 'MeA', 'BAOT']):
    ax.annotate(plot_type, xy=(0, 0.5), xytext=(-20, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                ha='center', va='center', fontsize = 9, rotation = 90)


### y axis

# neurites
for ax in axs_norm[::3]:
    ax.set_yticks(ticks = np.arange(0, 0.30 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '', '30 %'])
    ax.set_rlabel_position(80)

# dendrites
for ax in axs_norm[1::3]:
    ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '10 %', '', '20 %', ''])
    ax.set_rlabel_position(80)

# axons
for ax in axs_norm[2::3]:
    ax.set_yticks(ticks = np.arange(0, 0.35 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '', '30 %', ''], va = 'top')
    ax.set_rlabel_position(-80)


### x axis
for ax in axs_norm:
    # x axis
    ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
    ax.set_xticklabels(orientation_labels)


    # grid
    ax.grid(True, alpha = 0.5, color = 'gray')
    ax.set_axisbelow(True)
    






set_font_sizes(9, 9)


# save_figures(fig_norm, figure_name = 'population_polar_plots-normed_per_type-regions', save_dir = join(cell_morph_plots_dir, 'polar_plots_population'),
#              figure_format= 'both')



# %% exports occus to excel sheet

for occu_df, region_label in zip([ALL_pp_normed_totype_mean_occs, MeA_pp_normed_totype_mean_occs, BAOT_pp_normed_totype_mean_occs], ['All', 'MeA', 'BAOT']):

    occu_df.to_excel(join(cell_morph_plots_dir, 
                          'polar_plots_population', 
                          'population_polar_plots-normed_totype-regions-' + region_label + '-means.xlsx'), 
                     index = orientation_labels)
    
    
for occu_df, region_label in zip([ALL_pp_normed_toneurites_mean_occs, MeA_pp_normed_toneurites_mean_occs, BAOT_pp_normed_toneurites_mean_occs], ['All', 'MeA', 'BAOT']):

    occu_df.to_excel(join(cell_morph_plots_dir, 
                          'polar_plots_population', 
                          'population_polar_plots-normed_toneurite-regions-' + region_label + '-means.xlsx'), 
                     index = orientation_labels)    
    
