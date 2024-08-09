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

orientation_labels = ['p', '', 'd', '', 'a', '', 'v', '']

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

# initialise lists
neurite_types = ['neurites', 'dendrites', 'axons']
regions = ['all', 'MeA', 'BAOT']


# %% set dicts for plot (region, neurite_types)
occs_dict = {'all' : {'neurites'  : ALL_pp_normed_toneurites_mean_occs,
                      'dendrites' : ALL_pp_normed_totype_mean_occs['dendrites'],
                      'axons'     : ALL_pp_normed_totype_mean_occs['axons']},
              'MeA' : {'neurites'  : MeA_pp_normed_toneurites_mean_occs,
                      'dendrites' : MeA_pp_normed_totype_mean_occs['dendrites'],
                        'axons'    : MeA_pp_normed_totype_mean_occs['axons']},
              'BAOT': {'neurites'  : BAOT_pp_normed_toneurites_mean_occs,
                      'dendrites' : BAOT_pp_normed_totype_mean_occs['dendrites'],
                        'axons'    : BAOT_pp_normed_totype_mean_occs['axons']}}


# set color dict for neurite_types and regions
neurites_regions_color_dict = {'all'  : {'dendrites' : 'grey',                'axons' : 'lightgrey'},
                                'MeA'  : {'dendrites' : region_colors['MeA'],  'axons' : colors_dict['MeA_lighter']},
                                'BAOT' : {'dendrites' : region_colors['BAOT'], 'axons' : colors_dict['BAOT_lighter']}}

# set dict for titles of subplots
alpha_labels_dict = {'all' : {'neurites'   : '$\mathregular{A_{i}}$',
                              'dendrites'  : '$\mathregular{A_{ii}}$',
                              'axons'      : '$\mathregular{A_{iii}}$'},
                      'MeA' : {'neurites'   : '$\mathregular{B_{i}}$',
                              'dendrites'  : '$\mathregular{B_{ii}}$',
                              'axons'      : '$\mathregular{B_{iii}}$'},
                      'BAOT': {'neurites'   : '$\mathregular{C_{i}}$',
                              'dendrites'  : '$\mathregular{C_{ii}}$',
                              'axons'      : '$\mathregular{C_{iii}}$'},
                      'n_terminals' : {'neurites'   : '$\mathregular{D_{i}}$',
                                        'dendrites'  : '$\mathregular{D_{ii}}$',
                                        'axons'      : '$\mathregular{D_{iii}}$'}}

# set spaces dict for title of axes 
spaces_dict = {'neurites' : '     ', 'dendrites' : '      ', 'axons' : '       '}

# %% figure

# polar histogram 
fig_norm, axs_norm = plt.subplots(nrows = 4,
                                  ncols = 3,
                                  layout = 'constrained',
                                  figsize = get_figure_size(width = 150, height = 220),
                                  dpi = 600)

# flatten axes array
axs_norm = axs_norm.flatten()

### polar plots ###

# change projection for polar plots
for ax in axs_norm[:9]:
    change_projection(fig_norm, axs_norm, ax, projection = 'polar')
    
# set polar_plot idx
polar_idx = 0

# loop through regions (rows)
for region in regions:
 
    # loop through neurite types (cols)
    for neurite_type in neurite_types:
 
        # set axis
        ax = axs_norm[polar_idx]
        
        # set axis title region label
        if region == 'all':
            region_label = 'Both regions'
            
            # set title label
            title_label = f"{alpha_labels_dict[region][neurite_type]}: {region_label}\n{spaces_dict[neurite_type]}{neurite_type}"
            
            # set axis title
            ax.set_title(title_label,
                          fontsize=12, 
                          loc='left',
                          x = -0.4)           
            
        else:
            region_label = region
            # set axis title
            ax.set_title(alpha_labels_dict[region][neurite_type] + ': ' + region_label + ' ' + neurite_type,
                          fontsize=12, 
                          loc='left',
                          x = -0.4)
        
        # iterate on polar index
        polar_idx += 1
        
        # set plotting data
        occu = occs_dict[region][neurite_type]
        
        # plot differently for neurites and dendrites / axons (stack bar plot for neurites)
        if neurite_type == 'neurites':
            
            # dendrites normed to all neurites (bottom)
            ax.bar(bins_angles, occs_dict[region]['neurites']['dendrites'],
                    width = resul_binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = neurite_color_dict[region]['dendrites'],
                    zorder = 2)
            
            # axons normed to all neurites (top)
            ax.bar(bins_angles, occs_dict[region]['neurites']['axons'],
                    bottom = occs_dict[region]['neurites']['dendrites'],
                    width = resul_binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = neurite_color_dict[region]['axons'],
                    zorder = 2)
            
        else:
            ax.bar(bins_angles, occu,
                    width = resul_binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = neurite_color_dict[region][neurite_type],
                    zorder = 2)  
           

### y axis

# neurites
for ax in axs_norm[:9:3]:
    ax.set_yticks(ticks = np.arange(0, 0.30 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '', '30 %'], zorder = 3)
    ax.set_rlabel_position(80)

# dendrites
for ax in axs_norm[1:9:3]:
    ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '10 %', '', '20 %', ''], zorder = 3)
    ax.set_rlabel_position(80)

# axons
for ax in axs_norm[2:9:3]:
    ax.set_yticks(ticks = np.arange(0, 0.35 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '', '30 %', ''], va = 'top', zorder = 3)
    ax.set_rlabel_position(-80)


### x axis
for ax in axs_norm[:9]:
    # x axis
    ax.set_xticks(ticks = np.arange(0, np.pi*2, np.pi / 4), labels = orientation_labels, zorder = 3)

    # grid
    ax.grid(True, alpha = 0.5, color = 'gray', zorder = 0)
    ax.set_axisbelow(True)     
        

### number of terminal points ###

# set neurite_type lookup dict for terminals
neurite_terminal_types = {'neurites' : 'neuritic_terminals',
                          'dendrites' : 'dendritic_terminals',
                          'axons' : 'axonic_terminals'}

for neurite_idx, neurite_type in enumerate(neurite_types):
    
    # set axis
    ax = axs_norm[neurite_idx + 9]  
    
    # set axis title
    ax.set_title(alpha_labels_dict['n_terminals'][neurite_type] + ': ' + neurite_type.title(),
                  fontsize=12, 
                  loc='left',
                  x = -0.4)
    
    # set violin data
    # get data for both regions
    violin_data_bothregions = n_terminal_points[n_terminal_points['neurite_type'] == neurite_terminal_types[neurite_type]]
    violin_data_bothregions.loc[:, 'Region'] = ['both_regions'] * violin_data_bothregions.shape[0]

    # get data for individual regions
    violin_data = n_terminal_points[n_terminal_points['neurite_type'] == neurite_terminal_types[neurite_type]]
    violin_data = violin_data[violin_data['Region'] != 'BAOT/MeA']
    
    # concat data
    violin_data = pd.concat([violin_data, violin_data_bothregions], axis = 0)
      
    # violin plots   
    violins = sbn.violinplot(data = violin_data,
                              x = 'Region', 
                              y = 'n_terminal_points',
                              hue = True, hue_order=[True, False], split = True,
                              ax = ax, 
                              scale='width',
                              inner='quart',
                              linewidth = 1,
                              order = ['both_regions', 'MeA', 'BAOT'])
    
    # edit lines of quarts
    for r_idx, region in enumerate(['all', 'MeA', 'BAOT']):
        
        for l_idx in np.arange(0, 3):
        
            all_l_idx = r_idx * 3 + l_idx
            
            # violins.lines[all_l_idx].set_color(neurite_color_dict[region][neurite_type])
            violins.lines[all_l_idx].set_color(colors_dict['primecolor'])
            violins.lines[all_l_idx].set_linestyle('solid')
        
        # set edge color of violin
        violins.collections[r_idx].set_edgecolor(neurite_color_dict[region][neurite_type])
        # violins.collections[r_idx].set_edgecolor('None')
        
        # set facecolor of violin
        # violins.collections[r_idx].set_facecolor('None')
        violins.collections[r_idx].set_facecolor(neurite_color_dict[region][neurite_type])
    
    
    # swarmplots    
    sbn.swarmplot(data = violin_data,
                  x = 'Region', 
                  y = 'n_terminal_points',
                  ax = ax,
                  s = 1.5,
                  color=colors_dict['primecolor'],
                  order = ['both_regions', 'MeA', 'BAOT'])
        
    # errorbar
    for r_idx, region in enumerate(['all', 'MeA', 'BAOT']):
        ax.errorbar(x = r_idx+0.3,
                    y = n_terminal_points_means.at[region, neurite_terminal_types[neurite_type]],
                    yerr = n_terminal_points_stds.at[region, neurite_terminal_types[neurite_type]],
                    fmt='_', 
                    markersize = 6,
                    markerfacecolor = 'none',
                    capsize = 2,
                    color=neurite_color_dict[region][neurite_type],
                    linewidth = 1,
                    label = '_nolegend_')
    
    
    # edit axes
    # x
    ax.set_xlim(0-0.5, 2+0.5)
    ax.set_xticks(ticks=np.arange(0, 2 + 0.5, 1),
                  labels=['Both\nregions', 'MeA', 'BAOT'], rotation=60)
    ax.spines['bottom'].set_bounds([0, 2])
    ax.set_xlabel('')
    
    # edit seaborn legend
    ax.legend().set_visible(False)
    
    
    # y
    if neurite_type == 'neurites' or neurite_type == 'dendrites':
        ymin = 0
        ymax = 45
        ystep = 10
        ystepminor = 5
    
    elif neurite_type == 'axons':
        ymin = 0
        ymax = 10
        ystep = 5
        ystepminor = 1
    
    
    ypad = (ymax - ymin) * 0.05
    ax.set_ylim(ymin-ypad, ymax+ypad)
    ax.set_yticks(ticks=np.arange(ymin, ymax+1, ystep))
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystepminor, ystepminor), minor=True)
    ax.spines['left'].set_bounds([ymin, ymax])
    
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # axis labels
    ax.set_ylabel('Number of terminal\nbranches [#]')
    ax.set_xlabel('Region')

# align labels
fig_norm.align_labels()

# show figure
plt.show()

# save figure
save_figures(fig_norm, figure_name = 'population_polar_plots-new_orientation-figure', save_dir = join(cell_morph_plots_dir, 'polar_plots_population'),
              darkmode_bool = darkmode_bool, figure_format= 'both')






# %% exports occus to excel sheet

# for occu_df, region_label in zip([ALL_pp_normed_totype_mean_occs, MeA_pp_normed_totype_mean_occs, BAOT_pp_normed_totype_mean_occs], ['All', 'MeA', 'BAOT']):

#     occu_df.to_excel(join(cell_morph_plots_dir, 
#                           'polar_plots_population', 
#                           'population_polar_plots-normed_totype-regions-' + region_label + '-means.xlsx'), 
#                      index = orientation_labels)
    
    
# for occu_df, region_label in zip([ALL_pp_normed_toneurites_mean_occs, MeA_pp_normed_toneurites_mean_occs, BAOT_pp_normed_toneurites_mean_occs], ['All', 'MeA', 'BAOT']):

#     occu_df.to_excel(join(cell_morph_plots_dir, 
#                           'polar_plots_population', 
#                           'population_polar_plots-normed_toneurite-regions-' + region_label + '-means.xlsx'), 
#                      index = orientation_labels)    
    
