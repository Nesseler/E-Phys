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


# %% drop cells E-126 & E-158

for neurite_type in neurite_types:
    
    polar_occu_dict[neurite_type].drop(columns = ['E-126', 'E-158'], inplace = True)


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

# set neurite_type lookup dict for terminals
neurite_terminal_types = {'neurites' : 'neuritic_terminals',
                          'dendrites' : 'dendritic_terminals',
                          'axons' : 'axonic_terminals'}

# load excel sheet
n_terminal_points_df = pd.read_excel(join(cell_morph_descrip_dir, 'n_terminal_points.xlsx'), index_col = 'cell_ID')

# set number of axonic terminal points from 0 to nan
n_terminal_points_df.replace(to_replace = 0,
                             value = np.nan,
                             inplace = True)

# melt df to plot
n_terminal_points = pd.melt(n_terminal_points_df, var_name= 'neurite_type', value_name= 'n_terminal_points', ignore_index=False)

# get Regions
n_terminal_points['Region'] = MetaData.loc[n_terminal_points_df.index, 'Region']



# calc mean and std of number of terminal points
n_terminal_points_means = pd.DataFrame(columns = ['neuritic_terminals', 'dendritic_terminals', 'axonic_terminals'],
                                       index = ['MeA', 'BAOT'])
n_terminal_points_median = pd.DataFrame(columns = ['neuritic_terminals', 'dendritic_terminals', 'axonic_terminals'],
                                        index = ['MeA', 'BAOT'])
n_terminal_points_stds = pd.DataFrame(columns = ['neuritic_terminals', 'dendritic_terminals', 'axonic_terminals'],
                                      index = ['MeA', 'BAOT'])

# get cell_IDs in regions
# all_cell_IDs = n_terminal_points_df.index.to_list()
# region_cell_IDs = {'all' : all_cell_IDs,
#                     'MeA' : MetaData.loc[all_cell_IDs, 'Region'][MetaData.loc[all_cell_IDs, 'Region'] == 'MeA'].index.to_list(),
#                     'BAOT': MetaData.loc[all_cell_IDs, 'Region'][MetaData.loc[all_cell_IDs, 'Region'] == 'BAOT'].index.to_list()}

# calc per region
for region in regions:
    
    for neurite_type_terminals in neurite_terminal_types.values():
        
        # get current cell_IDs
        if neurite_type_terminals != 'axons':
            cur_region_cell_IDs = region_cell_IDs[region]
        elif neurite_type_terminals == 'axons':
            cur_region_cell_IDs = region_cell_IDs_w_axon[region]
            
        n_terminal_points_means.at[region, neurite_type_terminals] = n_terminal_points_df.loc[cur_region_cell_IDs, neurite_type_terminals].mean()
        n_terminal_points_stds.at[region, neurite_type_terminals] = n_terminal_points_df.loc[cur_region_cell_IDs, neurite_type_terminals].std()
        n_terminal_points_median.at[region, neurite_type_terminals] = n_terminal_points_df.loc[cur_region_cell_IDs, neurite_type_terminals].median()


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
                                  figsize = get_figure_size(width = 160, height = 175),
                                  dpi = 600,
                                  height_ratios=[0.2, 0.2, 0.3])

# flatten axes array
axs_norm = axs_norm.flatten()

### polar plots ###

# change projection for polar plots
for ax in axs_norm[:6]:
    change_projection(fig_norm, axs_norm, ax, projection = 'polar')
    
# set polar_plot idx
polar_idx = 0

# loop through regions (rows)
for region in regions:
 
    # loop through neurite types (cols)
    for neurite_type in neurite_types:
 
        # set axis
        ax = axs_norm[polar_idx]
        
        # iterate on polar index
        polar_idx += 1
        
        # set axis title
        ax.set_title(alpha_labels_dict[region][neurite_type] + ': ' + region + ' ' + neurite_type,
                      fontsize=12, 
                      loc='left',
                      x = -0.4)
    

        
        # set indices for dataframe
        polar_occu_cols = [f'{neurite_type_ls}-{region}' for neurite_type_ls in neurite_types]
        
        # plot differently for neurites and dendrites / axons (stack bar plot for neurites)
        if neurite_type == 'neurites':
            
            # dendrites normed to all neurites (bottom)
            ax.bar(bins_angles, polar_occu_toNeurites[f'dendrites-{region}'],
                    width = resul_binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = neurite_color_dict[region]['dendrites'],
                    zorder = 2)
            
            # axons normed to all neurites (top)
            ax.bar(bins_angles, polar_occu_toNeurites[f'axons-{region}'],
                    bottom = polar_occu_toNeurites[f'dendrites-{region}'],
                    width = resul_binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = neurite_color_dict[region]['axons'],
                    zorder = 2)
            
        else:
            ax.bar(bins_angles, polar_occu_toType[f'{neurite_type}-{region}'],
                    width = resul_binsize, 
                    align = 'edge',
                    edgecolor = 'none',
                    color = neurite_color_dict[region][neurite_type],
                    zorder = 2)  
           

### y axis

# neurites
for ax in axs_norm[:6:3]:
    ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '25 %'], zorder = 3)
    ax.set_rlabel_position(80)

# dendrites
for ax in axs_norm[1:6:3]:
    ax.set_yticks(ticks = np.arange(0, 0.25 + 0.01, 0.05), labels = ['', '', '', '15 %', '', '25 %'], zorder = 3)
    ax.set_rlabel_position(80)

# axons
for ax in axs_norm[2:6:3]:
    ax.set_yticks(ticks = np.arange(0, 0.40 + 0.01, 0.1), labels = ['', '', '20 %', '', '40 %'], va = 'top', zorder = 3)
    ax.set_rlabel_position(-80)


### x axis
for ax in axs_norm[:9]:
    # x axis
    ax.set_xticks(ticks = np.arange(0, np.pi*2, np.pi / 4), labels = orientation_labels, zorder = 3)

    # grid
    ax.grid(True, alpha = 0.5, color = 'gray', zorder = 0)
    ax.set_axisbelow(True)     
        

### number of terminal points ###



for neurite_idx, neurite_type in enumerate(neurite_types):
    
    # set axis
    ax = axs_norm[neurite_idx + 6]  
    
    # set axis title
    ax.set_title(alpha_labels_dict['n_terminals'][neurite_type] + ': ' + neurite_type.title(),
                 fontsize=12, 
                 loc='left',
                 x = -0.4)
    
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
                             hue_order = regions,
                             split = True,
                             gap = .1,
                             ax = ax, 
                             inner = None,
                             zorder = 0,
                             linewidth = 1,
                             density_norm = 'width')
    
    # edit lines of quarts
    for r_idx, region in enumerate(regions):
        
        for l_idx in np.arange(0, 3):
        
            all_l_idx = r_idx * 3 + l_idx
            
            # violins.lines[all_l_idx].set_color(neurite_color_dict[region][neurite_type])
            # violins.lines[all_l_idx].set_color(colors_dict['primecolor'])
            # violins.lines[all_l_idx].set_linestyle('solid')
        
        # set edge color of violin
        violins.collections[r_idx].set_edgecolor(neurite_color_dict[region][neurite_type])
        # violins.collections[r_idx].set_edgecolor('None')
        
        # set facecolor of violin
        # violins.collections[r_idx].set_facecolor('None')
        violins.collections[r_idx].set_facecolor(neurite_color_dict[region][neurite_type])
    
    
    # swarmplots    
    sbn.swarmplot(data = violin_data,
                  x = 'X', 
                  y = 'n_terminal_points',
                  hue = 'Region', 
                  hue_order = regions,
                  dodge = True,
                  ax = ax,
                  s = 2,
                  color=colors_dict['primecolor'],
                  zorder = 1)
        
    # errorbar
    for region_x, region in zip([-0.1, 0.1], regions):
        
        # get mean and std
        mean = n_terminal_points_means.at[region, neurite_terminal_types[neurite_type]]
        median = n_terminal_points_median.at[region, neurite_terminal_types[neurite_type]]
        std = n_terminal_points_stds.at[region, neurite_terminal_types[neurite_type]]
        
        ax.errorbar(x = region_x,
                    y = mean,
                    yerr = std,
                    fmt='_', 
                    markersize = 7,
                    markerfacecolor = 'none',
                    capsize = 3,
                    color=colors_dict['primecolor'],
                    linewidth = 1,
                    label = '_nolegend_',
                    zorder = 2)
        
        ax.scatter(x = region_x,
                   y = median,
                   marker='D', 
                   s = 5,
                   color=colors_dict['primecolor'],
                   linewidth = 1,
                   label = '_nolegend_',
                   zorder = 3)
    
    
    # edit axes
    # x
    xshift = 0.2
    
    ax.set_xlim(0-0.5, 0+0.5)
    ax.set_xticks(ticks=np.arange(-xshift, xshift+0.1, xshift*2),
                  labels=['MeA', 'BAOT'], rotation=0)
    ax.spines['bottom'].set_bounds([-xshift, xshift])
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
        ymax = 15
        ystep = 5
        ystepminor = 1
    
    
    ypad = (ymax - ymin) * 0.05
    ax.set_ylim(ymin-ypad, ymax+ypad)
    ax.set_yticks(ticks=np.arange(ymin, ymax+1, ystep))
    ax.set_yticks(ticks=np.arange(ymin, ymax + ystepminor, ystepminor), minor=True)
    ax.spines['left'].set_bounds([ymin, ymax])
    
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # axis label 
    ax.set_xlabel('')
    ax.set_ylabel('')
    
    # set grid
    ax.grid(False)

    
# axis labels
axs_norm[6].set_ylabel('Number of terminal\nbranches per cell [#]')

# align labels
fig_norm.align_labels()

# show figure
plt.show()

# save figure
figure_dir = join(cell_morph_plots_dir, 'figure-polar_population')
save_figures(fig_norm, 
             figure_name = 'population_polar_plots-new_figure', 
             save_dir = figure_dir,
             darkmode_bool = darkmode_bool, 
             figure_format= 'both')


# %% save dataframes of data in figure

# save dataframes
polar_occu_toNeurites.to_excel(join(figure_dir, 'polar_occurrances-population_fig-normedToNeurites.xlsx'), index_label='orientation_rad')
polar_occu_toType.to_excel(join(figure_dir, 'polar_occurrances-population_fig-normedToType.xlsx'), index_label='orientation_rad')

n_terminal_points.to_excel(join(figure_dir, 'polar_occurrances-population_fig-n_terminals.xlsx'), index_label='cell_ID')

# merge means, medians and stds into one dataframe
mean_median_std_terminal_points = pd.DataFrame(columns = ['mean', 'median', 'std'])

for label, measurement in zip(['mean', 'median', 'std'], [n_terminal_points_means, n_terminal_points_median, n_terminal_points_stds]):
    
    for region in regions:
        
        for neurite_type_terminals in neurite_terminal_types.values():
            
            mean_median_std_terminal_points.at[f'{neurite_type_terminals}-{region}', label] = measurement.at[region, neurite_type_terminals]

mean_median_std_terminal_points.to_excel(join(figure_dir, 'polar_occurrances-population_fig-n_terminals-means_medians_stds.xlsx'), index_label='neurite_type-region')


# %% statistical tests

from scipy.stats import normaltest, mannwhitneyu

# test for normal distribution per neurite type and region

# create dataframe for p_values
n_terminals_normaltest_pvalues = pd.DataFrame()

for neurite_type_terminals in neurite_terminal_types.values():

    for region in regions:
        
        n_terminal_points_per_type_n_per_region = n_terminal_points_df.loc[region_cell_IDs[region], neurite_type_terminals]

        ntest_stats, ntest_pvalue = normaltest(n_terminal_points_per_type_n_per_region, nan_policy= 'omit')
        
        n_terminals_normaltest_pvalues.at[f'{neurite_type_terminals}-{region}', 'statistic'] = ntest_stats
        n_terminals_normaltest_pvalues.at[f'{neurite_type_terminals}-{region}', 'pvalue'] = ntest_pvalue
        
# save normaltest p_values
n_terminals_normaltest_pvalues.to_excel(join(figure_dir, 'population_fig-n_terminals-ntest_pvalues.xlsx'), index_label='neurite_type-region')


# test for statistical difference between regions

# create dataframe
n_terminals_mannwhitneyu_pvalues = pd.DataFrame()

for neurite_type_terminals in neurite_terminal_types.values():
     
    # get n_terminal_points
    n_terminal_points_per_type_MeA = n_terminal_points_df.loc[region_cell_IDs['MeA'], neurite_type_terminals]
    n_terminal_points_per_type_BAOT = n_terminal_points_df.loc[region_cell_IDs['BAOT'], neurite_type_terminals]

    # apply test
    mannwhitneyu_stats, mannwhitneyu_pvalue = mannwhitneyu(x = n_terminal_points_per_type_MeA,
                                                           y = n_terminal_points_per_type_BAOT,
                                                           nan_policy= 'omit')
        
    # write to dataframe
    n_terminals_mannwhitneyu_pvalues.at[f'{neurite_type_terminals}', 'mannwhitneyu_statistic'] = mannwhitneyu_stats
    n_terminals_mannwhitneyu_pvalues.at[f'{neurite_type_terminals}', 'mannwhitneyu_pvalue'] = mannwhitneyu_pvalue


# save mannwhitneyu p_values
n_terminals_mannwhitneyu_pvalues.to_excel(join(figure_dir, 'population_fig-n_terminals-mannwhitneyu_pvalues.xlsx'), index_label='neurite_type')
