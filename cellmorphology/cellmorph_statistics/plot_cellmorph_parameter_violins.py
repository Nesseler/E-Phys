# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 11:22:10 2024

@author: nesseler
"""

from cellmorphology.cellmorph_packages import mtl, plt, sbn, pd, np, join

from parameters.directories_win import table_file, cell_morph_descrip_dir, cell_morph_plots_dir

from cellmorphology.cellmorph_parameters import cell_IDs_toDrop

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# set list of neurite types
neurite_types = ['neurites', 'dendrites', 'axons']

# %% load & clean data

n_primary = pd.read_excel(join(cell_morph_descrip_dir, 'n_primary_points.xlsx'), index_col = 'cell_ID')
n_terminal = pd.read_excel(join(cell_morph_descrip_dir, 'n_terminal_points.xlsx'), index_col = 'cell_ID')
bifurcation_ratios = pd.read_excel(join(cell_morph_descrip_dir, 'bifurcation_ratios.xlsx'), index_col = 'cell_ID')
total_cable_length = pd.read_excel(join(cell_morph_descrip_dir, 'total_cable_length.xlsx'), index_col = 'cell_ID')

# write parameters to one dataframec

cellmorph_parameters = pd.concat([n_primary, n_terminal, bifurcation_ratios, total_cable_length], axis = 1)

# rename columns to unify naming scheme
for neurite_prefix, neurite_type in zip(['neuritic', 'dendritic', 'axonic'], neurite_types):

    cellmorph_parameters.rename(columns = {f'{neurite_prefix}_primaries'         : f'{neurite_type}-n_primaries', 
                                           f'{neurite_prefix}_terminals'         : f'{neurite_type}-n_terminals',
                                           f'bifurcation_ratio_{neurite_prefix}' : f'{neurite_type}-bifurcation_ratio',
                                           f'total_cable_length-{neurite_type}'  : f'{neurite_type}-total_cable_length'},
                                           inplace = True)

# remove cells that do not contain all analysed values
cellmorph_parameters.drop(index = cell_IDs_toDrop, inplace = True)

# remove axonless cells from 
for col in cellmorph_parameters.columns:
    if 'axon' in col:
        cellmorph_parameters.replace(to_replace = {col : 0},
                                     value = np.nan,
                                     inplace = True)
        

# set list of parameters
parameters = ['n_primaries', 'n_terminals', 'bifurcation_ratio', 'total_cable_length']

# add region
cellmorph_parameters['Region'] = MetaData.loc[cellmorph_parameters.index, 'Region']


# %% reorder dataframe for plotting

# create dataframe for plotting
cellmorph_parameters_plotting = pd.DataFrame()

for p_idx, parameter in enumerate(parameters):
    
    # create interim dataframe
    cellmorph_parameters_plotting_perParameter = pd.DataFrame()
    
    for neurite_type in neurite_types:      
        
        # define interim dataframe to fill
        cellmorph_parameters_plotting_perType = pd.DataFrame()
        
        # write to Dataframe
        cellmorph_parameters_plotting_perType[parameter] = cellmorph_parameters[f'{neurite_type}-{parameter}']
        cellmorph_parameters_plotting_perType['neurite_type'] = [neurite_type] * len(cellmorph_parameters[f'{neurite_type}-{parameter}'])
        
        # concat to perParameter
        cellmorph_parameters_plotting_perParameter = pd.concat([cellmorph_parameters_plotting_perParameter, cellmorph_parameters_plotting_perType], axis = 0)
    
    # concat to plotting
    cellmorph_parameters_plotting = pd.concat([cellmorph_parameters_plotting, cellmorph_parameters_plotting_perParameter], axis = 1)
    
# remove duplicate columns
cellmorph_parameters_plotting = cellmorph_parameters_plotting.loc[:,~cellmorph_parameters_plotting.columns.duplicated()].copy()

# add region
cellmorph_parameters_plotting['Region'] = MetaData.loc[cellmorph_parameters_plotting.index, 'Region']

# remove cells without region
cellmorph_parameters_plotting = cellmorph_parameters_plotting.query('Region != "BAOT/MeA"')


# %% calc mean, median, std

cellmorph_parameters_mean_median_std = pd.DataFrame()
    
for region in ['MeA', 'BAOT']:
    
    # set dataframe
    cellmorph_parameters_perRegion = cellmorph_parameters[cellmorph_parameters['Region'] == region]
    
    cellmorph_parameters_mean_median_std[f'{region}-mean'] = cellmorph_parameters_perRegion.mean(numeric_only = True)
    cellmorph_parameters_mean_median_std[f'{region}-median'] = cellmorph_parameters_perRegion.median(numeric_only = True)
    cellmorph_parameters_mean_median_std[f'{region}-std'] = cellmorph_parameters_perRegion.std(numeric_only = True)


# %% initialize plotting

from functions.functions_plotting import save_figures, get_colors, get_figure_size
from cellmorphology.cellmorph_colors import neurite_color_dict

# set colors
darkmode_bool = False
colors_dict, region_colors = get_colors(darkmode_bool)

# set font size
mtl.rcParams.update({'font.size': 9})


# %% plotting

# set subplot titles
subplot_labels = {'n_primaries' :       {'neurites'  : '$\mathregular{A_{ii}}$: Neurites',
                                         'dendrites' : '$\mathregular{A_{iii}}$: Dendrites',
                                         'axons'     : '$\mathregular{A_{iv}}$: Axons'},
                  'n_terminals' :       {'neurites'  : '$\mathregular{B_{ii}}$: Neurites',
                                         'dendrites' : '$\mathregular{B_{iii}}$: Dendrites',
                                         'axons'     : '$\mathregular{B_{iv}}$: Axons'},
                  'bifurcation_ratio' : {'neurites'  : '$\mathregular{C_{ii}}$: Neurites',
                                         'dendrites' : '$\mathregular{C_{iii}}$: Dendrites',
                                         'axons'     : '$\mathregular{C_{iv}}$: Axons'},
                  'total_cable_length': {'neurites'  : '$\mathregular{D_{ii}}$: Neurites',
                                         'dendrites' : '$\mathregular{D_{iii}}$: Dendrites',
                                         'axons'     : '$\mathregular{D_{iv}}$: Axons'}}


for parameter in parameters:

    # create figure
    fig, axs = plt.subplots(nrows = 1, 
                            ncols = 3,
                            layout = 'constrained',
                            figsize = get_figure_size(width = 103.8, height = 58),
                            dpi = 600)
    
    # iterate through neurite_types
    for n_idx, neurite_type in enumerate(neurite_types):
        
        # set axis
        ax = axs[n_idx]
        
        # set title
        ax.set_title(subplot_labels[parameter][neurite_type], 
                     fontsize = 12,
                     loc = 'left')
        
        # limit data to perType
        cellmorph_parameters_plotting_perType = cellmorph_parameters_plotting[cellmorph_parameters_plotting['neurite_type'] == neurite_type]
        
        # set neurites color dict
        neuriteType_color_dict = {region : neurite_color_dict[region][neurite_type] for region in ['MeA', 'BAOT']}
        
        # violinplot
        violins = sbn.violinplot(data = cellmorph_parameters_plotting_perType,
                                 y = parameter,
                                 hue = 'Region',
                                 hue_order = ['MeA', 'BAOT'],
                                 palette = neuriteType_color_dict,
                                 split = True,
                                 gap = 0.1,
                                 linewidth = 1,
                                 density_norm = 'width',
                                 inner = None,
                                 ax = ax, 
                                 zorder = 1)
    
        # edit violins
        if not (neurite_type == 'axons' and parameter == 'n_primaries'):        
            for r_idx, region in enumerate(['MeA', 'BAOT']):
                violins.collections[r_idx].set_edgecolor(neuriteType_color_dict[region])
                violins.collections[r_idx].set_facecolor(neuriteType_color_dict[region])
    
        # swarmplots
        swarm = sbn.swarmplot(data = cellmorph_parameters_plotting_perType,
                              y = parameter,
                              hue = 'Region',
                              hue_order = ['MeA', 'BAOT'],
                              palette = ['k', 'k'],
                              ax=ax,
                              size=1.5,
                              dodge=True,
                              zorder = 2)
        
        # edit seaborn legend
        ax.legend().set_visible(False)
        
        # set x offset
        x_offset = [-0.1, +0.1]
        
        # plot mean, median, std
        for r_idx, region in enumerate(['MeA', 'BAOT']):
        
            # plot errorbar
            ax.errorbar(x = 0 + x_offset[r_idx],
                        y = cellmorph_parameters_mean_median_std.at[f'{neurite_type}-{parameter}', f'{region}-mean'],
                        yerr = cellmorph_parameters_mean_median_std.at[f'{neurite_type}-{parameter}', f'{region}-std'],
                        fmt='_', 
                        markersize = 6,
                        markerfacecolor = 'none',
                        capsize = 3,
                        color= colors_dict['primecolor'],
                        linewidth = 1,
                        elinewidth = 1,
                        label = '_nolegend_',
                        zorder = 3)
                
            # plot median
            ax.scatter(x = 0 + x_offset[r_idx],
                       y = cellmorph_parameters_mean_median_std.at[f'{neurite_type}-{parameter}', f'{region}-median'],
                       marker='D', 
                       s = 5,
                       color=colors_dict['primecolor'],
                       linewidth = 1,
                       label = '_nolegend_',
                       zorder = 4)    
        
        # x
        ax.set_xlim(0-0.5, 0+0.5)
        ax.set_xticks(ticks=np.arange(-0.2, 0.2+0.1, 0.4), labels=[region for region in ['MeA', 'BAOT']], rotation = 25)
        ax.spines['bottom'].set_bounds([-0.2, 0.2])
        axs[1].set_xlabel('Region')
    
        # y
        if parameter == 'n_primaries':
            ylabel = 'Number of primary\nbranches per cell [#]'
                
            if n_idx == 2:
                ymin = 0
                ymax = 2
                ystep = 1
                ystepminor = 1
            else:
                ymin = 0
                ymax = 10
                ystep = 5
                ystepminor = 1
    
        elif parameter == 'n_terminals':
            ylabel = 'Number of terminal\nbranches per cell [#]'
    
            if n_idx == 2:
                ymin = 0
                ymax = 12
                ystep = 5
                ystepminor = 1
            else:
                ymin = 0
                ymax = 46
                ystep = 20
                ystepminor = 2
                
        elif parameter == 'bifurcation_ratio':
            ylabel = 'Bifurcation ratio\nper cell [#]'
                
            ymin = 0
            ymax = 20
            ystep = 10
            ystepminor = 1
            
        elif parameter == 'total_cable_length':
            ylabel = 'Total cable length\nper cell [µm]'
                
            if n_idx == 2:
                ymin = 0
                ymax = 3000
                ystep = 1000
                ystepminor = 250
            else:
                ymin = 0
                ymax = 6500
                ystep = 2000
                ystepminor = 500
            
        # ylabel
        if n_idx == 0:
            ax.set_ylabel(ylabel)
        else:
            ax.set_ylabel('')
            
        # define ypad relative to range
        ypad = (ymax - ymin) * 0.05
    
        # apply y axis settings
        ax.set_ylim(ymin - ypad, ymax + ypad)
        ax.set_yticks(ticks=np.arange(ymin, ymax + ystepminor, ystep))
        ax.set_yticks(ticks=np.arange(ymin, ymax + ystepminor, ystepminor), minor=True)
        ax.spines['left'].set_bounds([ymin, ymax])
    
        [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # show plot
    plt.show()
    
    # save figure
    figure_dir = join(cell_morph_plots_dir, 'figure-n_primaries-n_terminals-bifurcation_ratio-total_cable_length')
    
    save_figures(fig,
                  figure_name = f'figure-{parameter}', 
                  save_dir = figure_dir,
                  darkmode_bool=darkmode_bool, 
                  figure_format='both')
    

# %% statistics

from scipy.stats import normaltest, mannwhitneyu, ttest_ind, ks_1samp
from scipy import stats

# create dataframe
cellmorph_normaltest = pd.DataFrame()
    
for region in ['MeA', 'BAOT']:
    
    # set dataframe
    cellmorph_parameters_perRegion = cellmorph_parameters[cellmorph_parameters['Region'] == region]
    
    for col in cellmorph_parameters.drop(columns = ['Region']).columns:
    
        # run rest
        kstest_res = ks_1samp(cellmorph_parameters_perRegion[col],
                              stats.norm.cdf,
                              nan_policy='omit')
        
        # write to dataframe
        cellmorph_normaltest.at[f'{col}-{region}', 'kstest_statistic'] = kstest_res.statistic
        cellmorph_normaltest.at[f'{col}-{region}', 'kstest_pvalue'] = kstest_res.pvalue
        
        if kstest_res.pvalue > 0.05 :  
            cellmorph_normaltest.at[f'{col}-{region}', 'kstest-normally_distributed'] = True
        else:
            cellmorph_normaltest.at[f'{col}-{region}', 'kstest-normally_distributed'] = False

    
# create dataframe for each region
cellmorph_parameters_MeA = cellmorph_parameters[cellmorph_parameters['Region'] == 'MeA'].drop(columns = ['Region'])
cellmorph_parameters_BAOT = cellmorph_parameters[cellmorph_parameters['Region'] == 'BAOT'].drop(columns = ['Region'])
    
# create dataframe for statistics
cellmorph_parameters_between_regions = pd.DataFrame(index = cellmorph_parameters_MeA.columns,
                                                    columns = ['mannwhitneyu_stats', 'mannwhitneyu_pvalue', 'ttest_ind_stats', 'ttest_ind_pvalue'])

for col in cellmorph_parameters_MeA:
    
    # look for normal distribution in both distributions
    MeA_normal_pvalue = cellmorph_normaltest.at[f'{col}-MeA', 'kstest_pvalue']
    BAOT_normal_pvalue = cellmorph_normaltest.at[f'{col}-BAOT', 'kstest_pvalue']
    
    if (MeA_normal_pvalue > 0.05 and BAOT_normal_pvalue > 0.05):
        # run rest
        ttest_res = ttest_ind(cellmorph_parameters_perRegion[col], nan_policy='omit', alternative='two-sided')
        
        # write to dataframe
        cellmorph_parameters_between_regions.at[f'{col}', 'ttest_ind_stats'] = ttest_res.statistic
        cellmorph_parameters_between_regions.at[f'{col}', 'ttest_ind_pvalue'] = ttest_res.pvalue

    else:
        # run rest
        mannwhitneyu_res = mannwhitneyu(cellmorph_parameters_MeA[col], cellmorph_parameters_BAOT[col],
                                        alternative='two-sided', nan_policy='omit')
        
        # write to dataframe
        cellmorph_parameters_between_regions.at[f'{col}', 'mannwhitneyu_stats'] = mannwhitneyu_res.statistic
        cellmorph_parameters_between_regions.at[f'{col}', 'mannwhitneyu_pvalue'] = mannwhitneyu_res.pvalue
        
        # write boolean
        if mannwhitneyu_res.pvalue < 0.05:
            cellmorph_parameters_between_regions.at[f'{col}', 'mannwhitneyu-statistical_difference'] = True
        else:
            cellmorph_parameters_between_regions.at[f'{col}', 'mannwhitneyu-statistical_difference'] = False
            
            
# Bonferroni correction

from statsmodels.stats.multitest import multipletests

# get rows to apply correction to
for parameter in parameters:

    # get rows
    parameter_rows = [neurite_type + "-" + parameter for neurite_type in neurite_types]
    
    # get pvalues of mannwhitneyu test
    parameter_pvalues = cellmorph_parameters_between_regions.loc[parameter_rows, 'mannwhitneyu_pvalue'].to_list()
    
    # apply correction
    rejects, pvals_corrected, _, _ = multipletests(parameter_pvalues,
                                                   alpha=0.05, 
                                                   method='bonferroni')
    
    # write to dataframe
    cellmorph_parameters_between_regions.loc[parameter_rows, 'bonferroni_pvalue'] = pvals_corrected
    
    # write boolean
    cellmorph_parameters_between_regions.loc[parameter_rows, 'bonferroni-statistical_difference'] = cellmorph_parameters_between_regions.loc[parameter_rows, 'bonferroni_pvalue'] < 0.05
    


cellmorph_normaltest.to_excel(join(figure_dir, 'normaltest.xlsx'), index_label= 'neurite_type-parameter-region')
cellmorph_parameters_between_regions.to_excel(join(figure_dir, 'statistical_difference.xlsx'), index_label= 'neurite_type-parameter')
 

