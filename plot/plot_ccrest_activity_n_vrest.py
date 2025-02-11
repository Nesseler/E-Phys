# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 18:12:32 2024

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir, cell_descrip_syn_dir, figure_dir, table_file

# init plotting
from functions.initialize_plotting import *

# exp = '-' 
exp = '-Syn-' 

activity_df = pd.read_excel(join(cell_descrip_syn_dir, 'cc_rest' + exp + 'activity.xlsx'), index_col = 'cell_ID')
activity_df2 = pd.read_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col = 'cell_ID')

activity_df = pd.concat([activity_df, activity_df2])

cell_IDs = list(activity_df.index)

regions = ['BAOT/MeA', 'MeA', 'BAOT']

# %% load Metadata

MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

activity_df = pd.concat([activity_df, MetaData.loc[:, ['Region', 'Sex']]], axis = 1)

# %%

# rewrite region column as categorical data with order
activity_df['Region'] = pd.Categorical(activity_df['Region'], categories = regions, ordered = True)


# sort dataframe
activity_df = activity_df.sort_values('n_spikes', ascending = True)

cell_IDs = list(activity_df.index)
n_cells = len(cell_IDs)

tick_size = 0.9



# %% initialize plotting

from functions.initialize_plotting import *



# %% EVENTPLOT + V_REST FIGURE in Regions

# get n cells per region and calc percentage
cells_of_region_perc = [activity_df['Region'].value_counts()[region] / n_cells for region in regions]


# initialize figure
ax_keys = ['A', 'B', 'C', 'D']

fig_regions, axs_regions = plt.subplot_mosaic('A;B;C;D', 
                                               layout = 'constrained',
                                              figsize = get_figure_size(width = 154.335, height = 150.5),
                                              # width_ratios = [2, 1.5],
                                              dpi = 600,
                                              height_ratios = cells_of_region_perc + [2.5])

fig_regions.set_constrained_layout_pads(w_pad=50./72., 
                                        h_pad=4./72.,
                                        hspace=0./72.,
                                        wspace=9./72.)
    
# eventplot

for idx_region, region in enumerate(regions):
    
    # set axis
    ax = axs_regions[ax_keys[idx_region]]
    
    # get cell_IDs of dataframe per region
    cell_IDs_region = activity_df[activity_df['Region'] == region].index.to_list()
    
    # get number of cells in region
    n_cells_region = len(cell_IDs_region)
    
    for idx_cell, cell_ID in enumerate(cell_IDs_region):
        
        # read time points of spikes as string representation of a list
        t_spikes = activity_df.at[cell_ID, 't_spikes'].strip('][').split(', ')
        
        # convert individual string elements to floats
        # check for empty list not to be converted to float
        if len(t_spikes) > 1:
            t_spikes = [float(t) for t in t_spikes]
        else:
            t_spikes = []
        
        ax.eventplot(t_spikes,
                     orientation = 'horizontal', 
                     lineoffsets=idx_cell, 
                     linewidth = 1.5,
                     linelengths=0.9, 
                     color = region_colors[MetaData.at[cell_ID, 'Region']])
    
    # set title
    # ax.set_title(region, fontsize = 14)
    
    # set shared x axis
    ax.set_xlim([0-0.5, 30+0.5])
    ax.set_xticks(np.arange(0, 30+1, 10))
    ax.set_xticks(np.arange(0, 30+1, 1), minor = True)
    
    # set y axis for each subplot  
    ax.set_ylim([0-(tick_size/2), n_cells_region-1+(tick_size/2)])
    
    n_cells_stepsize = 20
    ticks = np.arange(n_cells_stepsize - 1, n_cells_region, n_cells_stepsize)
    labels = ticks + 1
    
    ax.set_yticks(ticks = ticks, labels = labels)
    ax.set_yticks(ticks = np.arange(1, n_cells_region, 2),  minor = True)
    ax.set_ylabel(f'{region}\ncells [#]', fontsize = 12)
    
    # edit spines of eventplots
    ax.spines['left'].set_bounds([0, n_cells_region-1])
    ax.spines['bottom'].set_bounds([0, 30])
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# set eventplots x axes

# remove ticks between first two activity plots
for i in ['A', 'B']:
    axs_regions[i].set_xticks(ticks = [], labels = [])
    axs_regions[i].spines['bottom'].set_visible(False)
    axs_regions[i].tick_params(axis = 'x', size = 0)
    axs_regions[i].tick_params(axis = 'x', which = 'minor', size = 0)

axs_regions['C'].set_xlabel('Time [s]')



# # # v_rest # # #

# set axis
ax = axs_regions['D']

swarms = sbn.swarmplot(data = activity_df, 
                       x = "Region", 
                       y = "v_rest",
                       hue = 'Region',
                       palette = region_colors,
                       size = 3, 
                       ax = ax, 
                       color = colors_dict['primecolor'],
                       dodge = True,
                       zorder = 1)

# set offsets of violins
# v_positions = [-0.26667, 0, 0.26667]

v_positions = [-0.26667, 1, 2.26667]


# for a_idx, activity in enumerate(['silent', 'spiking']):
    
for r_idx, region in enumerate(regions):
    
    # get data
    violin_data = activity_df[activity_df['Region'] == region]
    violin_data = violin_data['v_rest'].to_list()
    
    # get mean, median, and std
    mean = np.mean(violin_data)
    median = np.median(violin_data)
    std = np.std(violin_data)
    n = len(violin_data)

    # set offset
    offset = 0.25

    # set bandwidth
    bandwidth = n**(-1./(1+2))
    # print(n**(-1./(1+4)), bandwidth)

    # plot half violin
    plot_half_violin(data = violin_data, 
                     ax = ax,
                     v_kde_cutoff = 0.01,
                     v_direction = 1,
                     v_bandwidth = bandwidth,
                     v_position = v_positions[r_idx],# + a_idx,
                     v_offset = offset,
                     v_width = 0.3,
                     v_baseline = False,
                     v_color = region_colors[region],
                     v_zorder = 0,
                     v_lw = 1.)
    
    # set x position of errorbar
    e_position = v_positions[r_idx] + offset #+ a_idx
    
    ax.errorbar(x = e_position,
                y = mean,
                yerr = std,
                fmt='_', 
                markersize = 6,
                markerfacecolor = 'none',
                capsize = 2,
                color=colors_dict['primecolor'],
                linewidth = 1.0,
                label = '_nolegend_',
                zorder = 2)
    
    ax.scatter(x = e_position,
                y = median,
                marker='D', 
                s = 5,
                color=colors_dict['primecolor'],
                linewidth = 1,
                label = '_nolegend_',
                zorder = 3)
                
        
# remove seaborn legend
# ax.get_legend().remove()
ax.legend(prop={'size': 10})
ax.get_legend().get_frame().set_linewidth(0.0)


# y
ydict = {'ax_min' : -100,
         'ax_max' : -50,
         'pad' : None,
         'step' : 10,
         'stepminor' : 2,
         'label' : 'Resting membrane potential\n[mV]'}

apply_axis_settings(ax, axis = 'y', **ydict)

# x
# xdict = {'ax_min' : 0,
#          'ax_max' : 2,
#          'pad' : 1.5,
#          'step' : 1,
#          'stepminor' : 1,
#          'label' : 'Region'}

# apply_axis_settings(ax, axis = 'x', **xdict)

ax.set_xlim([0 - 0.85, 2 + 2])
ax.set_xticks(ticks = np.arange(0, 2+ 1, 1))
ax.set_xticks(ticks = np.arange(0, 2+ 1, 1), minor = True)
ax.spines['bottom'].set_bounds([0, 2])

# despine
[axs_regions['D'].spines[spine].set_visible(False) for spine in ['top', 'right']]

# set grid as False for all subplots
[axs_regions[ax].grid(False) for ax in axs_regions]

# align labels
fig_regions.align_labels()

plt.show()

path_figure = join(figure_dir, 'temp_figs')

save_figures(fig_regions, 'Resting_n_eventplot_nspikes+Regions' + exp[:-1], path_figure, darkmode_bool,
              figure_format= 'both',
              dataframe_to_save = activity_df, index_label = 'cell_ID', add_measures = True, axis_for_calcs = 0,
              groups_bool= True, groups= ['BAOT/MeA', 'MeA', 'BAOT'], groups_name= 'Region')



