# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 09:23:38 2025

@author: nesseler
"""

import pandas as pd
import numpy as np
from os.path import join
import scipy as sc
import warnings
from tqdm import tqdm

from parameters.directories_win import vplot_dir

# custom functions
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_useful import butter_filter

# define protocol
PGF = 'cc_rest_adaEk'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF + '_pre', sheet_name = 'PGFs_Syn')

# get number of cells
n_cells = len(cell_IDs)

# init plotting
from functions.initialize_plotting import * # analysis:ignore

# PGF = PGF + '_pre'


cc_rest_adaEk_parameters = {'t'       : 30,
                            'i_hold'  : 0 ,
                            'SR'      : 50000,
                            'n_pre'   : 2,
                            'n_washin': 6,
                            'n_post'  : 2}

conditions = ['pre', 'washin', 'post']


# prepare dataframes for loading data
t_s = {'pre'    : np.arange(0, cc_rest_adaEk_parameters['t']*cc_rest_adaEk_parameters['n_pre'], 1/cc_rest_adaEk_parameters['SR']),
       'washin' : np.arange(0, cc_rest_adaEk_parameters['t']*cc_rest_adaEk_parameters['n_washin'], 1/cc_rest_adaEk_parameters['SR']),
       'post'   : np.arange(0, cc_rest_adaEk_parameters['t']*cc_rest_adaEk_parameters['n_post'], 1/cc_rest_adaEk_parameters['SR'])}

v_dfs = {'pre'    : pd.DataFrame(columns=cell_IDs, index = t_s['pre']),
         'washin' : pd.DataFrame(columns=cell_IDs, index = t_s['washin']),
         'post'   : pd.DataFrame(columns=cell_IDs, index = t_s['post'])}




# %%

print('loading ... ')

for cell_ID in tqdm(cell_IDs):

    for condition in conditions:
    
        # get the traceIndex and the file path string for data import functions
        traceIndex, file_path = get_traceIndex_n_file(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')
           
        # check for multiple protocols in list
        protocols = type(traceIndex[1])
        
        if protocols == int: #and (condition == 'pre' or condition == 'post'):
            
            # load and check for number of steps
            # get data with file path & trace index
            i, v, t, SR, n_step = get_cc_data(file_path, traceIndex, scale='s')
            
            # concatenate steps
            v = v.flatten('C')
            
            if n_step == 1:
                # pads with nans
                v = np.pad(v,
                            pad_width = (0, 1500000),
                            mode = 'constant',
                            constant_values = (np.nan,np.nan))
                
            # filter all data with 1kHz cutoff
            vf = butter_filter(v, order=3, cutoff=1e3, sampling_rate=SR)
        
            # replace first values with nans to eliminate filter artifact
            vf[:100] = np.nan
        
            # write to dataframe
            v_dfs[condition][cell_ID] = vf 
              
            
        elif protocols == list:
            
            # get number of protocols
            n_protocols = len(traceIndex[1])
            
            # check if number of protocols and expected match
            # if so protocols can simpy be concatenated
            if n_protocols == cc_rest_adaEk_parameters['n_' + condition]:
            
                # create numpy array to concatenate to
                v_empty = np.empty(0)
                
                # iterate through steps
                for istep, step in enumerate(traceIndex[1]):
                    
                    # redefine traceIndex
                    local_TI = traceIndex
                    local_TI[1] = step
                    
                    # get data with file path & trace index
                    i, v, t, SR, n_step = get_cc_data(file_path, local_TI, scale='s')
                    
                    # simply append consecutively recorded protocols
                    v_empty = np.append(v_empty, v)
                
                # filter all data with 1kHz cutoff
                vf = butter_filter(v_empty, order=3, cutoff=1e3, sampling_rate=SR)
        
                # replace first values with nans to eliminate filter artifact
                vf[:100] = np.nan
                
                # write to dataframe
                v_dfs[condition][cell_ID] = vf
                
            # if not 30 sec intervals were between protocols
            elif n_protocols == cc_rest_adaEk_parameters['n_' + condition] / 2:
                
                # create numpy array filled with nans
                v_empty = np.zeros_like(t_s[condition])
                v_empty.fill(np.nan)
                
                # iterate through steps
                for istep, step in enumerate(traceIndex[1]):
                    
                    # redefine traceIndex
                    local_TI = traceIndex
                    local_TI[1] = step
                    
                    # get data with file path & trace index
                    i, v, t, SR, n_step = get_cc_data(file_path, local_TI, scale='s')
                    
                    # filter all data with 1kHz cutoff
                    vf = butter_filter(v, order=3, cutoff=1e3, sampling_rate=SR)
        
                    # create index array for insertion
                    start = 0  + (30*istep*2)
                    stop =  30 + (30*istep*2)
                    v_idc = np.arange(start * SR, stop * SR)
                    
                    # write to v_empty
                    v_empty[v_idc] = vf
                    
                    # replace first values with nans to eliminate filter artifact
                    v_empty[start*SR:start*SR+100] = np.nan
            
                # write to dataframe
                v_dfs[condition][cell_ID] = v_empty 


# %% get means

def get_means_by_interval(data, sampling_interval = 30, SR = 50000):

    # get number of possible means
    n_means = (len(data) / SR) / sampling_interval
    
    # initialise list
    means = [None] * int(n_means)
    times = [None] * int(n_means)
    
    i = 0
    while i < n_means:
        
        start = 0  + (sampling_interval * SR * i)
        stop  = sampling_interval * SR  + (sampling_interval * SR * i)
        
        means[i] = np.nanmean(data[start:stop])
        times[i] = (sampling_interval * i) + (sampling_interval / 2)
        
        i+=1
        
    return means, times

int_1 = 30
int_2 = 10

# 30s means
mean_vmem_30s = {'pre'    : pd.DataFrame(columns=cell_IDs, index = np.arange(int_1/2, 60, int_1)),
                 'washin' : pd.DataFrame(columns=cell_IDs, index = np.arange(int_1/2, 180, int_1)),
                 'post'   : pd.DataFrame(columns=cell_IDs, index = np.arange(int_1/2, 60, int_1))}

# 5s means
mean_vmem_int2 = {'pre'    : pd.DataFrame(columns=cell_IDs, index = np.arange(int_2/2, 60, int_2)),
                  'washin' : pd.DataFrame(columns=cell_IDs, index = np.arange(int_2/2, 180, int_2)),
                  'post'   : pd.DataFrame(columns=cell_IDs, index = np.arange(int_2/2, 60, int_2))}

for cell_ID in cell_IDs:
    for condition in conditions:
        
        # define data
        data = v_dfs[condition][cell_ID].to_numpy()
        
        # get means and timings
        means30, times30 = get_means_by_interval(data, int_1, 50000)
        means_int2, times_int2 = get_means_by_interval(data, int_2, 50000)
        
        # write to dataframe
        mean_vmem_30s[condition].loc[times30, cell_ID] = means30
        mean_vmem_int2[condition].loc[times_int2, cell_ID] = means_int2
    
         
# %% plot

from functions.initialize_plotting import *

print('creating figures ...')

for cell_ID in tqdm(cell_IDs):

    fig, axs = plt.subplots(nrows = 4,
                            ncols = 3,
                            figsize = get_figure_size(),
                            dpi = 300,
                            layout = 'constrained',
                            width_ratios= [2, 6, 2],
                            sharex = 'col',
                            sharey = 'row')
    
    fig.suptitle(cell_ID)
    
    # flatten array of axes
    axs = axs.flatten()
    
    axs[0].set_title('$E_k$ = -99 mV')
    axs[1].set_title('washin')
    axs[2].set_title('$E_k$ = -85 mV')
    
    
    # plot original traces
    for c_idx, condition in enumerate(conditions):
        for ax in axs[c_idx:6:3]:
            ax.plot(t_s[condition], v_dfs[condition][cell_ID], 
                    lw = 0.5,
                    c = colors_dict['primecolor'])
    
    lineplot_dict = {'ms' : 7,
                     'marker' : 'o',
                     'mew' : 1,
                     'mfc' : 'k',
                     'mec' : colors_dict['primecolor'], 
                     'ls' : 'dashed',
                     'dashes' : (3,4),
                     'lw' : 1.5,
                     'c' : colors_dict['primecolor']
                     }
    
    # plot means
    for c_idx, condition in enumerate(conditions):
        
        # define axis
        ax = axs[c_idx+6]
        
        # add means
        ax.plot(mean_vmem_int2[condition][cell_ID].index.to_numpy(), 
                mean_vmem_int2[condition][cell_ID].to_numpy(), 
                **lineplot_dict)
        
    # plot means
    for c_idx, condition in enumerate(conditions):
        
        # define axis
        ax = axs[c_idx+9]
        
        # add means
        ax.plot(mean_vmem_30s[condition][cell_ID].index.to_numpy(), 
                mean_vmem_30s[condition][cell_ID].to_numpy(),
                **lineplot_dict)
    
        
        
    # y
    ydict_full = {'ax_min' : -100,
                  'ax_max' : 50,
                  'pad' : 2,
                  'step' : 50,
                  'stepminor' : 5,
                  'label' : '$V_{mem}$ [mV]'}
    
    # edit axis
    apply_axis_settings(axs[0], axis = 'y', **ydict_full)
    
    
    # y
    ydict = {'ax_min' : -90,
             'ax_max' : -70,
             'pad' : 0.5,
             'step' : 10,
             'stepminor' : 1,
             'label' : '$V_{mem}$ [mV]'}
    
    # edit axis
    apply_axis_settings(axs[3], axis = 'y', **ydict)
    apply_axis_settings(axs[6], axis = 'y', **ydict)
    apply_axis_settings(axs[9], axis = 'y', **ydict)
    
    
    # x (pre & post)
    xdict_pp = {'ax_min' : 0,
                'ax_max' : 60,
                'pad' : 60 * 0.02,
                'step' : 30,
                'stepminor' : 5,
                'label' : ''}
    
    # edit axis
    for ax in axs[::3]:
        apply_axis_settings(ax, axis = 'x', **xdict_pp)
    
    for ax in axs[2::3]:
        apply_axis_settings(ax, axis = 'x', **xdict_pp)
    
    
    
    # x (washin)
    xdict_washin = {'ax_min' : 0,
                    'ax_max' : 180,
                    'pad' : 180*0.02,
                    'step' : 30,
                    'stepminor' : 5,
                    'label' : ''}
    
    # edit axis
    for ax in axs[1::3]:
        apply_axis_settings(ax, axis = 'x', **xdict_washin)
    
    # remove spines
    [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]
    
    for ax in axs[1::3]:
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis = 'y', size = 0)
        ax.tick_params(axis = 'y', which = 'minor', size = 0)
        
    for ax in axs[2::3]:
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis = 'y', size = 0)
        ax.tick_params(axis = 'y', which = 'minor', size = 0)
        
    
    # set x label
    for ax in axs[9:11]:
        ax.set_xlabel('Time [s]')
        
    # align labels
    fig.align_labels()
    
    # create saving path and save
    vplots_path_fig = join(vplot_dir, 'adaEk')
    save_figures(fig, f'{cell_ID}-adaEK_washin', vplots_path_fig, darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()


# %% dot plot 



pre_vmem = mean_vmem_30s['pre'].loc[15]

post_vmem = mean_vmem_30s['post'].loc[15]



fig, ax = plt.subplots(nrows = 1,
                        ncols = 1,
                        layout = 'constrained', 
                        figsize = get_figure_size(width = 165.5))

#
swarm_dict = {'size' : 7,
              'marker' : 'o', 
              'edgecolor' : colors_dict['primecolor'],
              'facecolor' : 'k',
              'linewidth' : 1}

# plot pre data points as swarm
sbn.swarmplot(y = pre_vmem,
              x = [0] * len(pre_vmem),
              ax = ax,
              **swarm_dict)

# get positions of pre data points
pre_positions = np.array(ax.collections[0].get_offsets())

# plot post data points as swarm
sbn.swarmplot(y = post_vmem,
              x = [1] * len(post_vmem),
              ax = ax,
              **swarm_dict)

# get positions of post data points
post_positions = np.array(ax.collections[1].get_offsets())


# remap positions and cell_IDs
pre_positions_df = pd.DataFrame(columns = ['x', 'y'], index = cell_IDs)
post_positions_df = pd.DataFrame(columns = ['x', 'y'], index = cell_IDs)


for cell_ID in cell_IDs:
    
    # get pre
    pre_v = pre_vmem[cell_ID]
    post_v = post_vmem[cell_ID]
    
    # find in positions
    x_pre_idx = np.where(pre_positions == pre_v)[0][0]
    x_post_idx = np.where(post_positions == post_v)[0][0]
    
    # get x position
    x_pre = pre_positions[x_pre_idx][0]
    x_post = post_positions[x_post_idx][0]
    
    # write to dataframe
    pre_positions_df.at[cell_ID, 'x'] = x_pre
    pre_positions_df.at[cell_ID, 'y'] = pre_v
    post_positions_df.at[cell_ID, 'x'] = x_post
    post_positions_df.at[cell_ID, 'y'] = post_v
    

    
pad_perc = 0.1


def get_lines_for_paired_data(point1, point2, pad_perc):
    """
    This function calculates the start and end points of line that connect
    paired data points.
    Parameters:
        point1, point2: list, x and y coordinates of data points
        pad_perc: float, padding percentage on x axis between data points and
                  connecting lines
    Return:
        two lists of x and y coordinates for plotting
    """
    
    # get x and y coordinates from both points
    (xs, ys) = zip(point1, point2)
    
    # calc distance in x between points to enforce padding
    x_distance = xs[1] - xs[0]
    
    # calc linestart x coordinate
    linestart_x = xs[0] + x_distance * pad_perc
    
    # use linestart in x to interpolate the corresponding y
    linestart_y = np.interp(linestart_x, xs, ys)
   
    # calc linestop x coordinate
    linestop_x = xs[1] - x_distance * pad_perc
    
    # use linestart in x to interpolate the corresponding y
    linestop_y = np.interp(linestop_x, xs, ys)
    
    return [linestart_x, linestop_x], [linestart_y, linestop_y]
    

for cell_ID in cell_IDs:

    # get points 
    point1 = pre_positions_df.loc[cell_ID, ['x', 'y']]
    point2 =  post_positions_df.loc[cell_ID, ['x', 'y']]
    
    # get coordinates for lines
    line_xs, line_ys = get_lines_for_paired_data(point1, point2, pad_perc)
    
    # plot connecting lines
    ax.plot(line_xs, line_ys,
            lw = 1,
            c = 'grey')
    
    
for vmem, direction, pos in zip([pre_vmem, post_vmem], [-1, 1], [0, 1]):
    
    # plot half violin pre
    plot_half_violin(data = vmem.to_list(), 
                     ax = ax,
                     v_direction = direction,
                     v_resolution = 0.1,
                     v_kde_cutoff = 0.1,
                     v_abs_cutoff = [-90, np.nan],
                     v_position = pos,
                     v_offset = 0.15 * direction,
                     v_width = 0.2,
                     v_baseline = False,
                     v_color = colors_dict['primecolor'],
                     v_zorder = 0)  

    e_position = pos + (0.125 * direction)

    # plot errorbar
    ax.errorbar(x = e_position,
                y = vmem.mean(),
                yerr = vmem.std(),
                fmt='_', 
                markersize = 7,
                markerfacecolor = 'none',
                capsize = 3,
                color=colors_dict['primecolor'],
                linewidth = 1,
                label = '_nolegend_',
                zorder = 2)
    
    # plot median
    ax.scatter(x = e_position,
               y = vmem.median(),
               marker='D', 
               s = 11,
               color=colors_dict['primecolor'],
               linewidth = 1,
               label = '_nolegend_',
               zorder = 3)

# y
ydict = {'ax_min' : -90,
         'ax_max' : -70,
         'pad' : 0.5,
         'step' : 10,
         'stepminor' : 1,
         'label' : '$V_{mem}$ [mV]'}

apply_axis_settings(ax, axis = 'y', **ydict)


# x
xdict = {'ax_min' : 0,
         'ax_max' : 1,
         'pad' : 0.5,
         'step' : 1,
         'stepminor' : 1,
         'label' : ''}

apply_axis_settings(ax, axis = 'x', **xdict)

# set x ticks
ax.set_xticklabels(['$E_k$ = -99 mV', '$E_k$ = -85 mV'])

# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

plt.show()