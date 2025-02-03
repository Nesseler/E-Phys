# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:43:18 2025

@author: nesseler
"""

# init plotting
from functions.initialize_plotting import *
from os.path import join
import numpy as np


# %% define axis settings

# y
ydict_v = {'ax_min' : -150,
           'ax_max' : -50,
           'pad' : None,
           'step' : 20,
           'stepminor' : 5,
           'label' : 'Membrane\npotential [mV]'}

# x
xdict_tfull = {'ax_min' : 0,
               'ax_max' : 1500,
               'pad' : None,
               'step' : 250,
               'stepminor' : 50,
               'label' : 'Time [ms]'}


# %% full protocoll plot

def plot_full_sag(cell_ID, t, v, i, i_input, gradient = True):
    '''
    This function creates a figure displaying the entire protocol, optionally 
    with color-coded steps if gradient bool is set.
    Parameters:
        cell_ID: str, like 'E-303', unifque cell identifier
        t: numpy array, time dimension
        v: numpy nd array, voltage traces for all steps
        i: numpy nd array, current traces for all steps
        i_input : numpy array, current input for each step
        gradient: bool, activates color-coded steps
    '''

    # t = t_full
    # v = v_full
    # i = i
    # i_input = i_input
    # gradient = True
    
    # get number of steps
    n_steps = v.shape[0]
    
    # init figure and axes
    fig, axs = plt.subplots(nrows = 2,
                            ncols = 1,
                            layout = 'constrained',
                            figsize = get_figure_size(width = 160, height = 100),
                            sharex = True,
                            height_ratios = [1, 5],
                            dpi = 300)
    
   
    # set axis title
    axs[0].set_title(f'{cell_ID} cc_IF',
                     fontsize=9, 
                     loc='left',
                     x = 0.02)
    
    # color code
    if gradient:
        # specify linecollection settings
        lc_dict = {'lw' : 0.5,
                   'linestyle' : 'solid',
                   'array' : i_input,
                   'cmap' : 'rainbow_r'}
    else:
        # specify linecollection settings
        lc_dict = {'lw' : 0.5,
                   'linestyle' : 'solid',
                   'color' : colors_dict['primecolor']}
        
        
    # define line collections
    v_collection = LineCollection([np.column_stack([t, v[step]]) for step in range(n_steps)],
                                  **lc_dict)
    
    i_collection = LineCollection([np.column_stack([t, i[step]]) for step in range(n_steps)], 
                                  **lc_dict) 
    
    # voltage
    ax = axs[1]
    
    # add line collection
    ax.add_collection(v_collection)
    
    apply_axis_settings(ax, axis = 'y', **ydict_v)
    
    # current
    ax = axs[0]
    
    # add line collection
    ax.add_collection(i_collection)
    
    # y
    ydict_i = {'ax_min' : -1000,
               'ax_max' : 100,
               'pad' : 10,
               'step' : 1000,
               'stepminor' : 100,
               'label' : 'Input\ncurrent [pA]'}
    
    apply_axis_settings(ax, axis = 'y', **ydict_i)
    
    # x 
    apply_axis_settings(axs[-1], axis = 'x', **xdict_tfull)
    remove_spines_n_ticks([axs[0]], axis = 'x')
    
    # remove spines
    [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]
    
    # align labels
    fig.align_labels()
    
    # create saving path and save
    from parameters.directories_win import vplot_dir
    path_fig = join(vplot_dir, 'cc_sag', 'traces')
    save_figures(fig, f'{cell_ID}-cc_sag', path_fig, darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()
    

# %% passive properties verification plot


def plot_passiv_props_calc(cell_ID, t_full, v_full, t_stim, v_stim, passive_props_calcs, SR):
    '''
    This function creates a figure displaying the hyperpolarising steps for 
    the calculation of passive membrane properties.
    Parameters:
        cell_ID: str, like 'E-303', unifque cell identifier
        t_full: numpy array, time dimension
        v_full: numpy nd array, voltage traces for all steps
        t_stim: numpy array, time dimension for the stimulation period
        v_stim: numpy array, numpy nd array, voltage traces for all steps and
                stimulation period
        passive_props_calcs: pandas dataframe, contains measurements for each step
        SR: float, sampling rate
    '''

    from functions.functions_useful import exp_func

    # init figure and axes
    fig, axs = plt.subplots(nrows = 2,
                            ncols = 3,
                            layout = 'constrained',
                            figsize = get_figure_size(width = 160, height = 100),
                            sharex = True,
                            sharey = True,
                            dpi = 300)
    
    # set figure title
    fig.suptitle(f'{cell_ID} passive properties measurements',
                 fontsize = 9)
    
    # iterate through useable steps
    for step in passive_props_calcs.index:
        
        # get popt for step
        popt = passive_props_calcs.at[step, 'popt']
        
        # index of voltag minimum
        idx_vmin = passive_props_calcs.at[step, 'idx_vmin']
        
        # flatten axes array
        axs = axs.flatten()
        
        # set axis
        ax = axs[step-1]
        
        # set axis title
        ax.set_title(f'step {step}',
                     fontsize=6, 
                     loc='left',
                     x = 0.03,
                     y = 0.88)
        
        # plot
        ax.plot(t_full, v_full[step],
                lw = 0.5,
                color = colors_dict['primecolor'])
        
        # plot trace to min
        ax.plot(t_stim[:idx_vmin], v_stim[step][:idx_vmin],
                lw = 0.75,
                color = colors_dict['color3'],
                linestyle = 'solid')
               
        # plot exponential fit
        ax.plot(t_stim, exp_func(np.arange(0, 1000 * (SR/1e3)), *popt),
                lw = 0.5,
                color = colors_dict['color2'],
                linestyle = 'dashed')
    
        # add inset
        # inset marker
        box_xmin   = 240
        box_width  = 250 
        box_ymin   = -115
        box_height = 35
           
        ## ([left, bottom, width, height]), percentages
        ax_inset = ax.inset_axes([0.10, 0.05, 0.6, 0.40],
                                 xlim=(box_xmin, box_xmin+box_width), 
                                 ylim=(box_ymin, box_ymin+box_height), 
                                 xticklabels=[], 
                                 yticklabels=[])
        
        # edit linewidth of inset axis and its ticks
        [ax_inset.spines[spine].set_linewidth(0.5) for spine in ['left', 'bottom']]
        ax_inset.tick_params(width=0.5)
        ax_inset.tick_params(which = 'minor', width=0.25)
    
        
        # add rectangle marker
        ax.add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                               width = box_width, 
                               height = box_height,
                               fill = False,
                               color = colors_dict['primecolor'],
                               linestyle = '--',
                               lw = 0.5,
                               alpha = 0.5))
        
        # set axis
        ax = ax_inset
        
        # plot
        ax.plot(t_full, v_full[step],
                lw = 0.5,
                color = colors_dict['primecolor'])
        
        # plot trace to min
        ax.plot(t_stim[:idx_vmin], v_stim[step][:idx_vmin],
                lw = 0.75,
                color = colors_dict['color3'],
                linestyle = 'solid')
            
        # plot exponential fit
        ax.plot(t_stim, exp_func(np.arange(0, 1000 * (SR/1e3)), *popt),
                lw = 0.5,
                color = colors_dict['color2'],
                linestyle = 'dashed')
        
        # x
        ax_inset.set_xticks(ticks = np.arange(box_xmin, box_xmin + box_width, 500), labels = [])
        ax_inset.set_xticks(ticks = np.arange(box_xmin, box_xmin + box_width, 50), labels = [], minor = True)
        ax_inset.set_xlim([box_xmin, box_xmin + box_width])
        
        # y
        ax_inset.set_yticks(ticks = np.arange(box_ymin, box_ymin + box_height, 20), labels = [])
        ax_inset.set_yticks(ticks = np.arange(box_ymin, box_ymin + box_height, 5), labels = [], minor = True)
        ax_inset.set_ylim([box_ymin, box_ymin + box_height])
         
        # remove inset spines
        [ax_inset.spines[spine].set_visible(False) for spine in ['top', 'right']]
        
        # TODO: add text of measurements
        
    
    # x
    xdict = {'ax_min' : 0,
             'ax_max' : 1500,
             'pad' : None,
             'step' : 500,
             'stepminor' : 100,
             'label' : ''}
    
    [apply_axis_settings(ax, axis = 'x', **xdict) for ax in axs]
    
    # y
    ydict = {'ax_min' : -150,
             'ax_max' : -75,
             'pad' : None,
             'step' : 20,
             'stepminor' : 5,
             'label' : ''}
    
    [apply_axis_settings(ax, axis = 'y', **ydict) for ax in axs]
    
    # remove spines
    [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]
    
    # align labels
    fig.align_labels()
    
    # TODO: saving
    
    # display figure
    plt.show()