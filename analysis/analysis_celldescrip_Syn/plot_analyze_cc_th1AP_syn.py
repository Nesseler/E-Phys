# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 18:21:39 2025

@author: nesseler
"""

# init plotting
from functions.initialize_plotting import *
from os.path import join
import numpy as np


# %% define axis settings

# y
ydict_v = {'ax_min' : -100,
           'ax_max' : 60,
           'pad' : None,
           'step' : 50,
           'stepminor' : 5,
           'label' : 'Membrane\npotential [mV]'}

# x
xdict_tfull = {'ax_min' : 0,
               'ax_max' : 500,
               'pad' : None,
               'step' : 100,
               'stepminor' : 50,
               'label' : 'Time [ms]'}

# %%


def plot_full_th1AP(cell_ID, t, v, i, i_input, gradient = True):
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

    # v = v_full
    # t = t_full
    # gradient = True
    # i = i_calc

    # get number of steps
    n_steps = v.shape[0]
    
    # init figure and axes
    fig, axs = plt.subplots(nrows = 2,
                            ncols = 1,
                            layout = 'constrained',
                            figsize = get_figure_size(width = 160, height = 100),
                            sharex = True,
                            height_ratios = [1, 5])
    
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
                   'cmap' : 'rainbow'}
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
    ydict_i = {'ax_min' : -100,
               'ax_max' : 1000,
               'pad' : 10,
               'step' : 500,
               'stepminor' : 100,
               'label' : 'Input\ncurrent [pA]',
               'start_at_0' : True}
    
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
    path_fig = join(vplot_dir, 'cc_th1AP-traces')
    save_figures(fig, f'{cell_ID}-cc_th1AP', path_fig, darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()
    
    
# %%

def plot_rheospike(cell_ID, t, v, dvdt, spike_t, spike_v, spike_dvdt):
    '''
    This function creates a figure displaying the entire protocol, optionally 
    with color-coded steps if gradient bool is set.
    Parameters:
        cell_ID: str, like 'E-303', unifque cell identifier
        t: numpy array, time dimension
        v: numpy array, for voltage trace of rheobase step
        dvdt : numpy array, for first derivate of voltage trace of rheobase step
        t_spike : numpy array, time dimension of spike
        v_spike : numpy array, for voltage trace of spike
        dvdt_spike : numpy array, for first derivate of voltage trace of spike
    '''

    # set plotting dict
    line_dict = {'lw' : 0.5,
                 'linestyle' : 'solid'}
            
    
    # init figure and axes
    fig, axs = plt.subplots(nrows = 1,
                            ncols = 2,
                            layout = 'constrained',
                            figsize = get_figure_size(width = 180, height = 100),
                            width_ratios = [1, 1.5],
                            dpi = 300)
    
    # set figure title
    fig.suptitle(f'{cell_ID} cc_th1AP')
    
    # voltage
    
    # set axis title
    axs[0].set_title(f'cc_th1AP',
                     fontsize=9, 
                     loc='left',
                     x = 0.02)
    
    axs[0].plot(t, v,
                color = colors_dict['primecolor'],
                **line_dict)
    
    axs[0].plot(spike_t, spike_v,
                color = '#FFEC9DFF',
                **line_dict)
    
    # edit axis
    apply_axis_settings(axs[0], axis = 'y', **ydict_v)
    
    # x
    xdict_t = {'ax_min' : 240,
               'ax_max' : 300,
               'pad' : None,
               'step' : 20,
               'stepminor' : 5,
               'label' : 'Time [ms]'}
    
    apply_axis_settings(axs[0], axis = 'x', **xdict_t)
    
    # phase plane
    
    # set axis title
    axs[1].set_title(f'phase plane',
                     fontsize=9, 
                     loc='left',
                     x = 0.02)
    
    axs[1].plot(v, dvdt,
                color = colors_dict['primecolor'],
                **line_dict)
    
    axs[1].plot(spike_v, spike_dvdt,
                color = '#FFEC9DFF',
                **line_dict)
    
    # y
    dvdt_dict = {'ax_min' : -150,
                 'ax_max' : 250,
                 'pad' : None,
                 'step' : 50,
                 'stepminor' : 10,
                 'label' : 'Rate of membrane\npotential change [mV/ms]'}
    
    apply_axis_settings(axs[1], axis = 'y', **dvdt_dict)
    
    # edit axis
    apply_axis_settings(axs[1], axis = 'x', **ydict_v)
    
    # inset
    
    # inset marker
    box_xmin   = -75
    box_width  = 25 
    box_ymin   = -10
    box_height = 20
    
    # add rectangle marker
    axs[1].add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                               width = box_width, 
                               height = box_height,
                               fill = False,
                               color = colors_dict['primecolor'],
                               linestyle = '--',
                               lw = 0.5))
    
    ## ([left, bottom, width, height]), percentages
    ax_inset = fig.add_axes([0.58, 0.73, 0.13, 0.145])
    
    
    # plot full trace
    ax_inset.plot(v, dvdt,
                  color = colors_dict['primecolor'],
                  lw = 0.5)
    
    # plot rheobase spike
    ax_inset.plot(spike_v, spike_dvdt,
                  color = '#FFEC9DFF',
                  lw = 1)
    
    # x
    ax_inset.set_xticks(ticks = np.arange(box_xmin, box_xmin + box_width, 20), labels = [])
    ax_inset.set_xticks(ticks = np.arange(box_xmin, box_xmin + box_width, 5), labels = [], minor = True)
    ax_inset.set_xlim([box_xmin, box_xmin + box_width])
    
    # y
    ax_inset.set_yticks(ticks = np.arange(box_ymin, box_ymin + box_height, 10), labels = [])
    ax_inset.set_yticks(ticks = np.arange(box_ymin, box_ymin + box_height, 2), labels = [], minor = True)
    ax_inset.set_ylim([box_ymin, box_ymin + box_height])
     
    # remove inset spines
    [ax_inset.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # remove spines
    [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]
    
    # align labels
    fig.align_labels()
    
    # create saving path and save
    from parameters.directories_win import vplot_dir
    path_fig = join(vplot_dir, 'cc_th1AP-rheobase_spike')
    save_figures(fig, f'{cell_ID}-cc_th1AP-rheobase_spike', path_fig, darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()