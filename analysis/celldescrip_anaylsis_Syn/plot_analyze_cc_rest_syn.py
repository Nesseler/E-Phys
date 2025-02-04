# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:20:02 2025

@author: nesseler
"""

# init plotting
from functions.initialize_plotting import *
from os.path import join
import numpy as np

def create_cc_rest_vplot(cell_ID, t, vf, vf_wo_spikes, t_spikes, v_rest):
    '''
    This functions create the verification plot for the cc_rest protocol analysis,
    which calulates the v_rest and excludes spikes that occurred.
    '''

    fig, ax = plt.subplots(nrows = 1,
                           ncols = 1,
                           figsize = get_figure_size(width = 160, height = 100),
                           layout = 'constrained')
    
    # set axis title
    ax.set_title(f'{cell_ID} cc_rest',
                 fontsize=9, 
                 loc='left',
                 x = 0.02)
    
    linedict = {'lw' : 0.5}
    
    ax.plot(t, vf, 
            color = 'grey',
            **linedict, label = 'recorded')
    
    if np.isnan(vf_wo_spikes).all():
        label_vf_wo_spikes = '_nolegend_'
        label_t_spikes = '_nolegend_'
    else:
        label_vf_wo_spikes = 'spikes excluded'
        label_t_spikes = 'spike times'
        
    ax.plot(t, vf_wo_spikes,
            color = colors_dict['primecolor'],
            **linedict, label = label_vf_wo_spikes)

    # eventplot
    ax.eventplot(t_spikes,
                 orientation = 'horizontal',
                 lineoffsets = 60-3,
                 linelengths = 3,
                 color = 'r',
                 lw = 1,
                 label = label_t_spikes)
    
    ax.hlines(xmin = 0,
              xmax = 30,
              y = v_rest,
              lw = 0.5,
              linestyle = 'dashed',
              color = 'r',
              zorder = 3,
              label = '$v_{rest}$')
    
    # y
    ydict = {'ax_min' : -100,
             'ax_max' : 60,
             'pad' : 10,
             'step' : 20,
             'stepminor' : 5,
             'label' : 'Membrane potential [mV]'}
    
    apply_axis_settings(ax, axis = 'y', **ydict)
    
    # x
    xdict = {'ax_min' : 0,
             'ax_max' : 30,
             'pad' : 0.5,
             'step' : 10,
             'stepminor' : 5,
             'label' : 'Time [s]'}
    
    apply_axis_settings(ax, axis = 'x', **xdict)
    
    # legend
    fig.legend(fontsize = 8, frameon = False)
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # create saving path and save
    from parameters.directories_win import vplot_dir
    path_fig = join(vplot_dir, 'cc_rest-v_rest')
    save_figures(fig, f'{cell_ID}-v_rest', path_fig, darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()