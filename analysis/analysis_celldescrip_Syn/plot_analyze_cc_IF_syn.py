# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 13:13:24 2025

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
           'label' : 'Membrane potential [mV]'}

# x
xdict_tfull = {'ax_min' : 0,
               'ax_max' : 1500,
               'pad' : 10,
               'step' : 500,
               'stepminor' : 50,
               'label' : 'Time [ms]'}


# %% full protocoll plot

def plot_full_IF(cell_ID, t, v, i, i_input, gradient = True):
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
    path_fig = join(vplot_dir, 'cc_IF-traces')
    save_figures(fig, f'{cell_ID}-cc_IF', path_fig, darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()
    
    
# %% step plot

def plot_IF_step_spike_detection(cell_ID, step, t, v, t_spikes, 
                                 freq, inst_freq, init_inst_freq):
    '''
    This function creates a figure displaying a single step of the IF protocol
    and marks the detected spikes in the given trace.
    Parameters:
        cell_ID: str, like 'E-303', unifque cell identifier
        step: int, step index of the given trace 
        t: numpy array, time dimension
        v: numpy array, voltage trace for single step
        t_spikes: numpy array, list of spike times (relative to step onset)
        freq, inst_freq, initial_inst_freq: float, calculated spiking frequencies
    '''
    
    # init figure and axes
    fig, ax = plt.subplots(nrows = 1,
                           ncols = 1,
                           layout = 'constrained',
                           figsize = get_figure_size(width = 160, height = 100))
    
    # set axis title
    ax.set_title(f'{cell_ID} cc_IF step no.: {step}',
                 fontsize=9, 
                 loc='left',
                 x = 0.02)
    
    # plot
    ax.plot(t, v,
            color = colors_dict['primecolor'],
            lw = 0.5)
    
    # add measurements as text
    ax.text(x = 950,
            y = -100,
            s = f'freq: {round(freq, 3)} Hz\ninst freq: {round(inst_freq, 3)} Hz\ninit inst freq: {round(init_inst_freq, 3)} Hz',
            va = 'bottom',
            ha = 'left',
            fontsize = 8)
    
    # detected spikes
    ax.eventplot(t_spikes+250, 
                 orientation = 'horizontal', 
                 lineoffsets = 58, 
                 linelengths= 4,
                 linewidths = 1,
                 color = 'r')
    
    # edit axis
    apply_axis_settings(ax, axis = 'y', **ydict_v)
    
    # x
    xdict = {'ax_min' : 250,
             'ax_max' : 1250,
             'pad' : 20,
             'step' : 250,
             'stepminor' : 50,
             'label' : 'Time [ms]'}
    
    apply_axis_settings(ax, axis = 'x', **xdict)
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # display figure
    plt.show()
    
    # # create saving path and save
    # from parameters.directories_win import vplot_dir
    # path_fig = join(vplot_dir, 'cc_IF', 'traces')
    # save_figures(fig, f'{cell_ID}-cc_IF', path_fig, darkmode_bool, figure_format='png')
    


        
# %% rheobase spike

def plot_rheobase(cell_ID, idx_rheo, t, v, dvdt, rheospike_t, rheospike_v, rheospike_dvdt, rheospike_params):
    '''
    This function creates a figure displaying the rheobase step of the IF 
    protocol and marks the first spike. 
    Parameters:
        cell_ID: str, like 'E-303', unifque cell identifier
        idx_rheo: int, step index of the rheobase step
        t: numpy array, time dimension
        v: numpy array, voltage trace for single step
        dvdt: numpy array, first derivative of voltage trace for single step
        rheospike_t: numpy array, time dimension of rheobase spikes
        rheospike_v: numpy array, voltage trace of rheobase spikes
        rheospike_dvdt: numpy array, first derivative of voltage trace of rheobase spikes
        rheospike_params: dataframe, measurements of rheobase spike
    '''
    

    # init figure and axes
    fig, axs = plt.subplots(nrows = 1,
                            ncols = 2,
                            layout = 'constrained',
                            figsize = get_figure_size(),
                            dpi = 300)
    
    # set axis
    ax = axs[0]
    
    # set axis title
    ax.set_title(f'{cell_ID} cc_IF rheobase step no. {idx_rheo}',
                  fontsize=9, 
                  loc='left',
                  x = 0.02)
    
    # plot full trace
    ax.plot(t, v,
            color = colors_dict['primecolor'],
            lw = 0.5)
    
    # plot rheobase spike
    ax.plot(rheospike_t, rheospike_v,
            color = '#FFEC9DFF',
            lw = 0.5)
    
    # trace axis
    apply_axis_settings(ax, axis = 'y', **ydict_v)
    apply_axis_settings(ax, axis = 'x', **xdict_tfull)
    
    
    # inset marker
    box_xmin = rheospike_params.at[0, 't_peaks'] - 5
    box_width = 20
    box_ymin = -90
    box_height = 150
    
    # add rectangle marker
    ax.add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                            width = box_width, 
                            height = box_height,
                            fill = False,
                            color = colors_dict['primecolor'],
                            linestyle = '--',
                            lw = 0.5))
    
    
    # add inset
    ## ([left, bottom, width, height]), percentages
    ax_inset = ax.inset_axes([0.75, 0.65, 0.22, 0.30])
    
    # edit linewidth of inset axis and its ticks
    [ax_inset.spines[spine].set_linewidth(0.5) for spine in ['left', 'bottom']]
    ax_inset.tick_params(width=0.5)
    ax_inset.tick_params(which = 'minor', width=0.25)
    
    
    # plot full trace
    ax_inset.plot(t, v,
                  color = colors_dict['primecolor'],
                  lw = 0.5)
    
    # plot rheobase spike
    ax_inset.plot(rheospike_t, rheospike_v,
                  color = '#FFEC9DFF',
                  lw = 0.5)
    
    # x
    ax_inset.set_xticks(ticks = np.arange(0, 1500, 250), labels = [])
    ax_inset.set_xticks(ticks = np.arange(0, 1500, 5), labels = [], minor = True)
    ax_inset.set_xlim([box_xmin, box_xmin + box_width])
    
    # y
    ax_inset.set_yticks(ticks = np.arange(-100, 60+5, 50), labels = [])
    ax_inset.set_yticks(ticks = np.arange(-100, 60, 5), labels = [], minor = True)
    ax_inset.set_ylim([box_ymin, box_ymin + box_height])
    
    # # # dvdt # # #
    
    # set axis
    ax = axs[1]
    
    # set axis title
    ax.set_title(f'{cell_ID} cc_IF rheobase phase plane',
                  fontsize=9, 
                  loc='left',
                  x = 0.02)
    
    # plot full trace
    ax.plot(v, dvdt,
            color = colors_dict['primecolor'],
            lw = 0.5)
    
    # plot rheobase spike
    ax.plot(rheospike_v, rheospike_dvdt,
            color = '#FFEC9DFF',
            lw = 0.5)
    
    # edit axes
    # x
    apply_axis_settings(ax, axis = 'x', **ydict_v)
    
    # y
    dvdt_dict = {'ax_min' : -150,
                  'ax_max' : 250,
                  'pad' : 4,
                  'step' : 50,
                  'stepminor' : 10,
                  'label' : 'Rate of membrane\npotential change [mV/ms]'}
    
    apply_axis_settings(ax, axis = 'y', **dvdt_dict)
    
    # inset marker
    box_xmin   = -75
    box_width  = 30 
    box_ymin   = -10
    box_height = 20
    
    # add rectangle marker
    ax.add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                            width = box_width, 
                            height = box_height,
                            fill = False,
                            color = colors_dict['primecolor'],
                            linestyle = '--',
                            lw = 0.5))
    
    ## ([left, bottom, width, height]), percentages
    ax_inset2 = fig.add_axes([0.632, 0.79, 0.13, 0.145])
    
    
    # plot full trace
    ax_inset2.plot(v, dvdt,
                  color = colors_dict['primecolor'],
                  lw = 0.5)
    
    # plot rheobase spike
    ax_inset2.plot(rheospike_v, rheospike_dvdt,
                    color = '#FFEC9DFF',
                    lw = 1)
    
    # x
    ax_inset2.set_xticks(ticks = np.arange(box_xmin, box_xmin + box_width, 20), labels = [])
    ax_inset2.set_xticks(ticks = np.arange(box_xmin, box_xmin + box_width, 5), labels = [], minor = True)
    ax_inset2.set_xlim([box_xmin, box_xmin + box_width])
    
    # y
    ax_inset2.set_yticks(ticks = np.arange(box_ymin, box_ymin + box_height, 10), labels = [])
    ax_inset2.set_yticks(ticks = np.arange(box_ymin, box_ymin + box_height, 2), labels = [], minor = True)
    ax_inset2.set_ylim([box_ymin, box_ymin + box_height])
     
    # remove inset spines
    [ax_inset.spines[spine].set_visible(False) for spine in ['top', 'right']]
    [ax_inset2.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # remove spines
    [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]
    
    # create saving path and save
    from parameters.directories_win import vplot_dir
    path_fig = join(vplot_dir, 'cc_IF-rheobase_spike')
    save_figures(fig, f'{cell_ID}-cc_IF-rheobase_spike', path_fig, darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()
    
    
# %% adaptation

def plot_adaptation(cell_ID, t_full, v_full,
                    idx_rheo, idx_maxfreq, idx_halfmax, idx_maxinitinstfreq,
                    IF, IF_inst_init,
                    i_rheo_abs, i_maxfreq, i_halfmax, i_maxinitinstfreq,
                    n_lastspikes,
                    adaptation_spikes,
                    rheobase_spikes, maxfreq_spikes, halfmax_spikes, maxinitinstfreq_spikes,
                    adaptation_ISIs, 
                    rheobase_ISIs, maxfreq_ISIs, halfmax_ISIs, maxinitinstfreq_ISIs, 
                    fst_ISI, lst_ISIs, lst_inst_freqs, freq_adaptation_ratio, popt_adapfreq, t_linfit,  freq_adaptation_incline_linearfit,
                    fst_spike_vamplitude, lst_spike_vamplitude, spike_amplitude_adaptation,
                    fst_spike_FWHM, lst_spike_FWHM, spike_FWHM_adaptation):
    '''
    This function creates the adaptation figure for the cc_IF protocol analysis.
    Parameters:
        cell_ID
        t_full
        v_full
        idx_rheo
        idx_maxfreq
        idx_halfmax
        idx_maxinitinstfreq
        IF
        IF_inst_init
        i_rheo_abs
        i_maxfreq
        i_halfmax
        i_maxinitinstfreq
        n_lastspikes
        adaptation_spikes
        rheobase_spikes
        maxfreq_spikes
        halfmax_spikes
        maxinitinstfreq_spikes
        adaptation_ISIs
        rheobase_ISIs
        maxfreq_ISIs
        halfmax_ISIs
        maxinitinstfreq_ISIs
        fst_ISI
        lst_ISIs
        lst_inst_freqs
        freq_adaptation_ratio
        popt_adapfreq
        t_linfit
        freq_adaptation_incline_linearfit
        fst_spike_vamplitude
        lst_spike_vamplitude
        spike_amplitude_adaptation
        fst_spike_FWHM
        lst_spike_FWHM
        spike_FWHM_adaptation
    '''
    
    from functions.functions_useful import calc_dvdt_padded,linear_func
    from parameters.parameters import adaptation_n_lastspikes, adaptation_popt_guess_linear_ISIs
    
    # Exter
    colors = ['#FFEC9DFF', '#FAC881FF', '#F4A464FF', '#E87444FF', '#D9402AFF',
              '#BF2729FF', '#912534FF', '#64243EFF', '#3D1B28FF', '#161212FF']
    
    ax_keys = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']
    
    fig, axs = plt.subplot_mosaic('ABIK;CDIK;EFJL;GHJL',
                                  layout='tight',
                                  figsize=get_figure_size(),
                                  width_ratios=[3, 1.2, 3, 3],
                                  dpi = 300)
    
    # set figure title
    fig.suptitle(f'{cell_ID} frequency adaptation')
    
    # set titles
    subtitles_dict = {'fontsize' : 8,
                      'loc' : 'left'}
    
    
    # # # traces # # # 
    # rheobase
    axs['A'].set_title('A: Rheobase', **subtitles_dict)
    axs['A'].plot(t_full, v_full[idx_rheo],
                  c=colors[0],
                  lw=0.75)
    
    # max freq
    axs['C'].set_title('C: Max frequency (number of spikes)', **subtitles_dict)
    axs['C'].plot(t_full, v_full[idx_maxfreq],
                  c=colors[1],
                  lw=0.75)
    
    # half max freq
    axs['E'].set_title('E: Halfmax frequency', **subtitles_dict)
    axs['E'].plot(t_full, v_full[idx_halfmax],
                  c=colors[2],
                  lw=0.75)
    
    # max freq
    axs['G'].set_title('G: Max initial instantaneous frequency', **subtitles_dict)
    axs['G'].plot(t_full, v_full[idx_maxinitinstfreq],
                  c=colors[3],
                  lw=0.75)
    
    # edit axis
    # x
    xdict_t = {'ax_min' : 0,
               'ax_max' : 1500,
               'pad' : 10,
               'step' : 250,
               'stepminor' : 50,
               'label' : ''}
    
    ydict = {'ax_min' : -100,
              'ax_max' : 75,
              'pad' : 1.75,
              'step' : 100,
              'stepminor' : 25,
              'label' : '',
              'limits_n_0' : True}
        
    
    for ax_key in ['A','C','E','G']:
        # trace axis
        apply_axis_settings(axs[ax_key], axis = 'y', **ydict)
        
        # x
        apply_axis_settings(axs[ax_key], axis = 'x', **xdict_t)
        
    # remove axis 
    for ax_key in ['A','C','E']:
        axs[ax_key].set_xticks([])
        remove_spines_n_ticks([axs[ax_key]], axis = 'x')
    
    axs['G'].set_xlabel('Time [ms]')
    
    # set y label for traces
    fig.supylabel('Membrane potential [mV]')
    
    
    # # # phase plane # # #
    # rheobase, maxfreq, halfmax, max init inst freq
    for ax_key, c_idx, freq_idx in zip(['B', 'D', 'F', 'H'], [0, 1, 2, 3], [idx_rheo, idx_maxfreq, idx_halfmax, idx_maxinitinstfreq]):
        
        # set title
        axs[ax_key].set_title(f'{ax_key}: phase plane', **subtitles_dict)
        
        # plot
        axs[ax_key].plot(v_full[freq_idx], calc_dvdt_padded(v_full[freq_idx], t_full),
                      c=colors[c_idx],
                      lw=0.75)
    
    # set dict for phase plane y axis
    dvdt_dict = {'ax_min' : -150,
                  'ax_max' : 250,
                  'pad' : 4,
                  'step' : 250,
                  'stepminor' : 20,
                  'label' : '',
                  'limits_n_0' : True}
    
    # edit axis
    for ax_key in ['B','D','F','H']:
        
        # trace axis
        apply_axis_settings(axs[ax_key], axis = 'y', **dvdt_dict)
    
    # remove axis 
    for ax_key in ['B','D','F']:
        axs[ax_key].set_xticks([])
        remove_spines_n_ticks([axs[ax_key]], axis = 'x')
    
    # trace axis
    apply_axis_settings(axs['H'], axis = 'x', **ydict)
    
    
    # define plotting dicts
    plot_dict = {'marker' : '.', 
                 'markersize' : 6, 
                 'ls' : '-', 
                 'lw' : 0.75,
                 'markerfacecolor' : 'k',
                 'markeredgewidth' : '0.75'}
    
    fit_dict = {'c' : 'grey', 
                'ls' : '--', 
                'lw' : 1}
    
    lines_dict = {'colors' : colors_dict['primecolor'], 
                  'linestyle' : '--', 'lw' : 1}
    
    
    # # # IF curves # # #  
    
    # set axis
    ax = axs['I']
    
    # set title
    ax.set_title('I: Input current - frequency', **subtitles_dict)
    
    # IF curves
    ax.plot(IF[cell_ID].dropna(), c=colors[1], lw = 0.75)
    ax.plot(IF_inst_init[cell_ID].dropna(), c=colors[3], lw = 0.75)
    
    # mark previouse steps in IF curves
    ax.arrow(x=i_rheo_abs       , y=-5, dx=0, dy=4, color=colors[0], lw = 0.5)
    ax.arrow(x=i_maxfreq        , y=-5, dx=0, dy=4, color=colors[1], lw = 0.5)
    ax.arrow(x=i_halfmax        , y=-5, dx=0, dy=4, color=colors[2], lw = 0.5)
    ax.arrow(x=i_maxinitinstfreq, y=-5, dx=0, dy=4, color=colors[3], lw = 0.5)
    
    # edit axis
    iinput_dict = {'ax_min' : -50,
                   'ax_max' : 1000,
                   'pad' : 10,
                   'step' : 200,
                   'stepminor' : 50,
                   'label' : 'Input current [pA]',
                   'start_at_0' : True}
    
    apply_axis_settings(ax, axis = 'x', **iinput_dict)
    
    freq_dict = {'ax_min' : -10,
                 'ax_max' : 140,
                 'pad' : 1,
                 'step' : 20,
                 'stepminor' : 5,
                 'label' : 'Firing frequency [Hz]',
                 'start_at_0' : True}
    
    apply_axis_settings(ax, axis = 'y', **freq_dict)
    
    
    # # # spike frequency adaptation # # #
    
    # set title
    axs['K'].set_title('K: Spike frequency adaptation', **subtitles_dict)
    
    # plot
    for idx_freq, ISIs in enumerate([rheobase_ISIs, maxfreq_ISIs, halfmax_ISIs, maxinitinstfreq_ISIs]):
        
        axs['K'].plot(ISIs['t_ISI'], ISIs['inst_freq'],
                      c = colors[idx_freq],
                      **plot_dict)
    
        
        
    # plot measurements
    
    
    if freq_adaptation_ratio != 1.0:
        # line for first AP
        axs['K'].hlines(y = fst_ISI,
                        xmin = adaptation_ISIs.at[1, 't_ISI'],
                        xmax = adaptation_ISIs.tail(n_lastspikes)['t_ISI'].iloc[-1],
                        **lines_dict)
        
        # vertical line connecting first and last
        axs['K'].vlines(x = adaptation_ISIs.tail(n_lastspikes)['t_ISI'].iloc[-1],
                        ymin = lst_ISIs,
                        ymax = fst_ISI,
                        **lines_dict)
        
        # line for average of last APs
        axs['K'].hlines(y = lst_ISIs,
                        xmin = adaptation_ISIs.tail(n_lastspikes)['t_ISI'].iloc[0],
                        xmax = adaptation_ISIs.tail(n_lastspikes)['t_ISI'].iloc[-1],
                        **lines_dict)
    
    # plot fit
    if len(lst_inst_freqs) > 3:
        axs['K'].plot(t_linfit, linear_func(t_linfit, *popt_adapfreq), **fit_dict)
    
    # edit axis
    apply_axis_settings(axs['K'], axis = 'y', **freq_dict)
    
    xdict_t = {'ax_min' : 250,
               'ax_max' : 1250,
               'pad' : 15,
               'step' : 250,
               'stepminor' : 50,
               'label' : 'Time [ms]'}
    
    for ax_key in ['K','J','L']:
        # x
        apply_axis_settings(axs[ax_key], axis = 'x', **xdict_t)
        
    # add text
    axs['K'].text(x = xdict_t['ax_max'],
                  y = freq_dict['ax_min'],
                  s = f'freq. adap.: {round(freq_adaptation_ratio, 2)}\nfreq. incline: {round(freq_adaptation_incline_linearfit, 2)} Hz/s',
                  ha = 'right', va = 'bottom',
                  fontsize = 8)
        
    
    # # # spike amplitude adaptation # # #
    
    # set title
    axs['J'].set_title('J: Spike amplitude adaptation', **subtitles_dict)
    
    # plot
    for idx_freq, spikes in enumerate([rheobase_spikes, maxfreq_spikes, halfmax_spikes, maxinitinstfreq_spikes]):
        
        axs['J'].plot(spikes['t_peaks'], spikes['v_amplitude'],
                      c = colors[idx_freq],
                      **plot_dict)
    
    # plot measurements
    # line for first AP
    axs['J'].hlines(y = fst_spike_vamplitude,
                    xmin = adaptation_spikes.at[0, 't_peaks'],
                    xmax = adaptation_spikes.tail(n_lastspikes)['t_peaks'].iloc[-1],
                    **lines_dict)
    
    # vertical line connecting first and last
    axs['J'].vlines(x = adaptation_spikes.tail(n_lastspikes)['t_peaks'].iloc[-1],
                    ymin = lst_spike_vamplitude,
                    ymax = fst_spike_vamplitude,
                    **lines_dict)
    
    # line for average of last APs
    axs['J'].hlines(y = lst_spike_vamplitude,
                    xmin = adaptation_spikes.tail(n_lastspikes)['t_peaks'].iloc[0],
                    xmax = adaptation_spikes.tail(n_lastspikes)['t_peaks'].iloc[-1],
                    **lines_dict)
    
    # edit axis
    ampl_dict = {'ax_min' : 0,
                 'ax_max' : 140,
                 'pad' : 1.4,
                 'step' : 40,
                 'stepminor' : 5,
                 'label' : 'Spike amplitude [mV]'}
    
    apply_axis_settings(axs['J'], axis = 'y', **ampl_dict)
    
    # add text
    axs['J'].text(x = xdict_t['ax_max'],
                  y = ampl_dict['ax_min'],
                  s = "spike ampl. adap.: %.2f" % spike_amplitude_adaptation,
                  ha = 'right', va = 'bottom',
                  fontsize = 8)
    
    
    # # # spike FWHM adaptation # # #
    
    # set title
    axs['L'].set_title('L: Spike width adaptation', **subtitles_dict)
    
    # plot
    for idx_freq, spikes in enumerate([rheobase_spikes, maxfreq_spikes, halfmax_spikes, maxinitinstfreq_spikes]):
        
        axs['L'].plot(spikes['t_peaks'], spikes['FWHM'],
                      c = colors[idx_freq],
                      **plot_dict)
    
    # plot measurements
    # line for first AP
    axs['L'].hlines(y = fst_spike_FWHM,
                    xmin = adaptation_spikes.at[0, 't_peaks'],
                    xmax =  adaptation_spikes.tail(n_lastspikes)['t_peaks'].iloc[-1],
                    **lines_dict)
    
    # vertical line connecting first and last
    axs['L'].vlines(x =  adaptation_spikes.tail(n_lastspikes)['t_peaks'].iloc[-1],
                    ymin = lst_spike_FWHM,
                    ymax = fst_spike_FWHM,
                    **lines_dict)
    
    # line for average of last APs
    axs['L'].hlines(y = lst_spike_FWHM,
                    xmin = adaptation_spikes.tail(n_lastspikes)['t_peaks'].iloc[0],
                    xmax = adaptation_spikes.tail(n_lastspikes)['t_peaks'].iloc[-1],
                    **lines_dict)
    
    # edit axis
    FWHM_dict = {'ax_min' : 0,
                 'ax_max' : 7,
                 'pad' : 0.07,
                 'step' : 1,
                 'stepminor' : 0.5,
                 'label' : 'Spike FWHM [ms]'}
    
    apply_axis_settings(axs['L'], axis = 'y', **FWHM_dict)
    
    # add text
    axs['L'].text(x = xdict_t['ax_max'],
                  y = FWHM_dict['ax_min'],
                  s = "FWHM. adap.: %.2f" % spike_FWHM_adaptation,
                  ha = 'right', va = 'bottom',
                  fontsize = 8)
    
    # remove spines
    [axs[ax_key].spines[spine].set_visible(False) for ax_key in ax_keys for spine in ['top', 'right']]
    
    # align labels
    fig.align_labels()
    
    # create saving path and save
    from parameters.directories_win import vplot_dir
    path_fig = join(vplot_dir, 'cc_IF-adaptation')
    save_figures(fig, f'{cell_ID}-cc_IF-adaptation', path_fig, darkmode_bool, figure_format='png')
        
    # display figure
    plt.show()
    
