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
           'pad' : 10,
           'step' : 20,
           'stepminor' : 5,
           'label' : 'Membrane\npotential [mV]'}

# x
xdict_tfull = {'ax_min' : 0,
               'ax_max' : 1500,
               'pad' : 10,
               'step' : 250,
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
    path_fig = join(vplot_dir, 'cc_IF', 'traces')
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
                            figsize = get_figure_size(width = 200, height = 100))
    
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
            color = colors_dict['color2'],
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
    ax_inset = fig.add_axes([0.385, 0.59, 0.10, 0.3])
    
    
    # plot full trace
    ax_inset.plot(t, v,
                  color = colors_dict['primecolor'],
                  lw = 0.5)
    
    # plot rheobase spike
    ax_inset.plot(rheospike_t, rheospike_v,
                  color = colors_dict['color2'],
                  lw = 0.5)
    
    # x
    ax_inset.set_xticks(ticks = np.arange(0, 1500, 250), labels = [])
    ax_inset.set_xticks(ticks = np.arange(0, 1500, 5), labels = [], minor = True)
    ax_inset.set_xlim([box_xmin, box_xmin + box_width])
    
    # y
    ax_inset.set_yticks(ticks = np.arange(-100, 60+5, 20), labels = [])
    ax_inset.set_yticks(ticks = np.arange(-100, 60, 5), labels = [], minor = True)
    ax_inset.set_ylim([box_ymin, box_ymin + box_height])
    
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
            color = colors_dict['color2'],
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
    
    # remove inset spines
    [ax_inset.spines[spine].set_visible(False) for spine in ['top', 'right']]
    
    # remove spines
    [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

    # create saving path and save
    from parameters.directories_win import vplot_dir
    path_fig = join(vplot_dir, 'cc_IF', 'rheobase_1stAP')
    save_figures(fig, f'{cell_ID}-cc_IF-rheobase_1stAP', path_fig, darkmode_bool, figure_format='png')
    
    # display figure
    plt.show()
    
    
# %% adaptation

'''
Parameters:
    t_full
    v_full
    idx_rheo,
    idx_maxfreq
    idx_halfmax
    idx_maxinitinstfreq
    IF
    IF_initinstfreq
    i_rheo
    i_maxfreq
    i_halfmax
    i_maxinstinitialfreq
'''

from functions.functions_useful import calc_dvdt_padded

# Exter
colors = ['#FFEC9DFF', '#FAC881FF', '#F4A464FF', '#E87444FF', '#D9402AFF',
          '#BF2729FF', '#912534FF', '#64243EFF', '#3D1B28FF', '#161212FF']

ax_keys = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']

fig, axs = plt.subplot_mosaic('ABIK;CDIK;EFJL;GHJL',
                              layout='tight',
                              figsize=get_figure_size(),
                              width_ratios=[3, 1.2, 3, 3]
                              )

# set figure title
fig.suptitle(f'{cell_ID} frequency adaptation')


# # # traces # # # 
# rheobase
axs['A'].set_title('A: Rheobase', fontsize=8, loc='left')
axs['A'].plot(t_full, v_full[idx_rheo],
              c=colors[0],
              lw=0.75)

# max freq
axs['C'].set_title('C: Max frequency (number of spikes)', fontsize=8, loc='left')
axs['C'].plot(t_full, v_full[idx_maxfreq],
              c=colors[1],
              lw=0.75)

# half max freq
axs['E'].set_title('E: Halfmax frequency', fontsize=8, loc='left')
axs['E'].plot(t_full, v_full[idx_halfmax],
              c=colors[2],
              lw=0.75)

# max freq
axs['G'].set_title('G: Max initial instantaneous frequency', fontsize=8, loc='left')
axs['G'].plot(t_full, v_full[idx_maxinitinstfreq],
              c=colors[3],
              lw=0.75)

# edit axis
for ax_key in ['A','C','E','G']:
    
    ydict = {'ax_min' : -100,
              'ax_max' : 75,
              'pad' : 1.75,
              'step' : 100,
              'stepminor' : 25,
              'label' : '',
              'limits_n_0' : True}
    
    # trace axis
    apply_axis_settings(axs[ax_key], axis = 'y', **ydict)
    
# remove axis 
for ax_key in ['A','C','E']:
    axs[ax_key].set_xticks([])
    remove_spines_n_ticks([axs[ax_key]], axis = 'x')

# x
apply_axis_settings(axs['G'], axis = 'x', **xdict_tfull)

# set y label for traces
fig.supylabel('Membrane potential [mV]')


# # # phase plane # # #
# rheobase
axs['B'].plot(v_full[idx_rheo], calc_dvdt_padded(v_full[idx_rheo], t_full),
              c=colors[0],
              lw=0.75)

# maxfreq
axs['D'].plot(v_full[idx_maxfreq], calc_dvdt_padded(v_full[idx_maxfreq], t_full),
              c=colors[1],
              lw=0.75)

# halfmax
axs['F'].plot(v_full[idx_halfmax], calc_dvdt_padded(v_full[idx_halfmax], t_full),
              c=colors[2],
              lw=0.75)

# max init inst freq
axs['H'].plot(v_full[idx_maxinitinstfreq], calc_dvdt_padded(v_full[idx_maxinitinstfreq], t_full),
              c=colors[3],
              lw=0.75)

# edit axis
for ax_key in ['B','D','F','H']:
    
    dvdt_dict = {'ax_min' : -150,
                  'ax_max' : 250,
                  'pad' : 4,
                  'step' : 250,
                  'stepminor' : 20,
                  'label' : '',
                  'limits_n_0' : True}
    
    # trace axis
    apply_axis_settings(axs[ax_key], axis = 'y', **dvdt_dict)

# remove axis 
for ax_key in ['B','D','F']:
    axs[ax_key].set_xticks([])
    remove_spines_n_ticks([axs[ax_key]], axis = 'x')

# trace axis
apply_axis_settings(axs['H'], axis = 'x', **ydict)


# # # IF curves # # #  

# set axis
ax = axs['I']

# IF curves
ax.plot(IF[cell_ID].dropna(), c=colors[1])
ax.plot(IF_inst_init[cell_ID].dropna(), c=colors[3])

# mark previouse steps in IF curves
ax.arrow(x=i_rheo_abs       , y=-5, dx=0, dy=4, color=colors[0])
ax.arrow(x=i_maxfreq        , y=-5, dx=0, dy=4, color=colors[1])
ax.arrow(x=i_halfmax        , y=-5, dx=0, dy=4, color=colors[2])
ax.arrow(x=i_maxinitinstfreq, y=-5, dx=0, dy=4, color=colors[3])

plot_dict = {'marker' : '.', 
              'markersize' : 5, 
              'ls' : '-', 
              'lw' : 1}

fit_dict = {'c' : 'grey', 
            'ls' : '--', 
            'lw' : 1}

lines_dict = {'colors' : colors_dict['primecolor'], 
              'linestyle' : '--', 'lw' : 1}

iinput_dict = {'ax_min' : -50,
                'ax_max' : 1000,
                'pad' : 10,
                'step' : 200,
                'stepminor' : 50,
                'label' : 'Input current [pA]',
                'start_at_0' : True}

apply_axis_settings(ax, axis = 'x', **iinput_dict)


freq_dict = {'ax_min' : -10,
              'ax_max' : 100,
              'pad' : 1,
              'step' : 20,
              'stepminor' : 5,
              'label' : 'Firing frequency [Hz]',
              'start_at_0' : True}

apply_axis_settings(ax, axis = 'y', **freq_dict)


# remove spines
[axs[ax_key].spines[spine].set_visible(False) for ax_key in ax_keys for spine in ['top', 'right']]
    
# display figure
plt.show()






# %%




# # max freq
# axs['C'].plot(t_step, v_maxfreq,
#               c=colors[1],
#               lw=1)

# axs['D'].plot(v_maxfreq, dvdt_maxfreq,
#               c=colors[1],
#               lw=1)

# # half max freq
# axs['E'].plot(t_step, v_halfmax,
#               c=colors[2],
#               lw=1)

# axs['F'].plot(v_halfmax, dvdt_halfmax,
#               c=colors[2],
#               lw=1)

# # max inst initial freq
# axs['G'].plot(t_step, v_maxinstinitialfreq,
#               c=colors[3],
#               lw=1)

# axs['H'].plot(v_maxinstinitialfreq, dvdt_maxinstinitialfreq,
#               c=colors[3],
#               lw=1)



# IF curves
axs['I'].plot(IF_df[cell_ID], c=colors[1])
axs['I'].plot(IF_inst_initial_df[cell_ID], c=colors[3])

# mark previouse steps in IF curves
axs['I'].arrow(x=i_rheobase, y=-5, dx=0, dy=4, color=colors[0])
axs['I'].arrow(x=i_maxfreq, y=-5, dx=0, dy=4, color=colors[1])
axs['I'].arrow(x=i_halfmax, y=-5, dx=0, dy=4, color=colors[2])
axs['I'].arrow(x=i_maxinstinitialfreq, y=-5, dx=0, dy=4, color=colors[3])

plot_dict = {'marker' : '.', 'markersize' : 5, 'ls' : '-', 'lw' : 1}
fit_dict = {'c' : 'grey', 'ls' : '--', 'lw' : 1}
lines_dict = {'colors' : colors_dict['primecolor'], 'linestyle' : '--', 'lw' : 1}


# v_amplitude
axs['J'].plot(rheobase_spikes['t_peaks'], rheobase_spikes['v_amplitude'], 
              c=colors[0], **plot_dict)
axs['J'].plot(maxfreq_spikes['t_peaks'], maxfreq_spikes['v_amplitude'], 
              c=colors[1], **plot_dict)
axs['J'].plot(halfmax_spikes['t_peaks'], halfmax_spikes['v_amplitude'], 
              c=colors[2], **plot_dict)
axs['J'].plot(maxinstinitialfreq_spikes['t_peaks'], maxinstinitialfreq_spikes['v_amplitude'], 
              c=colors[3], **plot_dict)

axs['J'].hlines(y = lst_spike_vamplitude,
                xmin = adaptation_spikes.at[len(adaptation_spikes)-n_last_spikes, 't_peaks'],
                xmax = adaptation_spikes.at[len(adaptation_spikes)-1, 't_peaks'],
                **lines_dict)

axs['J'].vlines(x = adaptation_spikes.at[int(len(adaptation_spikes)-n_last_spikes/2), 't_peaks'],
                ymin = lst_spike_vamplitude,
                ymax = fst_spike_vamplitude,
                **lines_dict)

axs['J'].hlines(y = fst_spike_vamplitude,
                xmin = adaptation_spikes.at[0, 't_peaks'],
                xmax = adaptation_spikes.at[len(adaptation_spikes)-1, 't_peaks'],
                **lines_dict)

# ISI
axs['K'].plot(rheobase_ISIs['t_ISI'], rheobase_ISIs['inst_freq'], 
              c=colors[0], **plot_dict)
axs['K'].plot(maxfreq_ISIs['t_ISI'], maxfreq_ISIs['inst_freq'], 
              c=colors[1], **plot_dict)
axs['K'].plot(halfmax_ISIs['t_ISI'], halfmax_ISIs['inst_freq'], 
              c=colors[2], **plot_dict)
axs['K'].plot(maxinstinitialfreq_ISIs['t_ISI'], maxinstinitialfreq_ISIs['inst_freq'], 
              c=colors[3], **plot_dict)

axs['K'].hlines(y = lst_ISIs,
                xmin = adaptation_ISIs.at[len(adaptation_ISIs)-n_last_ISIs, 't_ISI'],
                xmax = adaptation_ISIs.at[len(adaptation_ISIs)-1, 't_ISI'],
                **lines_dict)

axs['K'].vlines(x = adaptation_ISIs.at[int(len(adaptation_ISIs)-n_last_ISIs/2), 't_ISI'],
                ymin = lst_ISIs,
                ymax = fst_ISI,
                **lines_dict)

axs['K'].hlines(y = fst_ISI,
                xmin = adaptation_ISIs.at[1, 't_ISI'],
                xmax = adaptation_ISIs.at[len(adaptation_ISIs)-1, 't_ISI'],
                **lines_dict)

if len(lst_inst_freqs) > 3:
    axs['K'].plot(t_linfit, linear_func(t_linfit, *popt), **fit_dict)


# FWHM
axs['L'].plot(rheobase_spikes['t_peaks'], rheobase_spikes['FWHM'], 
              c=colors[0], **plot_dict)
axs['L'].plot(maxfreq_spikes['t_peaks'], maxfreq_spikes['FWHM'], 
              c=colors[1], **plot_dict)
axs['L'].plot(halfmax_spikes['t_peaks'], halfmax_spikes['FWHM'], 
              c=colors[2], **plot_dict)
axs['L'].plot(maxinstinitialfreq_spikes['t_peaks'], maxinstinitialfreq_spikes['FWHM'], 
              c=colors[3], **plot_dict)


axs['L'].hlines(y = lst_spike_FWHM,
                xmin = adaptation_spikes.at[len(adaptation_spikes)-n_last_spikes, 't_peaks'],
                xmax = adaptation_spikes.at[len(adaptation_spikes)-1, 't_peaks'],
                **lines_dict)

axs['L'].vlines(x = adaptation_spikes.at[int(len(adaptation_spikes)-n_last_spikes/2), 't_peaks'],
                ymin = lst_spike_FWHM,
                ymax = fst_spike_FWHM,
                **lines_dict)

axs['L'].hlines(y = fst_spike_FWHM,
                xmin = adaptation_spikes.at[0, 't_peaks'],
                xmax = adaptation_spikes.at[len(adaptation_spikes)-1, 't_peaks'],
                **lines_dict)



# format axis
v_range = [-100, 75]

axs['A'].set_title('A: Rheobase', fontsize='small', loc='left')
axs['C'].set_title('C: Max frequency (number of spikes)', fontsize='small', loc='left')
axs['E'].set_title('E: Halfmax frequency', fontsize='small', loc='left')
axs['G'].set_title('G: Max initial instantaneous spiking frequency', fontsize='small', loc='left')

for ax_idx in ['A', 'C', 'E']:
    # x
    axs[ax_idx].set_xlim([0, 1500])
    axs[ax_idx].set_xticks(np.arange(0, 1500+1, 250), minor=True)
    axs[ax_idx].set_xticks(np.arange(0, 1500+1, 500), labels=[])

for ax_idx in ['G', 'J', 'K', 'L']:
    # x
    axs[ax_idx].set_xlim([0, 1500])
    axs[ax_idx].set_xticks(np.arange(0, 1500+1, 250), minor=True)
    axs[ax_idx].set_xticks(np.arange(0, 1500+1, 500))
    axs[ax_idx].set_xlabel('Time [ms]')

for ax_idx in ['A', 'C', 'E', 'G']:
    # y
    axs[ax_idx].set_ylim(v_range)
    axs[ax_idx].set_yticks(
        np.arange(v_range[0], v_range[1]+1, 25), minor=True)

fig.supylabel('Voltage [mV]')
axs['G'].set_xlabel('Time [ms]')

for ax_idx in ['B', 'D', 'F', 'H']:
    # x dvdt
    axs[ax_idx].set_title(f'{ax_idx}: phase plane', fontsize='small', loc='left')        
    axs[ax_idx].set_xlim(v_range)
    axs[ax_idx].set_xticks(np.arange(v_range[0], v_range[1] + 1, 100))
    axs[ax_idx].set_xticks(np.arange(v_range[0], v_range[1] + 1, 25), minor=True)
    # y dvdt
    # axs[ax_idx].set_ylabel('Rate of membrane potential change [mV/ms]')
    axs[ax_idx].set_ylim([-150, 250])
    axs[ax_idx].set_yticks(np.arange(0, 250 + 1, 200))
    axs[ax_idx].set_yticks(np.arange(-150, 250 + 1, 50), minor=True)

axs['H'].set_xlabel('Voltage [mV]')
# axs['F'].set_ylabel('Rate of membrane potential change [mV/ms]')


# # combine all spikes df to find min or max
# all_spikes = pd.concat([rheobase_spikes, maxfreq_spikes, halfmax_spikes, maxinstinitialfreq_spikes], axis = 0)
# all_ISIs = pd.concat([rheobase_ISIs, maxfreq_ISIs, halfmax_ISIs, maxinstinitialfreq_ISIs], axis = 0)

# ### IF axes ###
# IF_xmin = round_to_base(IF_df[cell_ID].dropna().index[0], 50)
# IF_xmax = round_up_to_base(IF_df[cell_ID].dropna().index[-1], 50)

# IF_ymax = round_up_to_base(active_properties_df.at[cell_ID, 'max_inst_initial_freq'], 20)

# axs['I'].set_title('I: Input current - frequency',
#                    fontsize='small', loc='left')
# # x
# axs['I'].set_xlabel('Input current [pA]')
# axs['I'].set_xlim([IF_xmin, IF_xmax])
# axs['I'].set_xticks(np.arange(0, IF_xmax+1, 100))
# axs['I'].set_xticks(np.arange(IF_xmin, IF_xmax+1, 25), minor = True)

# # y
# axs['I'].set_ylabel('Firing frequency [Hz]')
# axs['I'].set_ylim([-10, IF_ymax])
# axs['I'].set_yticks(np.arange(0, IF_ymax+1, 20))
# axs['I'].set_yticks(np.arange(0, IF_ymax+1, 5), minor = True)    


### v_amplitude ###
amp_ymax = round_up_to_base(all_spikes['v_amplitude'].max(), 10)
amp_ymin = round_down_to_base(all_spikes['v_amplitude'].min(), 10)

axs['J'].set_title('J: Spike amplitude adaptation',
                    fontsize='small', loc='left')

# y
axs['J'].set_ylabel('Spike amplitude [mV]')
axs['J'].set_ylim([amp_ymin, amp_ymax])
axs['J'].set_yticks(np.arange(amp_ymin, amp_ymax+1, 20))
axs['J'].set_yticks(np.arange(amp_ymin, amp_ymax+1, 10), minor = True) 

axs['J'].text(x = 1450,
              y = amp_ymin + (amp_ymax - amp_ymin)*0.05,
              s = "freq. adap.:\n %.2f" % spike_amplitude_adaptation,
              ha = 'right', va = 'bottom',
              fontsize = 8)


### ISI / Frequency ###
freq_ymax = round_up_to_base(all_ISIs['inst_freq'].max(), 10)
#freq_ymin = round_down_to_base(rheobase_ISIs['inst_freq'].min(), 10)
freq_ymin = 0

axs['K'].set_title('K: Spike frequency adaptation',
                    fontsize='small', loc='left')

    
axs['K'].text(x = 1450,
              y = freq_ymin + (freq_ymax - freq_ymin)*0.05,
              s = "freq. adap.: %.2f" % freq_adaptation_ratio + '\n' + 'freq. incline: %.2f' % freq_adaptation_incline_linearfit + ' Hz/s',
              ha = 'right', va = 'bottom',
              fontsize = 8)

# y
axs['K'].set_ylabel('instantaneous spike\nfrequency [Hz]')
axs['K'].set_ylim([freq_ymin, freq_ymax])
axs['K'].set_yticks(np.arange(freq_ymin, freq_ymax+1, 20))
axs['K'].set_yticks(np.arange(freq_ymin, freq_ymax+1, 10), minor = True) 


### FWHM ###
FWHM_ymax = round_up_to_base(all_spikes['FWHM'].max(), 1) + .5
FWHM_ymin = round_down_to_base(all_spikes['FWHM'].min(), 1)

axs['L'].set_title('L: Spike FWHM adaptation',
                    fontsize='small', loc='left')

axs['L'].text(x = 1450,
              y = 0.75,
              s = "FWHM. adap.:\n%.2f" % spike_FWHM_adaptation,
              ha = 'right', va = 'bottom',
              fontsize = 8)

# y
axs['L'].set_ylabel('Spike FHWM [ms]')
axs['L'].set_ylim([0.5, FWHM_ymax])
axs['L'].set_yticks(np.arange(1, FWHM_ymax+.1, 1))
axs['L'].set_yticks(np.arange(0.5, FWHM_ymax+.1, .25), minor = True)

fig.align_labels()

vplots_path_fig = join(vplot_dir, 'cc_IF', 'freq_adaptation')
save_figures(fig, f'{cell_ID}-frequency_adaptation', vplots_path_fig, darkmode_bool, figure_format='both')

plt.show()