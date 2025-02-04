# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:01:03 2024

@author: nesseler
"""

import warnings
warnings.warn('Script not up-to-date!')
# %%

from os.path import join
import pandas as pd
import matplotlib as mtl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import seaborn as sbn


from parameters.directories_win import table_file, quant_data_dir, cell_descrip_dir, vplot_dir

from parameters.PGFs import cc_IF_parameters
from parameters.parameters import min_peak_prominence_ccIF, min_max_peak_width_ccIF, min_peak_distance_ccIF, popt_guess_linear_ISIs

from functions.functions_useful import calc_time_series, calc_dvdt_padded, butter_filter, round_to_base, round_up_to_base, round_down_to_base, linear_func
from functions.functions_constructors import construct_current_array
from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_extractspike import get_AP_parameters
from getter.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes

vplot_bool = True
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)


# %% load passive properties

# cc_rest
activity_df = pd.read_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col='cell_ID')

# cc_IF
passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col='cell_ID')
active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col='cell_ID')
fstAP_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_col='cell_ID')
IF_step_idc = pd.read_excel(join(cell_descrip_dir, 'ccIF-step_indices.xlsx'), index_col='cell_ID')

IF_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_col='i_input')
IF_inst_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_col='i_input')
IF_inst_initial_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst_initial.xlsx'), index_col='i_input')


# calc delta v_thres to v_rest
delta_vrest_to_vthres_df = fstAP_df['v_threshold'] - activity_df['v_rest']

# %% analyzable cells

# protocol
PGF = 'cc_IF'

# get cell IDs
cell_IDs = get_cell_IDs_one_protocol(PGF)

# cell_IDs.reverse()

# cell_IDs = ['E-175']

for cell_ID in cell_IDs:
    # test if desired cells in passive properties
    if cell_ID not in passiv_properties_df.index.to_list():
        raise ValueError(
            f'Passiv properties not calculated for provided cell_ID: {cell_ID}!')

    # test if desired cells in passive properties
    if cell_ID not in activity_df.index.to_list():
        raise ValueError(
            f'Resting membrane properties not calculated for provided cell_ID: {cell_ID}!')

# get hold current as table
I_hold_table = pd.read_excel(
    table_file, sheet_name="V_or_I_hold", index_col='cell_ID').loc[cell_IDs, :]


# %% initialize export dataframes

adaptation_df = pd.DataFrame(index = cell_IDs)

# %% load steps

for cell_ID in cell_IDs:

    print(f'Started: {cell_ID}')

    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)

    # get IF data form file
    i, v, t, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')

    # sampling rate in ms
    SR_ms = SR / 1e3

    # concatenate individual steps
    n_points = int(np.shape(i)[0] * np.shape(i)[1])

    v_concat = v.flatten()

    t_ms = calc_time_series(v_concat, SR)
    t_s = calc_time_series(v_concat, SR, scale='s')

    # filter voltage (to vf)
    vf = butter_filter(v_concat,
                       order=3,
                       cutoff=1e3,
                       sampling_rate=SR)

    # construct current dataframe
    i_hold = I_hold_table.at[cell_ID, PGF]

    # calculate current steps relative to I_hold
    # rounded to nearest 5
    i_hold_rounded = round_to_base(i_hold, 5)

    # get current arrays and list of input current relative to i_hold
    i, i_input = construct_current_array(i_hold=i_hold_rounded,
                                         n_steps=n_steps,
                                         parameters_dict=cc_IF_parameters,
                                         SR_ms=SR_ms)

    step_dur = cc_IF_parameters['t_pre'] + \
        cc_IF_parameters['t_stim'] + cc_IF_parameters['t_post']
    step_points = step_dur * SR_ms

    dvdt = [None] * n_steps

    t_step = np.arange(step_dur, step=1/SR_ms)

    # loop through steps to limit voltage trace
    for step_idx in np.arange(0, n_steps, 1):

        # calc start and stop indices for step
        start_idx = int(step_points * step_idx)
        stop_idx = int(start_idx + step_points)

        # set voltage trace of step
        v_step = vf[start_idx:stop_idx]
        dvdt_step = calc_dvdt_padded(v_step, t_step)

        v[step_idx] = v_step
        dvdt[step_idx] = dvdt_step

    # %% step choice

    # rheobase
    idx_rheobase = IF_step_idc.at[cell_ID, 'rheobase_step_idx']
    v_rheobase = v[idx_rheobase]
    dvdt_rheobase = dvdt[idx_rheobase]
    i_rheobase = IF_df[cell_ID].dropna().index[idx_rheobase]

    # max freq
    # idx_maxfreq = IF_step_idc.at[cell_ID, 'maxfreq_step_idx']
    idx_maxfreq = IF_df[cell_ID].dropna().argmax()
    v_maxfreq = v[idx_maxfreq]
    dvdt_maxfreq = dvdt[idx_maxfreq]
    i_maxfreq = IF_df[cell_ID].dropna().index[idx_maxfreq]

    # max inst initial freq
    # idx_maxinstinitialfreq = IF_step_idc.at[cell_ID, 'maxinstinitialfreq_step_idx']
    idx_maxinstinitialfreq = IF_inst_initial_df[cell_ID].dropna().argmax()
    i_maxinstinitialfreq = IF_inst_initial_df[cell_ID].dropna().index[idx_maxinstinitialfreq]
    
    idx_step_maxinstinitial_freq = IF_df[cell_ID].dropna().index.to_list().index(i_maxinstinitialfreq)
    v_maxinstinitialfreq = v[idx_step_maxinstinitial_freq]
    dvdt_maxinstinitialfreq = dvdt[idx_step_maxinstinitial_freq]
    
    
    # step half-way between rheobase and max inst initial freq
    idx_step_halfmax = int(((idx_step_maxinstinitial_freq - idx_rheobase) / 2) + idx_rheobase)
    # get voltage & dvdt
    v_halfmax = v[idx_step_halfmax]
    i_halfmax = IF_df[cell_ID].dropna().index[idx_step_halfmax]
    dvdt_halfmax = dvdt[idx_step_halfmax]
    

    # %% spike adaptations

    def get_spike_params_df(t, v, dvdt):
        # find peaks
        idc_peaks, dict_peak = sc.signal.find_peaks(v,
                                                    prominence=min_peak_prominence_ccIF,
                                                    distance=min_peak_distance_ccIF *
                                                    (SR_ms),
                                                    width=np.multiply(min_max_peak_width_ccIF, SR_ms))

        # limit spike indices to stimulation time period
        pre_points = int(cc_IF_parameters['t_pre'] * SR_ms)
        pre_n_stim_points = int(
            (cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim']) * SR_ms)
        idc_peaks = [idx_peak for idx_peak in idc_peaks if idx_peak >
                     pre_points and idx_peak <= pre_n_stim_points]

        # get AP parameters of first spike
        spike_params_df, _ = get_AP_parameters(t, v, dvdt, idc_peaks)

        ISI_df = pd.DataFrame()
        ISI_df['ISI'] = spike_params_df['t_peaks'].diff()
        ISI_df['inst_freq'] = (1 / ISI_df['ISI']) * 1e3
        ISI_df['t_ISI'] = spike_params_df['t_peaks'].shift(
            1) + (spike_params_df['t_peaks'].diff() / 2)

        return spike_params_df, ISI_df

    # get all spike parameters
    rheobase_spikes, rheobase_ISIs = get_spike_params_df(t_step, v_rheobase, dvdt_rheobase)
    maxfreq_spikes, maxfreq_ISIs = get_spike_params_df(t_step, v_maxfreq, dvdt_maxfreq)
    halfmax_spikes, halfmax_ISIs = get_spike_params_df(t_step, v_halfmax, dvdt_halfmax)
    maxinstinitialfreq_spikes, maxinstinitialfreq_ISIs = get_spike_params_df(t_step, v_maxinstinitialfreq, dvdt_maxinstinitialfreq)


    # %% adaptation
    
    n_last_spikes = 4
    n_last_ISIs = n_last_spikes - 1
    
    ## spike amplitude
    # take max number of spikes
    adaptation_spikes = maxfreq_spikes
    
    # get amplitude of first spike
    fst_spike_vamplitude = adaptation_spikes.at[0, 'v_amplitude']
    
    # get amplitude of last 4 spikes & calc mean
    lst_spike_vamplitude = adaptation_spikes.tail(n_last_spikes)['v_amplitude'].mean()
    
    # calc ratio of first to last spikes
    spike_amplitude_adaptation = lst_spike_vamplitude / fst_spike_vamplitude
    
    
    ## spike FWHM
    fst_spike_FWHM = adaptation_spikes.at[0, 'FWHM']
    lst_spike_FWHM = adaptation_spikes.tail(n_last_spikes)['FWHM'].mean()
    spike_FWHM_adaptation = lst_spike_FWHM / fst_spike_FWHM
    

    ## spike frequency
    adaptation_ISIs = maxfreq_ISIs
    fst_ISI = adaptation_ISIs.at[1, 'inst_freq']
    lst_ISIs = adaptation_ISIs.tail(n_last_ISIs)['inst_freq'].mean()
    freq_adaptation_ratio = lst_ISIs / fst_ISI


    ## linear fit to ISIs between 500 and 1000 ms
    # limit dataframe to ISIs to ISIs later in step
    ISIs_for_linear_fit = adaptation_ISIs.query('t_ISI > 500')
    
    # linear fit
    lst_tISIs = ISIs_for_linear_fit['t_ISI']
    lst_inst_freqs = ISIs_for_linear_fit['inst_freq']
    
    # plt.scatter(adaptation_spikes['t_peaks'], adaptation_spikes['v_amplitude'])
    
    if len(lst_inst_freqs) > 3:
        # fit linear curve
        popt, pcov = sc.optimize.curve_fit(linear_func, lst_tISIs, lst_inst_freqs, 
                                           p0 = popt_guess_linear_ISIs, 
                                           maxfev = 5000)
        
        # plot exponential fit
        t_linfit = np.arange(start = ISIs_for_linear_fit['t_ISI'].iat[0], 
                             stop = ISIs_for_linear_fit['t_ISI'].iat[-1], 
                             step = 1)
        
        # get incline of linear fit
        freq_adaptation_incline_linearfit = popt[0] # Hz / ms
    
        # convert to Hz / s
        freq_adaptation_incline_linearfit = freq_adaptation_incline_linearfit * 1e3
    
    else:
        print(f'{cell_ID} linear fit not achieved. Number of spikes after 250 ms to low.')
        freq_adaptation_incline_linearfit = np.nan


    ## save values to dataframe
    adaptation_df.at[cell_ID, 'spike_amplitude_adaptation_ratio'] = spike_amplitude_adaptation
    adaptation_df.at[cell_ID, 'spike_FWHM_adaptation_ratio'] = spike_FWHM_adaptation
    adaptation_df.at[cell_ID, 'freq_adaptation_ratio'] = freq_adaptation_ratio
    adaptation_df.at[cell_ID, 'freq_adaptation_steadystate'] = freq_adaptation_incline_linearfit


# %% figure 

    # Exter
    colors = ['#FFEC9DFF', '#FAC881FF', '#F4A464FF', '#E87444FF', '#D9402AFF',
              '#BF2729FF', '#912534FF', '#64243EFF', '#3D1B28FF', '#161212FF']

    ax_keys = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']

    fig, axs = plt.subplot_mosaic('ABIK;CDIK;EFJL;GHJL',
                                  layout='tight',
                                  figsize=get_figure_size(),
                                  width_ratios=[3, 1, 3, 3]
                                  )

    fig.suptitle(f'{cell_ID} frequency adaptation')

    # rheobase
    axs['A'].plot(t_step, v_rheobase,
                  c=colors[0],
                  lw=1)

    axs['B'].plot(v_rheobase, dvdt_rheobase,
                  c=colors[0],
                  lw=1)

    # max freq
    axs['C'].plot(t_step, v_maxfreq,
                  c=colors[1],
                  lw=1)

    axs['D'].plot(v_maxfreq, dvdt_maxfreq,
                  c=colors[1],
                  lw=1)

    # half max freq
    axs['E'].plot(t_step, v_halfmax,
                  c=colors[2],
                  lw=1)

    axs['F'].plot(v_halfmax, dvdt_halfmax,
                  c=colors[2],
                  lw=1)

    # max inst initial freq
    axs['G'].plot(t_step, v_maxinstinitialfreq,
                  c=colors[3],
                  lw=1)

    axs['H'].plot(v_maxinstinitialfreq, dvdt_maxinstinitialfreq,
                  c=colors[3],
                  lw=1)

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

    
    # combine all spikes df to find min or max
    all_spikes = pd.concat([rheobase_spikes, maxfreq_spikes, halfmax_spikes, maxinstinitialfreq_spikes], axis = 0)
    all_ISIs = pd.concat([rheobase_ISIs, maxfreq_ISIs, halfmax_ISIs, maxinstinitialfreq_ISIs], axis = 0)

    ### IF axes ###
    IF_xmin = round_to_base(IF_df[cell_ID].dropna().index[0], 50)
    IF_xmax = round_up_to_base(IF_df[cell_ID].dropna().index[-1], 50)
    
    IF_ymax = round_up_to_base(active_properties_df.at[cell_ID, 'max_inst_initial_freq'], 20)

    axs['I'].set_title('I: Input current - frequency',
                       fontsize='small', loc='left')
    # x
    axs['I'].set_xlabel('Input current [pA]')
    axs['I'].set_xlim([IF_xmin, IF_xmax])
    axs['I'].set_xticks(np.arange(0, IF_xmax+1, 100))
    axs['I'].set_xticks(np.arange(IF_xmin, IF_xmax+1, 25), minor = True)
    
    # y
    axs['I'].set_ylabel('Firing frequency [Hz]')
    axs['I'].set_ylim([-10, IF_ymax])
    axs['I'].set_yticks(np.arange(0, IF_ymax+1, 20))
    axs['I'].set_yticks(np.arange(0, IF_ymax+1, 5), minor = True)    


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
    
    

# %% export dataframe

adaptation_df.to_excel(join(cell_descrip_dir, 'ccIF-spike_adaptation.xlsx'), index_label = 'cell_ID')    


print('Done!')


