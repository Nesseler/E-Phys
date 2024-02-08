# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:41:30 2024

@author: nesseler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
from copy import copy
from os.path import join, exists
from os import mkdir

# custom directories & parameters
from parameters.directories_win import table_dir, quant_data_dir, vplot_dir, cell_descrip_dir
from parameters.parameters import min_peak_distance, min_peak_prominence, min_spike_in_burst, bin_size_ISI_poisson, min_spikes_tobe_active
from parameters.PGFs import cc_cntrest_parameters

# custom functions
from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_useful import calc_time_series, butter_filter, calc_normed_hist, single_gaussian, double_gaussian, calc_dvdt_padded
from functions.functions_plotting import get_colors, save_figures
from functions.functions_spiketrains import calc_vmem_at_spiketrain

# %%

table = pd.read_excel(table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')

# loop to create string to include all frequencies in query
PGF = 'cc_cnt_rest'  

# limit lookup table
lookup_table = table.query(f'{PGF}.notnull()')

# cell IDs 
cell_IDs = lookup_table.query('cc_cnt_rest.notnull()').index.to_list()

# inititalize dataframe for single values for this recording
cnt_rest_values = pd.DataFrame(index = cell_IDs)

# local directories
trace_dir = join(quant_data_dir, 'cnt_rest', 'traces')
spikes_dir = join(quant_data_dir, 'cnt_rest', 'spikes')
ISIs_dir = join(quant_data_dir, 'cnt_rest', 'ISIs')
bursts_dir = join(quant_data_dir, 'cnt_rest', 'bursts')


# WARNING: significantly prolongs runtime
export_bool = False

# plotting specifications
darkmode_bool = False

vplot_bool = True

save_bool = True

colors_dict = get_colors(darkmode_bool)


# %%

for cell_ID in cell_IDs:

    vplot_dir_cell = join(vplot_dir, PGF, cell_ID)
    if not exists(vplot_dir_cell):
        mkdir(vplot_dir_cell)
    
    
    
    # get data 
        
    
    
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)
    
    # get IF data form file
    _, v, _, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')
    
    # sampling rate in ms
    SR_ms = SR / 1e3
    
    # concatenate individual steps
    n_points = int(np.shape(v)[0] * np.shape(v)[1])
    
    # reconstruct i array from parameters
    i = np.multiply(np.ones(n_points), cc_cntrest_parameters['i_hold'])
    
    v_concat = v.flatten() 

    print(f'Started: {cell_ID}')
    
    t_ms = calc_time_series(v_concat, SR)
    t_s = calc_time_series(v_concat, SR, scale = 's')
    t_total = len(t_s) / SR
    
    
    # filter voltage (to vf)
    vf = butter_filter(v_concat, 
                        order = 3,
                        cutoff = 1e3,
                        sampling_rate = SR)
    
    # calc dvdt for entire trace
    dvdt = calc_dvdt_padded(vf, t_ms)
    
    # exporting data takes time!
    if export_bool:
        # export v, vf, i, t_s, t_ms, dvdt, SR, SR_ms
        trace_df = pd.DataFrame({'v' : v_concat,
                                  'vf' : vf,
                                  'i' : i,
                                  't_ms' : t_ms,
                                  't_s' : t_s,
                                  'dvdt' : dvdt})
        
        trace_df.at[0, 'SR'] = SR
        trace_df.at[0, 'SR_ms'] = SR_ms
        
        trace_df.to_csv(join(trace_dir, f'{cell_ID}_trace_df.csv'))
    
    

    
    # %% find spikes
    idc_spikes, dict_peak = sc.signal.find_peaks(vf, 
                                                prominence = min_peak_prominence, 
                                                distance = min_peak_distance * (SR/1e3))
    
    # calculate spike times in seconds
    t_spikes = np.divide(idc_spikes, SR)
    n_spikes = len(t_spikes)
    
    # mean firing rate as number of spikes in total recorded time
    FR_mean = n_spikes / t_total
    cnt_rest_values.at[cell_ID, 'FR_mean'] = FR_mean
    
    # initialize spikes dataframe
    spikes_df = pd.DataFrame({'idc_spikes' : idc_spikes,
                              't_spikes' : t_spikes})
    
    # add to cnt_rest_values
    cnt_rest_values.at[cell_ID, 'spikes_n'] = n_spikes
    
    # categorize cell as active or not
    if n_spikes >= min_spikes_tobe_active:
        active_bool = 1
    else:
        active_bool = 0
    
    cnt_rest_values.at[cell_ID, 'active'] = active_bool
    
    # %% calc ISIs
    
    if active_bool:
        # calculate ISIs in seconds
        ISIs = np.diff(t_spikes)
        
        # calc t for ISIs as t between two spikes
        # as the t_spike + its ISI
        t_ISIs = np.add(t_spikes[:-1], np.divide(ISIs, 2))
    
    else:
        ISIs = []
        t_ISIs = []
        
    # add to spikes dataframe
    spikes_df.at[1: ,'ISIs'] = ISIs
    
    # create ISI_df
    ISI_df = pd.DataFrame({'t_ISIs' : t_ISIs, 
                            'ISIs' : ISIs})
    
    
    # %% trace and histograms
    
    
    # if active_bool:
    # calc instantanouse firing rate
    inst_FR = [1 / ISI for ISI in ISIs]
    
    # define time bins
    t_binsize = 5
    t_bins = np.arange(0, t_total + t_binsize, t_binsize)
    
    # calc binned firing rate
    binned_t_spikes = pd.DataFrame()
    
    for t_bin_idx, t_bin_u_edge in enumerate(t_bins[1:]):
        t_bin_l_edge = t_bins[t_bin_idx]
        
        bin_n_spikes = len([t_spike for t_spike in t_spikes if t_spike < t_bin_u_edge and t_spike > t_bin_l_edge])
        
        binned_t_spikes.at[t_bin_idx, 'n_spikes'] = bin_n_spikes
        binned_t_spikes.at[t_bin_idx, 'freq'] = bin_n_spikes / t_binsize
    
    
    
    if vplot_bool:
        # create figure
        fig_trace, ax_trace = plt.subplots(nrows = 3,
                                            ncols = 1,
                                            layout = 'constrained',
                                            sharex = 'col',
                                            sharey = 'row',
                                            gridspec_kw={'height_ratios': [1,1,1]})
                                            #figsize = get_figure_size())
        
        # set figure title
        fig_trace.suptitle(cell_ID)
        
        # plot voltage
        ax_trace[0].plot(t_s, vf,
                          c = colors_dict['primecolor'],
                          linewidth = 0.1)
        
        # spikes eventplots
        ax_trace[0].eventplot(spikes_df['t_spikes'],
                              orientation = 'horizontal', 
                              lineoffsets=60, 
                              linewidth = .1,
                              linelengths=10, 
                              color = [colors_dict['color3']])
            
        ax_trace[0].set_ylabel('Voltage [mV]')
        ax_trace[0].set_yticks(np.arange(-100, 70 + 1, 50))
        ax_trace[0].set_yticks(np.arange(-100, 70 + 1, 10), minor = True)
        
        
        # plot ISIs
        ax_trace[1].scatter(t_ISIs, ISIs,
                            marker = '.',
                            s = 5,
                            label = 'ISIs',
                            c = colors_dict['primecolor'])
        
        ax_trace[1].set_ylabel('ISI [s]')
        ax_trace[1].set_yscale('log', base = 10)
        ax_trace[1].set_yticks([0.001, 0.1, 10, 1000])
        
        
        # plot firing rate
        
        ax_trace[2].scatter(t_ISIs, inst_FR,
                            marker = '.',
                            s = 5,
                            color = colors_dict['primecolor'])
            
        ax_trace[2].bar(t_bins[:-1], binned_t_spikes['freq'], 
                        width = t_binsize, 
                        align = 'edge', 
                        label = 'inst. frequency',
                        facecolor = 'None',
                        edgecolor = colors_dict['color2'])
        
        ax_trace[2].set_ylabel('Firing rate [Hz]')
        
        
        # no grid
        [ax.grid(False) for ax in ax_trace]
        
        # despine panels
        [ax.spines[spine].set_visible(False) for ax in ax_trace for spine in ['top', 'right']]
        [ax.spines['bottom'].set_bounds([0, t_total]) for ax in ax_trace]
        
        ax_trace[0].spines['left'].set_bounds([-100, 70])
        
        
        ax_trace[0].set_xlim([-t_total * 0.05, t_total + t_total * 0.05])
        ax_trace[-1].set_xlabel('Time [s]')
        
        try:
            u_lim = 5 * round(np.ceil(np.nanmax(inst_FR)) / 5)
        except:
            u_lim = 20
        
        ax_trace[2].set_ylim([-0.5, u_lim])
        ax_trace[2].set_yticks(np.linspace(0, u_lim, int((u_lim / 5) + 1)))
        ax_trace[2].spines['left'].set_bounds([0, u_lim])
    
    
        # save figures
        if save_bool:
            save_figures(fig_trace, f'trace_ISIs_FRs-{cell_ID}', vplot_dir_cell, darkmode_bool)
    
    
    
    # %% ISI poisson distribution
    
    if active_bool:
        # set bins and bin size
        bin_size = bin_size_ISI_poisson # in ms
        bins = np.arange(0, (t_total / 6) + bin_size, bin_size)
        
        # calculate histogram for original ISIs
        ISI_hist, ISI_bin_edges = np.histogram(a = ISIs, bins = bins)
        
        # norm to maximal occurance
        ISI_hist_norm = calc_normed_hist(ISI_hist)
        
        # uniformly distribute the same number of spikes as measured
        # by randomly draw from uniform distribution between t0 and t_end
        uni_t_spikes = sc.stats.uniform.rvs(loc = 0, scale = t_total, size = n_spikes)
        
        # sort values to perserve time series aspect
        uni_t_spikes = np.sort(uni_t_spikes)
        
        # calculate ISIs from uniformyl distributed spikes
        ISIs_unif = np.diff(uni_t_spikes)
        
        # calculate mean ISI and mean inst. FR
        ISIs_unif_mean = np.mean(ISIs_unif)
        
        # calculate mean inst. FR as inverse of mean ISI
        FR_unif_mean = 1 / ISIs_unif_mean
        
        # histogram of uniform distributed spike times  
        uni_ISI_hist, _ = np.histogram(a = ISIs_unif, bins = bins)
        
        # norm to maximal occurance
        uni_ISI_hist_norm = calc_normed_hist(uni_ISI_hist)
        
        # inter-spike intervals of an poisson process are distributed according to an
        # exponential distribution. (Source: https://github.com/btel/python-in-neuroscience-tutorials/blob/master/poisson_process.ipynb)
        # p(x) = lambda * exp(-lambda * x)
        # p = probability
        # x = ISI (x axis)
        # lambda = firing rate
        
        # def exp_theo_pdf(x, lam):
        #     return lam * np.exp( - lam * x)
        
        # calculate the probability density function with bins given
        theo_ISIs_pdf = FR_unif_mean * np.exp( - FR_unif_mean * bins)
    
        # get maximal value
        theo_ISIs_pdf_max = np.max(theo_ISIs_pdf)
    
        # normalize pdf to maximal value
        theo_ISIs_pdf_normed = calc_normed_hist(theo_ISIs_pdf)
        
        # calculate the cumulatative density function 
        # as the inverse of the theoretical pdf
        theo_ISIs_cdf_normed = [1 - p for p in theo_ISIs_pdf_normed]
        
        # get median ISI interval (ISI_0.5)
        ISIs_theo_median = -np.log((theo_ISIs_pdf_max * 0.5) / FR_unif_mean) / FR_unif_mean 
    
        # write values to dataframe
        cnt_rest_values.at[cell_ID, 'ISIs_mean'] = np.mean(ISIs)
        cnt_rest_values.at[cell_ID, 'ISIs_median'] = np.median(ISIs)
        cnt_rest_values.at[cell_ID, 'ISIstd'] = np.std(ISIs)
        cnt_rest_values.at[cell_ID, 'ISIs_theo_median'] = ISIs_theo_median   
        
        # inst. FR quantifications
        cnt_rest_values.at[cell_ID, 'FR_unif_mean'] = FR_unif_mean
        cnt_rest_values.at[cell_ID, 'FR_inst_mean'] = np.mean(inst_FR)
        cnt_rest_values.at[cell_ID, 'FR_inst_median'] = np.median(inst_FR)
        cnt_rest_values.at[cell_ID, 'FR_inst_std'] = np.std(inst_FR)
        
        
        if vplot_bool:
            # ---------- ### spiketrain and ISIs figure ### ---------- #
            fig_hists, axs_hists = plt.subplots(nrows = 4,
                                                ncols = 1,
                                                layout = 'constrained',
                                                height_ratios = [2, 2, 2, 2])
            
            # set figure title
            fig_hists.suptitle(cell_ID)
            
            # eventplot of spikes
            axs_hists[0].eventplot(uni_t_spikes,
                                    color = colors_dict['primecolor'],
                                    linewidth = .5)
            
            axs_hists[0].eventplot(spikes_df['t_spikes'],
                                    color = colors_dict['color2'],
                                    lineoffsets = 2.5,
                                    linewidth = .5)
            
            axs_hists[0].set_xlim([-t_total * 0.05, t_total + t_total * 0.05])
            axs_hists[0].set_xlabel('Time [s]')
            
            # despine top panel
            [axs_hists[0].spines[spine].set_visible(False) for spine in ['left', 'top', 'right']]
            axs_hists[0].spines['bottom'].set_bounds([0, t_total])
            
            # ticks
            axs_hists[0].set_yticks(ticks = [1, 2.5], labels = ['Uniform', 'Original'], rotation = 30)
            
            ### histograms
            # original ISIs
            axs_hists[1].bar(ISI_bin_edges[:-1], ISI_hist_norm, 
                              width = bin_size, 
                              align = 'edge', 
                              label = 'ISIs original',
                              facecolor = 'None',
                              edgecolor = colors_dict['color2'])
            
            
            # uniform distributed spikes ISIs
            axs_hists[2].bar(ISI_bin_edges[:-1], uni_ISI_hist_norm, 
                              width = bin_size, 
                              align = 'edge', 
                              label = 'ISIs uniform spiketimes',
                              facecolor = 'None',
                              edgecolor = colors_dict['primecolor'])
            
            axs_hists[2].plot(bins, theo_ISIs_pdf_normed, 
                              c = colors_dict['primecolor'], 
                              label = 'Theoretical ISI PDF',
                              alpha = 0.5)
            
            axs_hists[2].plot(bins, theo_ISIs_cdf_normed, 
                              c = colors_dict['primecolor'], 
                              label = 'Theoretical ISI CDF',
                              alpha = 0.5,
                              linestyle = '--')
            
            axs_hists[2].vlines(ISIs_theo_median, 0, 1.0,
                                colors = 'r',
                                alpha = 0.5,
                                label = 'Theoretical median ISI')
            
            
            # normal distributed ISIs
            axs_hists[3].bar(ISI_bin_edges[:-1], ISI_hist_norm, 
                              width = bin_size, 
                              align = 'edge', 
                              label = 'ISIs original',
                              facecolor = 'None',
                              edgecolor = colors_dict['color2'])
            
            axs_hists[3].plot(bins, theo_ISIs_pdf_normed, 
                              c = colors_dict['primecolor'], 
                              label = 'Theoretical ISI PDF',
                              alpha = 0.5)
            
            axs_hists[3].plot(bins, theo_ISIs_cdf_normed, 
                              c = colors_dict['primecolor'], 
                              label = 'Theoretical ISI CDF',
                              alpha = 0.5,
                              linestyle = '--')
            
            axs_hists[3].vlines(ISIs_theo_median, 0, 1.0,
                                colors = 'r',
                                alpha = 0.5,
                                label = 'Theoretical median ISI')
            
            
            x_axis_max = np.ceil(ISIs_theo_median) * 4
            
            axs_hists[1].set_xlim([-0.2, x_axis_max + 0.2])
            axs_hists[1].set_ylim([-0.02, 1.02])
            
            [ax.set_ylabel('Normed\noccurance') for ax in axs_hists[1:]]
            axs_hists[3].set_xlabel('ISI [s]')
            
            # share x axis with second
            [ax.sharex(axs_hists[1]) for ax in axs_hists[2:]]
            [ax.sharey(axs_hists[1]) for ax in axs_hists[2:]]
            
            [ax.legend() for ax in axs_hists[1:]]
            
            # despine bottom panel
            [ax.spines[spine].set_visible(False) for ax in axs_hists[1:] for spine in ['top', 'right']]
            [ax.spines['bottom'].set_bounds([0, x_axis_max + 0.2]) for ax in axs_hists[1:]]
            [ax.spines['left'].set_bounds([0, 1]) for ax in axs_hists[1:]]
            
            # xaxis ticks
            axs_hists[3].set_xticks(np.linspace(0, x_axis_max, 5))
            
            # no grid
            [ax.grid(False) for ax in axs_hists]
            
            # save figures
            if save_bool:
                save_figures(fig_hists, f'spiketrain_n_ISI-{cell_ID}', vplot_dir_cell, darkmode_bool)
    
    
    
    # %% burst or not 
    
        # set non burst condition
        spikes_df['burst'] = [0] * len(spikes_df.index)
        spikes_df['burst_id'] = [np.nan] * len(spikes_df.index)
        
        burst_id = 0
        next_burst = False
        
        # moving window approach to burst classification
        for idx, t_spike in enumerate(t_spikes):
            
            if idx >= min_spike_in_burst:
                ISI_mov_window = ISIs[idx-4:idx+1]
                
                if all([ISI <= ISIs_unif_mean for ISI in ISI_mov_window]):
                    spikes_df.loc[idx-min_spike_in_burst:idx, 'burst'] = 1
                    spikes_df.loc[idx-min_spike_in_burst:idx, 'burst_id'] = burst_id
                    
                    next_burst = True          
                else:
                    if next_burst:
                        burst_id += 1
                        next_burst = False
        
        
        # get number of spikes within burst
        spikes_n_in_burst = spikes_df.query('burst == 1').shape[0]
        
        # calc ratio of spikes in burst to all spikes
        ratio_in_burst_to_all = spikes_n_in_burst / n_spikes
        
        # classification as burst or not
        if ratio_in_burst_to_all > 0.5:
            bursting_bool = 1
        elif ratio_in_burst_to_all <= 0.5:
            bursting_bool = 0
    
        
        # get dataframe with only spike in burst and get number of bursts
        burst_spikes_df = spikes_df.query('burst == 1')
        n_bursts = burst_spikes_df['burst_id'].max() + 1
    
        # write to dataframe
        cnt_rest_values.at[cell_ID, 'spikes_n_in_burst'] = spikes_n_in_burst
        cnt_rest_values.at[cell_ID, 'ratio_in_burst_to_all'] = ratio_in_burst_to_all
        cnt_rest_values.at[cell_ID, 'bursting'] = bursting_bool
        
        
        if vplot_bool:
            
            fig_burst, axs_burst = plt.subplots(nrows = 1,
                                                ncols = 1,
                                                layout = 'constrained')
            
            # set figure title
            if bursting_bool:
                b_str = 'bursting'
            else:
                b_str = 'not bursting'
                
            fig_burst.suptitle(f'{b_str} {cell_ID} spikes categorisation')
            
            axs_burst.eventplot(spikes_df.loc[spikes_df['burst'] == 1, 't_spikes'],
                                  orientation = 'horizontal', 
                                  lineoffsets=0, 
                                  linewidth = .1,
                                  linelengths=10, 
                                  color = [colors_dict['color3']],
                                  label = 'burst')
            
            axs_burst.eventplot(spikes_df.loc[spikes_df['burst'] == 0, 't_spikes'],
                                  orientation = 'horizontal', 
                                  lineoffsets=0, 
                                  linewidth = .1,
                                  linelengths=10, 
                                  color = [colors_dict['color2']],
                                  label = 'not burst')
    
            axs_burst.grid(False)
    
            axs_burst.set_xlim([-t_total * 0.05, t_total + t_total * 0.05])
            axs_burst.set_xlabel('Time [s]')
            
            # despine top panel
            [axs_burst.spines[spine].set_visible(False) for spine in ['left', 'top', 'right']]
            axs_burst.spines['bottom'].set_bounds([0, t_total])
            
            axs_burst.set_ylim([-20, 20])
            axs_burst.set_yticks([])
            
            axs_burst.legend()
            
            plt.show()
            
            if save_bool:
                save_figures(fig_burst, f'spike_classification-{cell_ID}', vplot_dir_cell, darkmode_bool)
    
    
    # %% burst quants
    
    if active_bool and bursting_bool:
    
        # burst DataFrame
        bursts_df = pd.DataFrame()
        
        t_all_1st_spikes = []
        t_all_lst_spikes = []
        
        burst_pad = 200e-3 # in s
    
    
        for burst_id in np.arange(n_bursts):
            
            # get spikes in burst
            burst_all_spikes = burst_spikes_df.loc[burst_spikes_df['burst_id'] == burst_id]
            burst_all_spike_idc = burst_all_spikes.index
            burst_t_all_spikes = burst_all_spikes['t_spikes'].to_list()
            burst_idc_all_spikes = [int(t * SR) for t in burst_t_all_spikes]
            
            # get first and last spike index
            burst_fst_spike_idx = burst_all_spike_idc[0]
            burst_lst_spike_idx = burst_all_spike_idc[-1]
            
            # get first and last spike time
            burst_t_1st_spike = burst_spikes_df.at[burst_fst_spike_idx, 't_spikes']
            burst_t_lst_spike = burst_spikes_df.at[burst_lst_spike_idx, 't_spikes']
            
            # add to list
            t_all_1st_spikes.append(burst_t_1st_spike)
            t_all_lst_spikes.append(burst_t_lst_spike)
        
            # get ISIs within burst (excluding first ISI)
            burst_ISIs = burst_all_spikes.loc[burst_fst_spike_idx+1:, 'ISIs'].to_list()
            burst_inst_FR = np.reciprocal(burst_ISIs)
            
            # calc mean FR
            burst_mean_FR = np.mean(burst_inst_FR)
            
            # add to burst dataframe
            bursts_df.at[burst_id, 'mean_FR'] = burst_mean_FR
        
            # get voltage trace from burst
            burst_start_idx = int((burst_t_1st_spike - burst_pad) * SR)
            burst_stop_idx = int((burst_t_lst_spike + burst_pad) * SR)
            
            # re-index spikes relative to burst
            burst_idc_all_spikes = [int(idx_spike - burst_start_idx) for idx_spike in burst_idc_all_spikes]
            
            # limit t, v, dvdt to include only burst
            burst_t = t_s[burst_start_idx:burst_stop_idx]
            burst_v = copy(vf[burst_start_idx:burst_stop_idx]) # copy to avoid parse by reference
            burst_dvdt = copy(dvdt[burst_start_idx:burst_stop_idx])
            
            # get average membrane voltage    
            bursts_df.at[burst_id, 'v_mem'] = calc_vmem_at_spiketrain(burst_t, burst_v, burst_dvdt, burst_idc_all_spikes, np.min(ISIs)*1e3*3, SR)
            
        
        burst_durations = np.subtract(t_all_lst_spikes, t_all_1st_spikes)
        
        bursts_df['burst_duration'] = burst_durations 
        
        IBIs = np.subtract(t_all_1st_spikes[1:], t_all_lst_spikes[:-1])
        
        bursts_df.at[1:, 'IBI'] = IBIs
        
        
        # write to cnt_rest_values dataframe
        cnt_rest_values.at[cell_ID, 'burst_n'] = n_bursts
        cnt_rest_values.at[cell_ID, 'burst_mean_inst_FR'] = bursts_df['mean_FR'].mean()
        cnt_rest_values.at[cell_ID, 'burst_mean_vmem'] = bursts_df['v_mem'].mean()
        cnt_rest_values.at[cell_ID, 'burst_mean_duration'] = bursts_df['burst_duration'].mean()
        cnt_rest_values.at[cell_ID, 'burst_median_duration'] = bursts_df['burst_duration'].median()
        cnt_rest_values.at[cell_ID, 'burst_mean_IBI'] = bursts_df['IBI'].mean()
        cnt_rest_values.at[cell_ID, 'burst_median_IBI'] = bursts_df['IBI'].median()
    
    # %% histogram
    
    # calc all points histogram
    
    bin_width = 0.1 # mV
    bins_min = -100
    bins_max = 40
    
    allp_bins = np.arange(bins_min, bins_max + bin_width, bin_width)
    
    allp_hist, allp_bins_edges = np.histogram(a = vf, bins = allp_bins)
    
    
    # calc normed histogram
    normed_allp_hist = calc_normed_hist(allp_hist)
    
    # calc center of bins 
    allp_bin_cens = [bin_e + bin_width / 2 for bin_e in allp_bins[:-1]]
    
    
    if active_bool and n_bursts > 3:
        # double gaussian fit
        # source blog: http://emilygraceripka.com/blog/16
        popt_dgauss, pcov_dgauss = sc.optimize.curve_fit(double_gaussian, allp_bin_cens, normed_allp_hist,
                                                          p0 = [1, -85, 3, 0.5, -70, 3])
        
        perr_2gauss = np.sqrt(np.diag(pcov_dgauss))
    
        popt_1stgauss = popt_dgauss[0:3]
        popt_2ndgauss = popt_dgauss[3:6]
        
        gauss_peak_1 = single_gaussian(allp_bin_cens, *popt_1stgauss)
        gauss_peak_2 = single_gaussian(allp_bin_cens, *popt_2ndgauss)
        gauss_both_peaks = double_gaussian(allp_bin_cens, *popt_dgauss)
        
        # save fitted center to v_up andn v_down states in DataFrame for cell
        v_up = popt_2ndgauss[1]
        v_down = popt_1stgauss[1]
        cnt_rest_values.at[cell_ID, 'v_up_allp_fgauss'] = v_up
        cnt_rest_values.at[cell_ID, 'v_down_allp_fgauss'] = v_down
    
    else:
        # single gaussian fit
        popt_sgauss, pcov_sgauss = sc.optimize.curve_fit(single_gaussian, allp_bin_cens, normed_allp_hist,
                                                          p0 = [1, -85, 3])
    
        gauss_speak = single_gaussian(allp_bin_cens, *popt_sgauss)
        v_rest_allp_gauss = popt_sgauss[1]
        
        # write to dataframe
        cnt_rest_values.at[cell_ID, 'v_rest_allp_gauss '] = v_rest_allp_gauss 
        

    
    if vplot_bool:
        # ---------- ### all points histograms (abs. and rel.) ### ---------- #
        
        fig_hist, ax_hist = plt.subplots(nrows = 1,
                                          ncols = 1,
                                          layout = 'constrained',
                                          sharex = 'col',
                                          sharey = 'row')
        
        fig_hist.suptitle(f'absolute histogram {cell_ID}')
        
        ax_hist.bar(allp_bins[:-1], allp_hist,
                    width = bin_width,
                    align = 'edge',
                    facecolor = 'None',
                    edgecolor = colors_dict['color1'])
        
        
        ax_hist.grid(False)
        
        ax_hist.set_ylim([-(n_points * 0.01), n_points])
        ax_hist.set_yticks(np.arange(0, n_points + 1, 5e6))
        ax_hist.ticklabel_format(axis = 'y', scilimits = (6, 6), useMathText = True)
        ax_hist.set_ylabel('Number of points')
        
        # despine panels
        [ax_hist.spines[spine].set_visible(False) for spine in ['top', 'right']]
        ax_hist.spines['left'].set_bounds([0, n_points])
        ax_hist.spines['bottom'].set_bounds([bins_min, bins_max])
        
        
        plt.show()
        
        if save_bool:
            save_figures(fig_hist, f'allp_histo_abs-{cell_ID}', vplot_dir_cell, darkmode_bool)
        
        
        ### normed histogram ###
        
        fig_relhist, ax_relhist = plt.subplots(nrows = 1,
                                                ncols = 1,
                                                layout = 'constrained')
        
        fig_relhist.suptitle(f'normed histogram {cell_ID}')
        
        
        
        ax_relhist.scatter(allp_bin_cens, normed_allp_hist,
                            color = colors_dict['primecolor'],
                            marker = '.',
                            s = 10,
                            label = 'all points histogram')
        
        if active_bool and n_bursts > 3:
            
            ax_relhist.plot(allp_bin_cens, gauss_both_peaks,
                            label = 'double gauss fit',
                            color = 'gray', 
                            linewidth = 1)
        
            ax_relhist.plot(allp_bin_cens, gauss_peak_1,
                            label = 'v_up gauss fit',
                            color = colors_dict['color2'], 
                            linewidth = 1,
                            linestyle = '--')
            
            ax_relhist.plot(allp_bin_cens, gauss_peak_2,
                            label = 'v_down gauss fit',
                            color = colors_dict['color3'], 
                            linewidth = 1,
                            linestyle = '--')
            
            ax_relhist.arrow(x = v_up,
                              y = 0.02,
                              dx = 0,
                              dy = -0.03,
                              length_includes_head = True,
                              head_width = 1,
                              head_length = 0.03,
                              color = colors_dict['color2'],
                              linewidth = 1,
                              label = '_nolegend_')
            
            ax_relhist.arrow(x = v_down,
                              y = 0.02,
                              dx = 0,
                              dy = -0.03,
                              length_includes_head = True,
                              head_width = 1,
                              head_length = 0.03,
                              color = colors_dict['color3'],
                              linewidth = 1,
                              label = '_nolegend_')
        
        else:
            ax_relhist.plot(allp_bin_cens, gauss_speak,
                            label = 'v_rest gauss fit',
                            color = colors_dict['color3'], 
                            linewidth = 1,
                            linestyle = '--')
            
            ax_relhist.arrow(x = v_rest_allp_gauss,
                              y = 0.02,
                              dx = 0,
                              dy = -0.03,
                              length_includes_head = True,
                              head_width = 1,
                              head_length = 0.03,
                              color = colors_dict['color3'],
                              linewidth = 1,
                              label = '_nolegend_')
            
        
        ax_relhist.grid(False)
        ax_relhist.legend()
        
        ax_relhist.set_ylim([-0.01, 1])
        ax_relhist.set_yticks(np.arange(0, 1 + 0.1, 0.2))
        ax_relhist.set_ylabel('Normed number of points')
        
        # despine panels
        [ax_relhist.spines[spine].set_visible(False) for spine in ['top', 'right']]
        ax_relhist.spines['left'].set_bounds([0, 1])
        ax_relhist.set_xlim([-102, -40])
        ax_relhist.spines['bottom'].set_bounds([-100, -40])
        
        
        if save_bool:
            save_figures(fig_relhist, f'allp_histo_normed-{cell_ID}', vplot_dir_cell, darkmode_bool)
        
        plt.show()
    
    



    # %% save dataframes

    spikes_df.to_excel(join(spikes_dir, f'cccntrest-spikes-{cell_ID}.xlsx'))

    if active_bool:        
        ISI_df.to_excel(join(ISIs_dir, f'cccntrest-ISIs-{cell_ID}.xlsx'))   
        
    if active_bool and bursting_bool:
        bursts_df.to_excel(join(bursts_dir, f'cccntrest-bursts-{cell_ID}.xlsx'))




cnt_rest_values.to_excel(join(cell_descrip_dir, 'cccntrest-values.xlsx'))








