# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:41:30 2024

@author: nesseler
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

# custom directories & parameters
from parameters.directories_win import table_dir
from parameters.parameters import min_peak_distance, min_peak_prominence, min_spike_in_burst

# custom functions
from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_useful import calc_time_series, butter_filter
from functions.functions_plotting import get_colors, get_figure_size


table = pd.read_excel(table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')


# loop to create string to include all frequencies in query
PGF = 'cc_cnt_rest'  

# limit lookup table
lookup_table = table.query(f'{PGF}.notnull()')

# cell IDs 
cell_ID = 'E-113'


# %%


    
print(f'Started: {cell_ID}')

# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)

# get IF data form file
i, v, t, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')

# sampling rate in ms
SR_ms = SR / 1e3

# concatenate individual steps
n_points = int(np.shape(i)[0] * np.shape(i)[1])

i_concat = i.flatten() 

v_concat = v.flatten() 

t_ms = calc_time_series(v_concat, SR)
t_s = calc_time_series(v_concat, SR, scale = 's')
t_total = len(t_s) / SR

# plt.plot(v_concat[0:15000])
# plt.show()


# filter voltage (to vf)
vf = butter_filter(v_concat, 
                   order = 3,
                   cutoff = 1e3,
                   sampling_rate = SR)



# find peaks
idx_peaks, dict_peak = sc.signal.find_peaks(vf, 
                                            prominence = min_peak_prominence, 
                                            distance = min_peak_distance * (SR/1e3))

t_peaks_s = np.divide(idx_peaks, SR)

spikes_df = pd.DataFrame({'t_spikes' : t_peaks_s})


# %% calc ISIs

# calculate ISIs in seconds
ISIs_s = np.diff(t_peaks_s)

# pad array with nan for first spike
ISIs_s = np.pad(ISIs_s,
                pad_width = (1, 0),
                mode = 'constant',
                constant_values = (np.nan,))

# write to spikes_df
spikes_df['ISIs'] = ISIs_s


# %% poisson distribution

# set bins and bin size
bin_size = 25e-3 # in ms
bins = np.arange(0, 20 + bin_size, bin_size)

# calculate histogram for original ISIs
ISI_hist, ISI_bin_edges = np.histogram(a = ISIs_s, bins = bins)
n_ISIs = np.sum(ISI_hist)
max_n_ISIs = np.max(ISI_hist)

# norm to maximal occurance
ISI_hist_norm = [n / max_n_ISIs for n in ISI_hist]



colors_dict = get_colors(False)

n_spikes = len(t_peaks_s)
n_ISIs = n_spikes -1

uni_t_spikes = sc.stats.uniform.rvs(loc = 0, scale = t_total, size = n_spikes)
uni_t_spikes = np.sort(uni_t_spikes)
uni_ISIs = np.diff(uni_t_spikes)
uni_mean_ISI = np.mean(uni_ISIs)
uni_mean_FR = 1 / uni_mean_ISI
uni_std_ISI = np.std(uni_ISIs)




# histogram of uniform distributed spike times

uni_ISI_hist, uni_ISI_bin_edges = np.histogram(a = uni_ISIs, bins = bins)
uni_n_ISIs = np.sum(uni_ISI_hist)
uni_max_n_ISIs = np.max(uni_ISI_hist)
uni_ISI_hist_norm = [n / uni_max_n_ISIs for n in uni_ISI_hist]

uni_ISI_hist_cumsum = np.cumsum(uni_ISI_hist)
uni_ISI_hist_cumsum_normed = [n / n_ISIs for n in uni_ISI_hist_cumsum]



# inter-spike intervals of an poisson process are distributed according to an
# exponential distribution. (Source: https://github.com/btel/python-in-neuroscience-tutorials/blob/master/poisson_process.ipynb)
# p(x) = lambda * exp(-lambda * x)
# p = probability
# x = ISI (x axis)
# lambda = firing rate

# def exp_theo_pdf(x, lam):
#     return lam * np.exp( - lam * x)


theo_ISIs_pdf = uni_mean_FR * np.exp( - uni_mean_FR * bins)
theo_ISIs_pdf_max = np.max(theo_ISIs_pdf)

theo_ISIs_pdf_normed = [o / theo_ISIs_pdf_max for o in theo_ISIs_pdf]

# calculate the cumulatative density function ### QUESTION! ###
# theo_ISIs_cdf = np.cumsum(theo_ISIs_pdf)
# theo_ISIs_cdf_max = np.max(theo_ISIs_cdf)

# theo_ISIs_cdf_normed = [o / theo_ISIs_cdf_max for o in theo_ISIs_cdf]
theo_ISIs_cdf_normed = [1 - p for p in theo_ISIs_pdf_normed]



# get median ISI interval (ISI_0.5)
theo_median_ISI = -np.log((theo_ISIs_pdf_max * 0.5) / uni_mean_FR) / uni_mean_FR



# histogram figure

fig_hists, axs_hists = plt.subplots(nrows = 4,
                                    ncols = 1,
                                    layout = 'constrained',
                                    height_ratios = [1, 2, 2, 2])

# eventplot of spikes
axs_hists[0].eventplot(uni_t_spikes,
                       color = colors_dict['primecolor'],
                       linewidth = .5)

axs_hists[0].eventplot(spikes_df['t_spikes'],
                       color = colors_dict['color2'],
                       lineoffsets = 2.5,
                       linewidth = .5)

axs_hists[0].set_xlim([-10, t_total + 10])
axs_hists[0].set_xlabel('Time [s]')

# despine top panel
[axs_hists[0].spines[spine].set_visible(False) for spine in ['left', 'top', 'right']]
axs_hists[0].spines['bottom'].set_bounds([0, t_total])


# no ticks
axs_hists[0].set_yticks(ticks = [1, 2.5], labels = ['Uniform', 'Original'], rotation = 30)


### histogram
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

axs_hists[2].vlines(theo_median_ISI, 0, 1.0,
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

axs_hists[3].vlines(theo_median_ISI, 0, 1.0,
                    colors = 'r',
                    alpha = 0.5,
                    label = 'Theoretical median ISI')


x_axis_max = np.ceil(theo_median_ISI) * 4

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




# %% burst or not 

# set non burst condition
spikes_df['burst'] = [0] * len(spikes_df.index)
spikes_df['burst_id'] = [np.nan] * len(spikes_df.index)


burst_id = 0
next_burst = False



# moving window approach to burst classification
for idx, t_spike in enumerate(t_peaks_s):
    
    if idx >= min_spike_in_burst:
        ISI_mov_window = ISIs_s[idx-3:idx+1]
        
        if all([ISI <= uni_mean_ISI for ISI in ISI_mov_window]):
            spikes_df.loc[idx-min_spike_in_burst:idx, 'burst'] = 1
            spikes_df.loc[idx-min_spike_in_burst:idx, 'burst_id'] = burst_id
            
            next_burst = True          
        else:
            if next_burst:
                burst_id += 1
                next_burst = False



# %%

t_binsize = 1
t_bins = np.arange(0, t_total + t_binsize, t_binsize)

fig_trace, ax_trace = plt.subplots(nrows = 3,
                                   ncols = 1,
                                   layout = 'constrained',
                                   sharex = 'col',
                                   sharey = 'row',
                                   gridspec_kw={'height_ratios': [1,1,1]})



# plot voltage
ax_trace[0].plot(t_s, vf)

# plot spikes in burst
ax_trace[0].eventplot(spikes_df.loc[spikes_df['burst'] == 1,'t_spikes'],
                      orientation = 'horizontal', 
                      lineoffsets=60, 
                      linewidth = .5,
                      linelengths=10, 
                      color = [colors_dict['color2']])

# plot spikes not in burst
ax_trace[0].eventplot(spikes_df.loc[spikes_df['burst'] == 0,'t_spikes'],
                      orientation = 'horizontal', 
                      lineoffsets=60, 
                      linewidth = .5,
                      linelengths=10, 
                      color = [colors_dict['color1']])

ax_trace[0].set_ylabel('Voltage [mV]')
ax_trace[0].set_yticks(np.arange(-100, 70 + 1, 50))
ax_trace[0].set_yticks(np.arange(-100, 70 + 1, 10), minor = True)


# plot ISIs
ax_trace[1].plot(t_peaks_s, ISIs_s, '-o')

ax_trace[1].set_ylabel('ISI [s]')



ax_trace[2].plot(t_peaks_s, [1 / ISI for ISI in ISIs_s])



# no grid
[ax.grid(False) for ax in ax_trace]

# despine panels
[ax.spines[spine].set_visible(False) for ax in ax_trace for spine in ['top', 'right']]
[ax.spines['bottom'].set_bounds([0, t_total]) for ax in ax_trace]

ax_trace[0].spines['left'].set_bounds([-100, 70])


ax_trace[0].set_xlim([-10, t_total + 10])


# %% histogram


fig_hist, ax_hist = plt.subplots(nrows = 2,
                                 ncols = 2,
                                 layout = 'constrained',
                                 sharex = 'col',
                                 sharey = 'row',
                                 gridspec_kw={'height_ratios': [1, 3], 'width_ratios': [2, 1]},
                                 figsize = get_figure_size())

fig_hist.suptitle(cell_ID)

bin_width = .5 # mV
bins_min = -100
bins_max = 40

ylim_border = 0.1

bins = np.arange(bins_min, bins_max + bin_width, bin_width)


N, bins, patches = ax_hist[1][0].hist(vf, bins = bins, density = True)


ax_hist[1][0].set_ylim([0, ylim_border])
ax_hist[1][0].set_xticks(np.arange(bins_min, bins_max + bin_width, 20))
ax_hist[1][0].set_xticks(bins, minor = True)



N, bins, patches = ax_hist[0][0].hist(vf, bins = bins, density = True)


ax_hist[0][0].set_ylim([ylim_border, 1])


N, bins, patches = ax_hist[0][1].hist(vf, bins = bins, density = True)
N, bins, patches = ax_hist[1][1].hist(vf, bins = bins, density = True)

ax_hist[1][1].set_xticks(np.arange(bins_min, bins_max + bin_width, 20))
ax_hist[1][1].set_xticks(bins, minor = True)
ax_hist[1][1].set_xlim([-90, -60])

[ax.grid(False) for axs in ax_hist for ax in axs]








