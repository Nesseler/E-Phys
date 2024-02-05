# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:41:30 2024

@author: nesseler
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn
import numpy as np
import scipy as sc
from copy import copy

# custom directories & parameters
from parameters.directories_win import table_dir
from parameters.parameters import min_peak_distance, min_peak_prominence, min_spike_in_burst, dvdt_threshold

# custom functions
from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_useful import calc_time_series, butter_filter, calc_normed_hist, single_gaussian, double_gaussian, calc_dvdt_padded
from functions.functions_plotting import get_colors, get_figure_size
from functions.functions_spiketrains import extract_spike, calc_vmem_at_spiketrain

table = pd.read_excel(table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')


# loop to create string to include all frequencies in query
PGF = 'cc_cnt_rest'  

# limit lookup table
lookup_table = table.query(f'{PGF}.notnull()')

# cell IDs 
cell_ID = 'E-113'

cnt_rest_parameters = pd.DataFrame(index = [cell_ID])


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

# calc t for ISIs as t between two spikes
# as the t_spike + its ISI
t_ISIs = np.add(t_peaks_s[:-1], np.divide(ISIs_s, 2))

# pad array with nan for first spike
padISIs_s = np.pad(ISIs_s,
                   pad_width = (1, 0),
                   mode = 'constant',
                   constant_values = (np.nan,))

spikes_df['ISIs'] = padISIs_s

# write to ISI_df
ISI_df = pd.DataFrame({'t_ISIs' : t_ISIs, 
                       'ISIs' : ISIs_s})


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

cnt_rest_parameters.at[cell_ID, 'theo_median_ISI'] = theo_median_ISI

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
        ISI_mov_window = ISIs_s[idx-4:idx+1]
        
        if all([ISI <= uni_mean_ISI for ISI in ISI_mov_window]):
            spikes_df.loc[idx-min_spike_in_burst:idx, 'burst'] = 1
            spikes_df.loc[idx-min_spike_in_burst:idx, 'burst_id'] = burst_id
            
            next_burst = True          
        else:
            if next_burst:
                burst_id += 1
                next_burst = False



n_spikes_in_burst = spikes_df.query('burst == 1').shape[0]

ratio_in_burst_to_all = n_spikes_in_burst / n_spikes

if ratio_in_burst_to_all > 0.5:
    cnt_rest_parameters.at[cell_ID, 'bursting'] = 1
elif ratio_in_burst_to_all <= 0.5:
    cnt_rest_parameters.at[cell_ID, 'bursting'] = 0


# %%

# calc instantanouse firing rate
inst_FR = [1 / ISI for ISI in ISIs_s]

# define time bins
t_binsize = 5
t_bins = np.arange(0, t_total + t_binsize, t_binsize)

# calc binned firing rate
binned_t_spikes = pd.DataFrame()

for t_bin_idx, t_bin_u_edge in enumerate(t_bins[1:]):
    t_bin_l_edge = t_bins[t_bin_idx]
    
    bin_n_spikes = len([t_spike for t_spike in t_peaks_s if t_spike < t_bin_u_edge and t_spike > t_bin_l_edge])
    
    binned_t_spikes.at[t_bin_idx, 'n_spikes'] = bin_n_spikes
    binned_t_spikes.at[t_bin_idx, 'freq'] = bin_n_spikes / t_binsize




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
ax_trace[1].scatter(t_ISIs, ISIs_s,
                    marker = '.',
                    s = 5,
                    label = 'ISIs')

ax_trace[1].set_ylabel('ISI [s]')
ax_trace[1].set_yscale('log', base = 10)
ax_trace[1].set_yticks([0.001, 0.1, 10, 1000])


# plot firing rate

ax_trace[2].scatter(t_ISIs, inst_FR,
                    marker = '.',
                    s = 5)
    
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


ax_trace[0].set_xlim([-10, t_total + 10])
ax_trace[-1].set_xlabel('Time [s]')

u_lim = 5 * round(np.ceil(np.nanmax(inst_FR)) / 5)

ax_trace[2].set_ylim([-0.5, u_lim])
ax_trace[2].set_yticks(np.linspace(0, u_lim, int((u_lim / 5) + 1)))
ax_trace[2].spines['left'].set_bounds([0, u_lim])







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
v_up = popt_1stgauss[1]
v_down = popt_2ndgauss[1]
cnt_rest_parameters.at[cell_ID, 'v_up_fgauss'] = v_up
cnt_rest_parameters.at[cell_ID, 'v_down_fgauss'] = v_down







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

ax_hist.set_ylim([-0.3e6, 30e6])
ax_hist.set_yticks(np.arange(0, 30e6 + 1, 5e6))
ax_hist.ticklabel_format(axis = 'y', scilimits = (6, 6), useMathText = True)
ax_hist.set_ylabel('Number of points')

# despine panels
[ax_hist.spines[spine].set_visible(False) for spine in ['top', 'right']]
ax_hist.spines['left'].set_bounds([0, 30e6])
ax_hist.spines['bottom'].set_bounds([bins_min, bins_max])


plt.show()




fig_relhist, ax_relhist = plt.subplots(nrows = 1,
                                       ncols = 1,
                                       layout = 'constrained')

fig_relhist.suptitle(f'normed histogram {cell_ID}')



ax_relhist.scatter(allp_bin_cens, normed_allp_hist,
                   color = colors_dict['primecolor'],
                   marker = '.',
                   s = 10,
                   label = 'all points histogram')

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



plt.show()



# %% burst parameters

plt.plot(t_s, vf)


burst_spikes_df = spikes_df.query('burst == 1')


# plot spikes in burst
plt.eventplot(spikes_df.loc[spikes_df['burst'] == 1,'t_spikes'],
                      orientation = 'horizontal', 
                      lineoffsets=40, 
                      linewidth = .5,
                      linelengths=10, 
                      color = [colors_dict['color2']])

# plot spikes not in burst
plt.eventplot(spikes_df.loc[spikes_df['burst'] == 0,'t_spikes'],
                      orientation = 'horizontal', 
                      lineoffsets=40, 
                      linewidth = .5,
                      linelengths=10, 
                      color = [colors_dict['color1']])


# plt.xlim([350, 380])

# plt.xlim([550, 600])

# plt.xlim([300, 350])

plt.xlim([0, 30])


# %% burst quants





# burst DataFrame
bursts_df = pd.DataFrame()

n_bursts = burst_spikes_df['burst_id'].max() + 1

t_all_1st_spikes = []
t_all_lst_spikes = []

burst_pad = 200e-3 # in s


# calc dvdt for entire trace
dvdt = calc_dvdt_padded(vf, t_ms)



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
    bursts_df.at[burst_id, 'v_mem'] = calc_vmem_at_spiketrain(burst_v, burst_dvdt, burst_idc_all_spikes, SR, dvdt_threshold)
    

burst_durations = np.subtract(t_all_lst_spikes, t_all_1st_spikes)

bursts_df['burst_duration'] = burst_durations 

IBIs = np.subtract(t_all_1st_spikes[1:], t_all_lst_spikes[:-1])

bursts_df.at[1:, 'IBI'] = IBIs 




# %% autcorrelation

def crosscorr(datax, datay, lag=0):
    """ Lag-N cross correlation.
    (Kindly provided by Dr. Christopher Wiesbrock.)
    Parameters:
        lag : int, default 0
        datax, datay : pandas.Series objects of equal length
    Returns:
        crosscorr : float
    """
    return datax.corr(datay.shift(lag))


vf_ds = pd.Series(vf)

lag_s = 600
lag_points = lag_s * SR

lag_steps_s = 10000e-3
lag_steps_points = lag_steps_s * SR

lag = np.arange(0, lag_points, lag_steps_points)
autocorr_values=np.zeros((len(lag)))

for i in range(len(autocorr_values)):
    autocorr_values[i]=crosscorr(vf_ds, vf_ds, lag=i)
    print(i / len(lag))


plt.scatter(lag, autocorr_values)


