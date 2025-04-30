# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 13:45:39 2025

@author: nesseler
"""


from scipy.ndimage import maximum_filter1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sc
import sys
# sys.path.append('../../core/')
from miniML import MiniTrace, EventDetection
from miniML_plot_functions import miniML_plots

%matplotlib inline

import matplotlib as mtl
# set font size
mtl.rcParams.update({'font.size': 9, 'font.family' : 'Arial'})


filename = 'C:/Users/nesseler/miniML/example_data/gc_mini_trace.h5'
scaling = 1e12
unit = 'pA'

# # load from h5 file
# trace = MiniTrace.from_h5_file(filename=filename,
#                                tracename='mini_data',
#                                scaling=scaling,
#                                sampling=2e-5,
#                                unit=unit)

filename = 'Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-021-blkr.dat'
rectype = 'vc-Erest-3min'
scaling = 1e12
unit = 'pA'
SR = 100_000

trace = MiniTrace.from_heka_file(filename=filename,
                                 rectype=rectype,
                                 group=1,
                                 exclude_series = np.delete(np.arange(0, 26), 13),
                                 scaling=scaling,
                                 unit=unit)

b, a = sc.signal.bessel(3, 750, fs=100e3)

ori_trace = trace.data

data_filtered = sc.signal.lfilter(b, a, ori_trace)


# %%
trace.plot_trace()

trace = MiniTrace(data=data_filtered,
                  sampling_interval=1/100e3,
                  y_unit='pA',
                  filename='None')

trace.plot_trace()

# %%

# factors = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 35, 40, 45, 50]

# # factors = [5, 19, 30]
# thresholds = [0.25, 0.5]

# detection_stats = pd.DataFrame(columns = ['n_events', 'u_score'],
#                                index = [str(f) + '_' + str(t) for t in thresholds for f in factors])

# predictions = pd.DataFrame(columns = [str(f) + '_' + str(t) for t in thresholds for f in factors],
#                            index = np.arange(0, data_filtered.shape[0]))

# event_peakidc = []

# for factor in factors:

#     for th in thresholds:    
        
#         # win_size = 600 * factor
#         direction = 'negative'
        
#         eventdetection_settings = {'window_size' : 600 * factor,
#                                    'model_threshold' : th,
#                                    'batch_size' : 512,
#                                    'event_detection_peakw' : 5,
#                                    'stride' : 30,
#                                    'rel_prom_cutoff' : 0.25,
#                                    'convolve_win' : 20 * factor}
    
#         detection = EventDetection(data=trace,
#                                    model_path='C:/Users/nesseler/miniML/models/GC_lstm_model.h5',
#                                    window_size = eventdetection_settings['window_size'],
#                                    model_threshold = eventdetection_settings['model_threshold'],
#                                    batch_size = eventdetection_settings['batch_size'],
#                                    event_direction=direction,
#                                    compile_model=True,
#                                    verbose=2)
    
        
    
#         detection.detect_events(eval=True,
#                                 stride = eventdetection_settings['stride'],
#                                 peak_w = eventdetection_settings['event_detection_peakw'],
#                                 rel_prom_cutoff = eventdetection_settings['rel_prom_cutoff'],
#                                 convolve_win = eventdetection_settings['convolve_win'],
#                                 # gradient_convolve_win = 40 * factor,
#                                 resample_to_600 = True)
        
#         # write detection object to dict
#         detection_stats.at[str(factor) + '_' + str(th), 'n_events'] = detection.event_stats.event_count
#         detection_stats.at[str(factor) + '_' + str(th), 'u_score'] = np.mean(detection.event_scores)
        
#         # prediction
#         filtered_prediction = maximum_filter1d(detection.prediction, size=int(5*detection.interpol_factor), origin=-2)
#         predictions[str(factor) + '_' + str(th)] = np.append(filtered_prediction, [np.nan]*(data_filtered.shape[0]-filtered_prediction.shape[0]))
        
#         # events (peak locations)
#         event_peakidc.append(detection.event_peak_locations)
    
#         # miniML plots
#         MiniPlots = miniML_plots(data=detection)
#         MiniPlots.plot_event_overlay()
#         MiniPlots.plot_event_histogram(plot='amplitude', cumulative=False)
    
    
    
# %%

fig, axs = plt.subplots(nrows = 2, ncols = 1,
                        sharex = True,
                        dpi = 300)

for th, c in zip(thresholds, ['gray', 'k']):
    rows = [str(f) + '_' + str(th) for f in factors]

    axs[0].plot(factors, detection_stats.loc[rows, 'n_events'],
                lw = 1,
                color = c,
                marker = '.',
                ms = 3,
                label = f'thres: {th}')
    
    axs[1].plot(factors, detection_stats.loc[rows, 'u_score'],
                lw = 1,
                color = c,
                marker = '.',
                ms = 3,
                label = f'thres: {th}')

axs[0].legend(frameon = False,
              fontsize = 7,
              loc = 'upper right')

axs[0].set_ylim([0, 550])
axs[0].set_ylabel('Number of\ndetected events [#]')
    
axs[1].legend(frameon = False,
              fontsize = 7,
              loc = 'lower right')

axs[1].set_ylim([0.70, 1])
axs[1].set_ylabel('Average\nevent score')
axs[1].set_xlim([0, 51])
axs[1].set_xlabel('Sliding window size [ms]')
axs[1].set_xticks(ticks = np.arange(0, 50+0.1, 10), 
                  labels = np.arange(0, (600/(SR/1e3))*(50+0.1), (600/(SR/1e3))*10, dtype = int))
axs[1].set_xticks(ticks = np.arange(0, 51+0.1, 1),
                  minor = True)

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

fig.align_labels()

plt.show()


# %%

factor = 19
factor_i = factors.index(factor)
th = 0.25
th_i = thresholds.index(th)
event_i = (factor_i*2)+th_i

parent_fig = plt.figure(layout='constrained',
                        figsize = (6, 3),
                        dpi=300)

parent_fig.suptitle(f'Trace\nwindow: {(600*factor)/(SR/1e3)} ms - thres: {th}',
                    fontsize = 9)


x_full = np.arange(0, data_filtered.shape[0] / SR, step=1/SR)
# x_predic = np.arange(0, predictions.loc[:, factor].dropna().shape[0] / SR, step=1/SR)

axs = parent_fig.subplots(nrows=2,
                          ncols=1,
                          sharex=True,
                          height_ratios=[1, 3])

# plot prediction
axs[0].plot(x_full, predictions.loc[:, str(factor) + '_' + str(th)],
            color='k',
            alpha=0.5,
            lw=0.5,
            label='filtered')

axs[0].hlines(xmin = 0, xmax = x_full[-1],
              y = th, 
              color = 'r',
              lw = 0.75,
              ls = 'dashed')

# plot data
axs[1].plot(x_full, ori_trace,
            color='k',
            alpha = 0.25,
            lw=0.5,
            label='data')

# plot data
axs[1].plot(x_full, data_filtered,
            color='k',
            lw=0.5,
            label='data (filt)')

# plot event peak indicators
axs[1].scatter(x = x_full[event_peakidc[event_i]],
               y = data_filtered[event_peakidc[event_i]],
               marker = 'o',
               color = 'r',
               s = 5,
               lw = 1)

# plot eventplot
axs[1].eventplot(positions = event_peakidc[event_i] / SR,
                 orientation = 'horizontal',
                 lineoffsets = -18,
                 linelengths = 1,
                 linewidths = 1,
                 color = 'r',
                 label = 'events')

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# axis 
axs[0].set_ylim([-0.05, 1.05])
axs[0].set_ylabel('Probability')

axs[1].set_xlim([1, 2])
axs[1].set_xlabel('Time [s]')

axs[1].set_ylim([-20, 10])
axs[1].set_ylabel('Current [pA]')

parent_fig.align_labels()

plt.show()


# %%
factors_old = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 35, 40, 45, 50]
factors = [5, 19, 30]
thresholds = [0.25, 0.5]

parent_fig = plt.figure(layout='constrained',
                        figsize = (6, 3),
                        dpi=300)

parent_fig.suptitle(f'Trace - Prediction - Events',
                    fontsize = 9)


x_full = np.arange(0, data_filtered.shape[0] / SR, step=1/SR)
# x_predic = np.arange(0, predictions.loc[:, factor].dropna().shape[0] / SR, step=1/SR)

plt_idc = np.arange(1 * SR, 2 * SR, step = 1, dtype = int)


axs = parent_fig.subplots(nrows=3,
                          ncols=1,
                          sharex=True,
                          height_ratios=[1, 2, 1])

# plot data
axs[0].plot(x_full[plt_idc], data_filtered[plt_idc],
            color='k',
            lw=0.5,
            label='data (filt)')

for fi, factor in enumerate(factors):
    
    axs[1].plot(x_full[plt_idc], predictions.loc[plt_idc, str(factor) + '_' + str(0.5)] - fi,
                color='k',
                alpha=0.5,
                lw=0.5,
                label='filtered')
    
    axs[1].text(x = x_full[plt_idc][1000],
                y = -fi + 0.1,
                s = f'w: {(600*factor)/(SR/1e3)} ms',
                fontsize = 6,
                color = 'gray')
    
    for ti, th in enumerate(thresholds): 
        
        axs[1].hlines(xmin = x_full[plt_idc][0], xmax = x_full[plt_idc][-1],
                      y = - fi+th, 
                      color = 'r',
                      lw = 0.25,
                      ls = 'dashed')
        
        factor_i = factors_old.index(factor)
        th_i = thresholds.index(th)

        axs[2].eventplot(positions = event_peakidc[(factor_i*2)+th_i] / SR,
                         orientation = 'horizontal',
                         lineoffsets = -(fi*2)+ti,
                         linelengths = 0.8,
                         linewidths = 1,
                         color = 'r',
                         label = 'events')
        
        axs[2].text(x = x_full[plt_idc][10000],
                    y = -(fi*2)+ti -0.5,
                    s = f'th: {th}',
                    fontsize = 6,
                    color = 'r')
        
    axs[2].text(x = x_full[plt_idc][1000],
                y = -(fi*2)+ti -1,
                s = f'w: {(600*factor)/(SR/1e3)} ms',
                fontsize = 6,
                color = 'gray')
    

axs[1].set_ylim([-2.05, 1.05])

axs[2].set_xlim([1, 2])
axs[2].set_xlabel('Time [s]')

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

parent_fig.align_labels()

plt.show()