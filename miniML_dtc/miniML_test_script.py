# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 16:00:56 2025

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


# %matplotlib inline

import matplotlib as mtl
# set font size
mtl.rcParams.update({'font.size': 9, 'font.family' : 'Arial'})




scaling = 1e12
unit = 'pA'
SR = 100_000


# rectype = 'vc-Ek-3min'

# E-238
# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-007-TTX.dat',
#                                   rectype=rectype,
#                                   group=2,
#                                   exclude_series = np.delete(np.arange(0, 100), 12-1),
#                                   scaling=scaling,
#                                   unit=unit)

# E-300
# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-021-blkr.dat',
#                                   rectype=rectype,
#                                   group=3,
#                                   exclude_series = np.delete(np.arange(0, 100), 15-1),
#                                   scaling=scaling,
#                                   unit=unit)


rectype = 'vc-Erest-3min'

# E-297
trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-020-TTX_Cd.dat',
                                  rectype=rectype,
                                  group=3,
                                  exclude_series = np.delete(np.arange(0, 100), 22-1),
                                  scaling=scaling,
                                  unit=unit)

# E-277
# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-016-blkr.dat',
#                                   rectype=rectype,
#                                   group=2,
#                                   exclude_series = np.delete(np.arange(0, 100), 17-1),
#                                   scaling=scaling,
#                                   unit=unit)

# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-021-blkr.dat',
#                                  rectype=rectype,
#                                  group=4,
#                                  exclude_series = np.delete(np.arange(0, 24), 9),
#                                  scaling=scaling,
#                                  unit=unit)

# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-021-blkr.dat',
#                                  rectype=rectype,
#                                  group=1,
#                                  exclude_series = np.delete(np.arange(0, 26), 13),
#                                  scaling=scaling,
#                                  unit=unit)


# E-300
# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-021-blkr.dat',
#                                   rectype=rectype,
#                                   group=3,
#                                   exclude_series = np.delete(np.arange(0, 100), 13-1),
#                                   scaling=scaling,
#                                   unit=unit)




# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-022-blkr.dat',
#                                  rectype=rectype,
#                                  group=1,
#                                  exclude_series = np.delete(np.arange(0, 11), 8),
#                                  scaling=scaling,
#                                  unit=unit)

# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-022-blkr.dat',
#                                  rectype=rectype,
#                                  group=2,
#                                  exclude_series = np.delete(np.arange(0, 24), 10),
#                                  scaling=scaling,
#                                  unit=unit)

# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-023-adaEk.dat',
#                                  rectype=rectype,
#                                  group=4,
#                                  exclude_series = np.delete(np.arange(0, 48), 9),
#                                  scaling=scaling,
#                                  unit=unit)

# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-023-adaEk.dat',
#                                  rectype=rectype,
#                                  group=5,
#                                  exclude_series = np.delete(np.arange(0, 49), 12),
#                                  scaling=scaling,
#                                  unit=unit)

# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-025-adaEk_TTX.dat',
#                                   rectype=rectype,
#                                   group=1,
#                                   exclude_series = np.delete(np.arange(0, 27), 13),
#                                   scaling=scaling,
#                                   unit=unit)

# trace = MiniTrace.from_heka_file(filename='Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/Syn-AMY-019-TTX_Cd.dat',
#                                  rectype=rectype,
#                                  group=1,
#                                  exclude_series = np.delete(np.arange(0, 13), 10),
#                                  scaling=scaling,
#                                  unit=unit)


b, a = sc.signal.bessel(3, 750, fs=100e3)

ori_trace = trace.data

data_filtered = sc.signal.lfilter(b, a, ori_trace)


# %% remove initial or end steps




# %%


# trace.plot_trace()

trace = MiniTrace(data=data_filtered,
                  sampling_interval=1/100e3,
                  y_unit='pA',
                  filename='None')

trace.plot_trace()

# %%

factor = 19
model_th = 0.5
win_size = 600 * factor
direction = 'negative'

# set settings for event detection
eventdetection_settings = {'window_size' : 600 * factor,
                           'model_threshold' : model_th,
                           'batch_size' : 512,
                           'event_detection_peakw' : 5,
                           'stride' : 30,
                           'rel_prom_cutoff' : 0.25,
                           'convolve_win' : 20 * factor}

# run prediction
detection = EventDetection(data=trace,
                           model_path='C:/Users/nesseler/miniML/models/GC_lstm_model.h5',
                           window_size = eventdetection_settings['window_size'],
                           model_threshold = eventdetection_settings['model_threshold'],
                           batch_size = eventdetection_settings['batch_size'],
                           event_direction=direction,
                           compile_model=True,
                           verbose=2)

# detect events
detection.detect_events(eval=True,
                        stride = eventdetection_settings['stride'],
                        peak_w = eventdetection_settings['event_detection_peakw'],
                        rel_prom_cutoff = eventdetection_settings['rel_prom_cutoff'],
                        convolve_win = eventdetection_settings['convolve_win'],
                        resample_to_600 = True)

MiniPlots = miniML_plots(data=detection)

MiniPlots.plot_prediction(include_data=True, plot_filtered_prediction=True,
                          plot_filtered_trace=True, plot_event_params=True)
MiniPlots.plot_event_overlay()
MiniPlots.plot_event_histogram(plot='amplitude', cumulative=False)
# MiniPlots.plot_gradient_search()

# %%


# MiniPlots.plot_event_overlay()


# plt.gcf().savefig('C:/Users/nesseler/Desktop/test.png', format = 'png')


# %% create dataframe

# # get number of events
# n_events = detection.event_locations.shape[0]


# events = pd.DataFrame(columns=['event_locations', 'event_scores', 'event_peak_locations', 'event_peak_times'],
#                       index=np.arange(n_events))
# events.index.name = 'event_id'

# events = events.assign(event_locations=detection.event_locations,
#                        event_scores=detection.event_scores,
#                        event_peak_locations=detection.event_peak_locations,
#                        event_peak_times=detection.event_peak_times)

# %% figure


# for event_i in range(21):
#     event_idx = events.iat[event_i, 2]
#     event_win_idc = np.arange(event_idx-1000, event_idx+6000, dtype=int)
#     start_point = detection.start_pnts[event_i] - event_idx
#     end_point = detection.end_pnts[event_i] - event_idx

#     plt.title(str(event_idx) + '\n amplitude ' + str(round(detection.event_stats.amplitudes[event_i], 2)) + '\n risetime ' + str(
#         round(detection.event_stats.risetimes[event_i]*1e3, 2)) + '\n score ' + str(round(detection.event_stats.event_scores[event_i], 2)))
#     plt.plot(detection.prediction[event_win_idc])
#     # plt.plot(ori_trace[event_win_idc])
#     plt.plot(data_filtered[event_win_idc])
#     plt.scatter(1000, data_filtered[event_idx])

#     # plt.plot(detection.events[event_i])
#     # plt.scatter(events.iat[event_i, 0] - event_idx +1000, data_filtered[events.iat[event_i, 0]])

#     plt.show()
#     # plt.plot(detection.smth_gradient[event_win_idc])
#     # plt.axhline(detection.grad_threshold)
#     # plt.show()

# %%

# print(detection.event_bsls)
# print(detection.bsl_starts)
# print(detection.event_start_times)
# print(detection.min_positions_rise)
# print(detection.half_decay)
# print(detection.decaytimes)


# # onset!
# print(detection.event_start)

# %%

# for event_i in range(21):
#     plt.plot(detection.events[event_i])
#     plt.show()


# %%


def exp_fit(x: np.ndarray, amp: float, tau: float, offset: float) -> np.ndarray:
    """
    Fits an exponential curve to the given data.

    Parameters:
        x (np.ndarray): The input data.
        amp (float): The amplitude of the exponential curve.
        tau (float): The time constant of the exponential curve.
        offset (float): The offset of the exponential curve.

    Returns:
        np.ndarray: The fitted exponential curve.
    """

    return amp * np.exp(-(x - x[0]) / tau) + offset


def round_to_base(number, base):
    return base * round(number/base)


def round_up_to_base(number, base):
    return base * np.ceil(number/base)


def round_down_to_base(number, base):
    return base * np.floor(number/base)


# %%


filtered_prediction = maximum_filter1d(
    detection.prediction, size=int(5*detection.interpol_factor), origin=-2)


def create_single_event_fig(parent_fig, event_i: int = 0):

    t_winsize = win_size / (SR/1e3)
    t_preevent = t_winsize * 0.35 # ms
    t_postevent = t_winsize * 0.75  # ms
    peak_idx = detection.event_peak_locations[event_i]
    start_idx = int(peak_idx - (t_preevent * SR / 1e3))
    stop_idx = int(peak_idx + (t_postevent * SR / 1e3))
    event_idc = np.arange(start_idx, stop_idx, 1, dtype=int)

    x_event = np.arange(0 - t_preevent, 0 + t_postevent, step=1/(SR/1e3))

    # shift x axis for prediction
    t_shift_winsize = win_size / (SR/1e3)
    event_idc_shifted = np.arange(
        start_idx-(win_size/2), stop_idx-(win_size/2), 1, dtype=int)

    # get half decay
    decay_idx = int(detection.half_decay[event_i])
    decay_time = detection.decaytimes[event_i] * 1e3

    # fit exponential function for tau estimation
    points_after = int(win_size/2) # + (win_size/3))

    decay_idc = np.arange(peak_idx, peak_idx+points_after, 1, dtype=int)

    decay = data_filtered[decay_idc]

    decay_x = np.arange(0, decay.shape[0] / (SR/1e3), step=1/(SR/1e3))

    peak = data_filtered[peak_idx]

    fit, _ = curve_fit(exp_fit, decay_x, -decay,
                       p0=[-peak, decay_time, 0],
                       bounds=([0, 0, -np.inf], [np.inf, 1e3, 1e3]))

    fitted_decay = -exp_fit(decay_x, *fit)

    axs = parent_fig.subplots(nrows=2,
                              ncols=1,
                              sharex=True,
                              height_ratios=[1, 3])

    # plot prediction
    axs[0].plot(x_event, detection.prediction[event_idc],
                color='k',
                alpha=0.5,
                lw=0.75,
                label='unfiltered')

    axs[0].plot(x_event, filtered_prediction[event_idc],
                color='k',
                alpha=0.5,
                lw=0.75,
                ls='dashed',
                label='filtered')

    axs[0].plot(x_event, filtered_prediction[event_idc_shifted],
                color='k',
                lw=0.75,
                label='filtered &\nshifted (winsize/2)')

    axs[0].legend(loc='upper right',
                  fontsize=7,
                  frameon=False)

    # plot data
    axs[1].plot(x_event, data_filtered[event_idc],
                lw=0.75,
                color='k',
                label='data (filt.)')

    # plot baseline
    axs[1].hlines(xmin=(detection.bsl_starts[event_i] - peak_idx) / (SR/1e3),
                  xmax=(detection.bsl_ends[event_i] - peak_idx) / (SR/1e3),
                  y=detection.event_bsls[event_i],
                  color='r',
                  lw=1,
                  label='baseline')

    # plot event onset
    axs[1].scatter((detection.event_start[event_i] - peak_idx) / (SR/1e3),
                   data_filtered[detection.event_start[event_i]],
                   color='r',
                   lw=1,
                   facecolors='none',
                   label='onset')

    # plot rise
    axs[1].plot([detection.min_positions_rise[event_i] * 1e3, detection.max_positions_rise[event_i]*1e3] - (peak_idx / (SR/1e3)),
                [detection.min_values_rise[event_i],
                    detection.max_values_rise[event_i]],
                ls='dashed',
                color='r',
                lw=1,
                label='rise (10-90%)')

    # plot rise start and stop indicators
    axs[1].hlines(xmin=[detection.min_positions_rise[event_i] * 1e3, detection.max_positions_rise[event_i]*1e3] - (peak_idx / (SR/1e3)) - 0.5,
                  xmax=[detection.min_positions_rise[event_i] * 1e3,
                        detection.max_positions_rise[event_i]*1e3] - (peak_idx / (SR/1e3)) + 0.5,
                  y=[detection.min_values_rise[event_i],
                      detection.max_values_rise[event_i]],
                  color='r',
                  lw=1,
                  label='_nolegend_')

    # indicate peak
    axs[1].scatter(0, data_filtered[peak_idx],
                   color='r',
                   label='peak')

    # plot half decay
    axs[1].scatter(decay_time, data_filtered[decay_idx],
                   color='r',
                   marker='x',
                   s=15,
                   label='half decay')

    # plot fitted decay
    axs[1].plot(decay_x, fitted_decay,
                color='r',
                ls='dotted',
                lw=1,
                label='fitted decay')

    # plot window size indicator
    axs[1].hlines(xmin = (detection.bsl_starts[event_i] - peak_idx) / (SR/1e3),
                  xmax = (detection.bsl_starts[event_i] - peak_idx) / (SR/1e3) + (win_size / (SR/1e3)),
                  y=np.max(data_filtered[event_idc]),
                  color='gray',
                  alpha=0.5,
                  lw=4,
                  label='moving window')

    # legend
    axs[1].legend(loc='lower right',
                  fontsize=7,
                  frameon=True,
                  facecolor = 'w',
                  ncols=1)

    # axis
    # x
    axs[1].set_xlabel('Time [ms]', fontsize=9)

    xmin_r = round_down_to_base(-t_preevent, 5)
    xmax_r = round_up_to_base(t_postevent, 5)
    xrange = xmax_r - xmin_r
    axs[0].set_xlim([xmin_r-(xrange*0.01), xmax_r+(xrange*0.01)])
    [ax.spines['bottom'].set_bounds([xmin_r, xmax_r]) for ax in axs]

    # y
    axs[0].set_ylim([-0.05, 1.05])
    axs[0].spines['left'].set_bounds([0, 1])
    axs[0].set_ylabel('Probability', fontsize=9)

    axs[1].set_ylabel('Current [pA]', fontsize=9)
    ymin_r = round_down_to_base(np.min(data_filtered[event_idc]), 1)
    ymax_r = round_up_to_base(np.max(data_filtered[event_idc]), 1)
    yrange = ymax_r - ymin_r
    axs[1].set_ylim([ymin_r-(yrange*0.05), ymax_r+(yrange*0.05)])
    axs[1].spines['left'].set_bounds([ymin_r, ymax_r])

    # add text label of measurements
    event_meas = {'Score: ': detection.event_scores[event_i],
                  'Amplitude [pA]: ': detection.event_stats.amplitudes[event_i],
                  'Risetime [ms]: ': detection.risetimes[event_i] * 1e3,
                  'Half decay time [ms]: ': detection.decaytimes[event_i] * 1e3,
                  'Tau [ms]: ': fit[0]}

    text_str = ''
    for k, v in event_meas.items():
        text_str = text_str + k + str(round(v, 2)) + '\n'

    axs[1].text(x=xmin_r,
                y=ymin_r,
                s=text_str,
                fontsize=7,
                va='bottom',
                ha='left')

    # remove spines
    [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

    # align labels
    parent_fig.align_labels()


# %%

# fig = plt.figure(layout='constrained',
#                  figsize=(9, 8),
#                  dpi=600)

# fig.suptitle('Detected events')

# subfigs = fig.subfigures(2, 2, wspace=0.07)

# # flatten subfigs array
# subfigs = subfigs.flatten()

# # set titles
# titles = ['A', 'B', 'C', 'D']

# for subfig_i, subfig in enumerate(subfigs):
#     subfig.suptitle(titles[subfig_i], x=0.1, ha='left')

# create_single_event_fig(parent_fig=subfigs[0], event_i=9)
# create_single_event_fig(parent_fig=subfigs[1], event_i=6)
# create_single_event_fig(parent_fig=subfigs[2], event_i=60)
# create_single_event_fig(parent_fig=subfigs[3], event_i=94)


# plt.show()


# %%


# for i in range(30):
    
#     try: 
#         fig = plt.figure(layout='constrained',
#                           dpi=600)
        
#         create_single_event_fig(parent_fig=fig, event_i=i)
        
        
#         plt.show()
#     except ValueError:
#         print(i)
        
# %%

# fig = plt.figure(layout='constrained',
#                  dpi=600)

# create_single_event_fig(parent_fig=fig, event_i=1)


# plt.show()

# %%

# parent_fig = plt.figure(layout='constrained',
#                         figsize = (6, 3),
#                         dpi=300)

# # parent_fig.suptitle('Full trace')


# x_full = np.arange(0, data_filtered.shape[0] / SR, step=1/SR)
# x_predic = np.arange(0, filtered_prediction.shape[0] / SR, step=1/SR)

# axs = parent_fig.subplots(nrows=2,
#                           ncols=1,
#                           sharex=True,
#                           height_ratios=[1, 3])

# # plot prediction
# axs[0].plot(x_predic, filtered_prediction,
#             color='k',
#             alpha=0.5,
#             lw=0.5,
#             label='filtered')

# axs[0].hlines(xmin = 0, xmax = x_predic[-1],
#               y = model_th, 
#               color = 'r',
#               lw = 0.5,
#               ls = 'dashed')

# # plot data
# # axs[1].plot(x_full, ori_trace,
# #             color='k',
# #             alpha = 0.5,
# #             lw=0.5,
# #             label='data')

# # plot data
# axs[1].plot(x_full, data_filtered,
#             color='k',
#             lw=0.5,
#             label='data')


# # plot event peak indicators
# axs[1].scatter(x = x_full[detection.event_peak_locations],
#                y = data_filtered[detection.event_peak_locations],
#                marker = 'o',
#                color = 'r',
#                s = 5,
#                lw = 1)

# # plot eventplot
# axs[1].eventplot(positions = detection.event_peak_locations / SR,
#                  orientation = 'horizontal',
#                  lineoffsets = 9.25,
#                  linelengths = 1.5,
#                  linewidths = 1,
#                  color = 'r',
#                  label = 'events')



# # remove spines
# [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# # axis 
# axs[0].set_ylim([-0.05, 1.05])
# axs[0].set_ylabel('Probability')

# t = 2
# axs[1].set_xlim([0+t, 2+t])
# axs[1].set_xlabel('Time [s]')

# axs[1].set_ylim([-20, 10])
# axs[1].set_ylabel('Current [pA]')

# parent_fig.align_labels()

# plt.show()


# %%


# for t in np.arange(0, 180, 1):
    
t = 1
trange = 179

ext_prediction = np.append(filtered_prediction, [np.nan]*(data_filtered.shape[0]-filtered_prediction.shape[0]))


parent_fig = plt.figure(layout='constrained',
                        figsize = (6, 3),
                        dpi=300)

parent_fig.suptitle(f'Prediction - Trace - Events',
                    fontsize = 9)


x_full = np.arange(0, data_filtered.shape[0] / SR, step=1/SR)

plt_idc = np.arange((0+t) * SR, (trange+t) * SR, step = 1, dtype = int)
plt_peak_idc = [idx for idx in detection.event_peak_locations if (idx > plt_idc[0] and idx < plt_idc[-1])]

axs = parent_fig.subplots(nrows=3,
                          ncols=1,
                          sharex=True,
                          height_ratios=[1, 0.3, 2])

# plot prediction
axs[0].plot(x_full[plt_idc], ext_prediction[plt_idc],
            color='k',
            alpha=0.5,
            lw=0.5,
            label='filtered')

axs[0].hlines(xmin = x_full[plt_idc][0], xmax = x_full[plt_idc][-1],
              y = model_th, 
              color = 'r',
              lw = 0.5,
              ls = 'dashed')

# plot eventplot
axs[1].eventplot(positions = x_full[plt_peak_idc],
                 orientation = 'horizontal',
                 lineoffsets = 0,
                 linelengths = 1,
                 linewidths = 1,
                 color = 'r',
                 label = 'events')

# plot data
axs[2].plot(x_full[plt_idc], data_filtered[plt_idc],
            color='k',
            lw=0.5,
            label='data (filt)')

# plot event peak indicators
axs[2].scatter(x = x_full[plt_peak_idc],
               y = data_filtered[plt_peak_idc],
               marker = 'o',
               color = 'r',
               s = 5,
               lw = 1)


# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# axis 
axs[0].set_ylim([-0.05, 1.05])
axs[0].set_ylabel('Probability')

axs[1].set_ylim([-1, 1])
axs[1].set_ylabel('Events')
axs[1].set_yticks(ticks = [0], labels = [])

axs[2].set_xlim([0+t-(trange*0.01), trange+t+(trange*0.01)])
[ax.spines['bottom'].set_bounds([0+t, trange+t]) for ax in axs]
axs[2].set_xlabel('Time [s]')

# axs[1].set_ylim([-20, 10])
axs[2].set_ylabel('Current [pA]')

parent_fig.align_labels()

plt.show()


