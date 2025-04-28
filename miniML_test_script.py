# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 18:40:21 2025

@author: nesseler
"""

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

%matplotlib inline


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
                                 exclude_series=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26],
                                 scaling=scaling,
                                 unit=unit)


b, a = sc.signal.bessel(3, 1000, fs=100e3)

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
factor = 19

win_size = 600 * factor
direction = 'negative'

detection = EventDetection(data=trace,
                           model_path='C:/Users/nesseler/miniML/models/GC_lstm_model.h5',
                           window_size=win_size,
                           model_threshold=0.5,
                           batch_size=512,
                           event_direction=direction,
                           compile_model=True,
                           verbose=2)

event_detection_peakw = 5

detection.detect_events(eval=True,
                        stride=win_size / 100,
                        peak_w=event_detection_peakw,
                        rel_prom_cutoff=0.25,
                        convolve_win=20 * factor,
                        # gradient_convolve_win = 40 * factor,
                        resample_to_600=True)

MiniPlots = miniML_plots(data=detection)



# %%

# factors = [5, 10, 15, 20, 25, 50, 100]
# n_events = []

# for factor in factors:
    
#     win_size = 600 * factor
#     direction = 'negative'

#     detection = EventDetection(data=trace,
#                                model_path='C:/Users/nesseler/miniML/models/GC_lstm_model.h5',
#                                window_size=win_size,
#                                model_threshold=0.5,
#                                batch_size=512,
#                                event_direction=direction,
#                                compile_model=True,
#                                verbose=2)

#     event_detection_peakw = 5

#     detection.detect_events(eval=True,
#                             stride=win_size / 100,
#                             peak_w=event_detection_peakw,
#                             rel_prom_cutoff=0.25,
#                             convolve_win=20 * factor,
#                             # gradient_convolve_win = 40 * factor,
#                             resample_to_600=True)

#     MiniPlots = miniML_plots(data=detection)
    
#     n_events.append(detection.event_stats.event_count)
    

# %%

MiniPlots.plot_prediction(include_data=True, plot_filtered_prediction=True,
                          plot_filtered_trace=True, plot_event_params=True)
MiniPlots.plot_event_overlay()
MiniPlots.plot_event_histogram(plot='amplitude', cumulative=False)
MiniPlots.plot_gradient_search()


detection.save_to_csv(filename='E-298-test.csv')


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
print(detection.half_decay)
print(detection.decaytimes)


# onset!
print(detection.event_start)

# %%

for event_i in range(21):
    plt.plot(detection.events[event_i])
    plt.show()


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
    t_preevent = t_winsize  # ms
    t_postevent = t_winsize * 1.5  # ms
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
    points_after = int(win_size + (win_size/3))

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
                label='data')

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
    axs[1].hlines(xmin=(detection.bsl_starts[event_i] - peak_idx) / (SR/1e3),
                  xmax=(
                      detection.bsl_starts[event_i] - peak_idx) / (SR/1e3) + (win_size / (SR/1e3)),
                  y=np.max(data_filtered[event_idc]),
                  color='gray',
                  alpha=0.5,
                  lw=4,
                  label='moving window')

    # legend
    axs[1].legend(loc='lower right',
                  fontsize=7,
                  frameon=False,
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
    axs[1].set_ylim([ymin_r-(yrange*0.01), ymax_r+(yrange*0.01)])
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
    [ax.spines[spine].set_visible(False)
     for ax in axs for spine in ['top', 'right']]

    # align labels
    parent_fig.align_labels()


# %%

fig = plt.figure(layout='constrained',
                 figsize=(9, 8),
                 dpi=600)

fig.suptitle('Detected events')

subfigs = fig.subfigures(2, 2, wspace=0.07)

# flatten subfigs array
subfigs = subfigs.flatten()

# set titles
titles = ['A', 'B', 'C', 'D']

for subfig_i, subfig in enumerate(subfigs):
    subfig.suptitle(titles[subfig_i], x=0.1, ha='left')

create_single_event_fig(parent_fig=subfigs[0], event_i=10)
create_single_event_fig(parent_fig=subfigs[1], event_i=5)
create_single_event_fig(parent_fig=subfigs[2], event_i=130)
create_single_event_fig(parent_fig=subfigs[3], event_i=131)


plt.show()


# %%

fig = plt.figure(layout='constrained',
                 dpi=600)

create_single_event_fig(parent_fig=fig, event_i=10)


plt.show()

# %%

%matplotlib qt

parent_fig = plt.figure(layout='constrained',
                        figsize = (6, 3),
                        dpi=300)

parent_fig.suptitle('Full trace')


x_full = np.arange(0, data_filtered.shape[0] / SR, step=1/SR)

axs = parent_fig.subplots(nrows=2,
                          ncols=1,
                          sharex=True,
                          height_ratios=[1, 3])

# plot prediction
axs[0].plot(x_full[:-12000], filtered_prediction,
            color='k',
            alpha=0.5,
            lw=0.5,
            label='filtered')

# plot data
axs[1].plot(x_full, data_filtered,
            color='k',
            lw=0.5,
            label='data')

# plot eventplot
axs[1].eventplot(positions = detection.event_peak_locations / SR,
                 orientation = 'horizontal',
                 lineoffsets = np.min(data_filtered) - 5,
                 linelengths = 2,
                 linewidths = 0.5,
                 color = 'k',
                 label = 'events')


plt.show()


# %%
%matplotlib inline