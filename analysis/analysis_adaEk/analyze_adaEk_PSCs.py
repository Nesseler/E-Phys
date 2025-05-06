# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 19:02:06 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir

# custom functions
from functions.functions_import import get_vc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_filter import butter_filter


PGF = 'vc-Erest-3min'

conditions = ['ctrl', 'adaEk']

SR = 100000

t = np.arange(0, (6*30), 1/SR)


# %% load event detection

import pickle
import gc

# init plotting
from functions.initialize_plotting import *


# import directories 
from parameters.directories_win import synaptic_dir
# miniML_path = synaptic_dir + '/miniML_dtc-validation'

# set cell_ID
cell_ID = 'E-317'
treatments = ['ctrl', 'adaEk']

# set dicts
traces = dict.fromkeys(['ctrl', 'adaEk'])
amplitudes = dict.fromkeys(['ctrl', 'adaEk'])
amplitudes_hist = dict.fromkeys(['ctrl', 'adaEk'])
risetimes = dict.fromkeys(['ctrl', 'adaEk'])
events = dict.fromkeys(['ctrl', 'adaEk'])
event_locations = dict.fromkeys(['ctrl', 'adaEk'])

# set histogram bins
bins = np.arange(-50, 0+1, 1)

# set x axis
x = np.arange(0, 180, 1/SR)

for treatment in ['ctrl', 'adaEk']:

    # set filename
    filename = f'miniMLdetect_{cell_ID}_Erest_3min_{treatment}' 
    
    # open a file, where you stored the pickled data
    file = open((synaptic_dir + f'/miniML_dtc-Erest-{treatment}/' + filename + '.pickle'), 'rb')
    
    # dump information to that file
    detection = pickle.load(file)
    
    # close and remove (from memory) the file
    file.close()
    del file 
    gc.collect()
    
    # write to dicts
    traces[treatment] = detection['mini_trace']
    events[treatment] = detection['events']
    amplitudes[treatment] = detection['individual_values']['amplitudes']
    risetimes[treatment] = detection['individual_values']['risetimes']
    event_locations[treatment] = detection['event_location_parameters']['event_peak_locations']
    
    # calc occurrances
    amplitudes_hist[treatment], _ = np.histogram(a = detection['individual_values']['amplitudes'], bins = bins)



# %% create figure per cell

adaEk_color = {'ctrl' : colors_dict['primecolor'],
               'adaEk' : 'goldenrod'}

# init figure
fig, axs = plt.subplot_mosaic(mosaic = [[0, 0, 0, 1, 1, 1], [2, 2, 2, 3, 3, 3], [4, 4, 4, 4, 5, 5], [4, 4, 4, 4, 6, 6]],#'AACC;BBDD;EEEF;EEEG', 
                             dpi = 300, 
                             layout = 'constrained',
                             figsize = get_figure_size(width = 159.2, height = 130),
                             height_ratios = [3, 0.5, 3, 3],
                             width_ratios =  [1, 1, 1, 1, 1, 1])

fig.set_constrained_layout_pads(w_pad=0./72., h_pad=4./72.,
                                hspace=0./72., wspace=0./72.)

# figure title
fig.suptitle(f'{cell_ID} - PSCs - adapted Ek')

axs_titles = {0 : 'A: PSCs ctrl', 1 : 'B: PSCs adaEk',
              2 : 'C: Events ctrl', 3: 'D: Events adaEk',
              4 : 'C: Ampltiude histogram',
              5 : 'D: Events ctrl', 6 : 'E: Events adaEk'}


# traces
for ti, treatment in enumerate(treatments):
    
    # traces
    ax = axs[ti]
    
    # title
    ax.set_title(axs_titles[ti], loc = 'left', fontsize = 9)
    
    # plot prediction
    ax.plot(x, traces[treatment],
            color= adaEk_color[treatment],
            lw=0.5,
            label= 'trace', 
            zorder = 1)
    
    ax.scatter(x[event_locations[treatment]], traces[treatment][event_locations[treatment]],
               color = 'r',
               alpha = 0.5,
               s = 3,
               lw=0.5,
               label='events',
               zorder = 0)
    
    ax.legend(frameon = False,
              fontsize = 9,
              loc = 'lower right',
              labelspacing = 0.1,
              borderpad = 0.1,
              borderaxespad = 0.1)
    
    # y axis
    ax.set_ylim([-40-0.4, 10+0.4])
    ax.spines['left'].set_bounds([-40, 10])
    ax.set_yticks(ticks = np.arange(-40, 10+1, 20))
    ax.set_yticks(ticks = np.arange(-40, 10+1, 5), minor = True)
        
    # Eventplots
    ax = axs[ti+2]

    # plot eventplot
    ax.eventplot(positions = x[event_locations[treatment]],
                  orientation = 'horizontal',
                  lineoffsets = 0,
                  linelengths = 1,
                  linewidths = 0.5,
                  color = adaEk_color[treatment],
                  label = 'events')

axs[0].set_ylabel('Current [pA]')
axs[2].set_ylabel('Events')

for axi in range(5):
    axs[axi].set_xlim([0-1.8, 180+1.8])
    axs[axi].spines['bottom'].set_bounds([0, 180])
    axs[axi].set_xticks(ticks = np.arange(0, 180+1, 60),
                        labels = [])
    axs[axi].set_xticks(ticks = np.arange(0, 180+1, 10), minor = True)
    

remove_spines_n_ticks([axs[0], axs[1]], axis = 'x')
remove_spines_n_ticks([axs[2], axs[3]], axis = 'y')

for axi in [2, 3]:
    axs[axi].set_yticks(ticks = [])
    axs[axi].set_xticks(ticks = np.arange(0, 180+1, 60),
                        labels = np.arange(0, 180+1, 60))
    axs[axi].set_xlabel('Time [ms]')

# set histogram axis
ax = axs[4]

# title
ax.set_title(axs_titles[4], loc = 'left', fontsize = 9)

ax.stairs(amplitudes_hist['ctrl'], bins, 
          fill = False,
          lw = 1,
          alpha = 0.7,
          color = colors_dict['primecolor'],
          label = 'ctrl')

# plot treatment
ax.stairs(amplitudes_hist['adaEk'], bins, 
          fill = False,
          lw = 1,
          color = 'goldenrod',
          label = 'adaEk')
    
ax.legend(title = 'Treatment', 
          frameon = False,
          loc = 'upper left',
          fontsize = 9)

# edit axis
ax.set_ylabel('Event count [#]')
ax.set_ylim([-2, 132])
ax.set_yticks(ticks = np.arange(0, 130+1, 20))
ax.set_yticks(ticks = np.arange(0, 130+1, 5), minor = True)
ax.spines['left'].set_bounds([0, 130])

ax.set_xlabel('Amplitude [pA]')
ax.set_xticks(ticks = np.arange(-50, 0+1, 10),
              labels = np.arange(-50, 0+1, 10))
ax.set_xticks(ticks = np.arange(-50, 0, 1), minor = True)
ax.set_xlim([-51, 1])
ax.spines['bottom'].set_bounds([-50, 0])


# average events

# traces
for ti, treatment in enumerate(treatments):
    
    # axis
    axi = ti+5
    ax = axs[axi]
    ax.set_title(axs_titles[axi], loc = 'left', fontsize = 9)

    # set data
    event_x = np.arange(0, events[treatment].shape[1]) * 1/(SR/1e3)
    event_average = np.mean( events[treatment], axis=0)
    
    # plot
    ax.plot(event_x, events[treatment].T,
            c = 'k',
            alpha = 0.3,
            lw = 1,
            label = '_nolegend_')
    
    ax.plot(event_x, event_average,
            c = 'r', 
            lw = 1.5,
            label = 'average event')
    
    ax.legend(frameon = False,
               fontsize = 9,
               loc = 'lower right')
    
    # edit axis
    ax.set_ylabel('Current [pA]')
    ax.set_ylim([-41, 10])
    ax.set_yticks(ticks = np.arange(-40, 10+1, 20))
    ax.set_yticks(ticks = np.arange(-40, 10+1, 5), minor = True)
    ax.spines['left'].set_bounds([-40, 10])
    
    ax.set_xlabel('Time [ms]')
    ax.set_xticks(ticks = np.arange(0, round(event_x[-1])+1, 50))
    ax.set_xticks(ticks = np.arange(0, round(event_x[-1])+1, 10), minor = True)
    ax.set_xlim([-5, round(event_x[-1])+5])
    ax.spines['bottom'].set_bounds([0, round(event_x[-1])])


# align labels
fig.align_labels()

# remove spines
[axs[i].spines[spine].set_visible(False) for i in range(7) for spine in ['top', 'right']]

# display figure
plt.show()

# %%

def create_uniform_events(n_events: int = 100, t: float = 180) -> np.array:
    """
    

    Parameters
    ----------
    n_events : int, optional
        DESCRIPTION. The default is 100.
    t : float, optional
        DESCRIPTION. The default is 180.

    Returns
    -------
    uni_events : TYPE
        DESCRIPTION.

    """
    # uniformly distribute the same number of events as measured
    # by randomly draw from uniform distribution between t0 and t_end
    uni_events = sc.stats.uniform.rvs(loc = 0, scale = t, size = n_events)
    
    # sort values to perserve time series aspect
    uni_events = np.sort(uni_events)
        
    return uni_events


treatment = 'ctrl'
total_dur = 180
n_events = event_locations['ctrl'].shape[0]
uniform_events = create_uniform_events(n_events, t=total_dur)
uniform_hist, _ = np.histogram(np.diff(uniform_events, prepend = 0), bins = bins)


def calc_PDF_fromTimepoints(data: np.array, bins: np.array):

    # calc inter-event interval (IEI)
    IEI = np.diff(data, prepend = 0)    

    # calc mean rate
    mean_rate = 1 / np.mean(IEI)
    
    # calculate the probability density function with bins given
    pdf = mean_rate * np.exp(- mean_rate * bins)

    return pdf


def calc_CDF_fromTimepoints(data: np.array, bins: np.array):
    
    pdf = calc_PDF_fromTimepoints(data, bins)

    cdf = [np.max(pdf) - p for p in pdf]

    return cdf


def minmax_norm(data: np.array, minzero: bool = False):
    
    
    # calc min max normalized version of array 
    if minzero:
        normed = (data - 0) / (np.max(data) - 0)
    else:
        normed = (data - np.min(data)) / (np.max(data)- np.min(data))
    
    return normed
    

pdf = calc_PDF_fromTimepoints(uniform_events, bins)
cdf = calc_CDF_fromTimepoints(uniform_events, bins)

# plt.plot(bins, pdf)
# plt.plot(bins, cdf)
# # plt.stairs(minmax_norm(uniform_hist), bins)


# %%

# plt.eventplot(positions = x[event_locations['ctrl']],
#               orientation = 'horizontal',
#               lineoffsets = 0,
#               linelengths = 1,
#               linewidths = 0.5,
#               color = adaEk_color['ctrl'],
#               label = 'events')

# plt.eventplot(positions = create_uniform_events(event_locations['ctrl'].shape[0]),
#               orientation = 'horizontal',
#               lineoffsets = 1.5,
#               linelengths = 1,
#               linewidths = 0.5,
#               color = adaEk_color['ctrl'],
#               label = 'events')

# plt.xlim([0,180])

# %% frequency histogram test

treatment = 'ctrl'
# total_dur = 180
# n_events = event_locations['ctrl'].shape[0]
# uniform_events = create_uniform_events(n_events, t=total_dur)


# set data and bins
# data = amplitudes['ctrl']

data = np.diff(event_locations[treatment] * 1/SR, prepend=0)
data_shuffled = np.diff(uniform_events, prepend=0)
bins = np.arange(0, 30+0.05, 0.05)

# for figure
trace = traces[treatment]
event_locations = event_locations[treatment]
x = np.arange(0, trace.shape[0] / SR, 1/SR)

n_events = event_locations.shape[0]
uniform_events = create_uniform_events(n_events, t=trace.shape[0] / SR)

x_max = trace.shape[0] / SR
trace_mean = np.mean(trace)

events_hist, _ = np.histogram(np.diff(event_locations / SR, prepend = 0), bins = bins)
uniform_hist, _ = np.histogram(np.diff(uniform_events, prepend = 0), bins = bins)
theo_pdf = calc_PDF_fromTimepoints(data = uniform_events, bins = bins)
theo_cdf = calc_CDF_fromTimepoints(data = uniform_events, bins = bins)


# %% 

from functions.functions_useful import round_up_to_base, round_down_to_base

# %%

# init figure
fig, axs = plt.subplots(nrows = 3,
                        ncols = 1, 
                        dpi = 300, layout = 'constrained', 
                        height_ratios = [2, 1, 7])

# figure title
fig.suptitle(f'{cell_ID} PSCs - Event frequency (IEI)')

# trace
ax = axs[0]

# title
# ax.set_title(axs_titles[ti], loc = 'left', fontsize = 9)

# plot prediction
ax.plot(x, trace,
        color = 'k',
        lw=0.5,
        label= 'trace', 
        zorder = 1)

ax.scatter(x[event_locations], trace[event_locations],
           color = 'r',
           alpha = 0.5,
           s = 3,
           lw=0.5,
           label='events',
           zorder = 0)

ax.legend(frameon = False,
          fontsize = 9,
          loc = 'lower right',
          labelspacing = 0.1,
          borderpad = 0.1,
          borderaxespad = 0.1)

apply_axis_settings(ax, axis = 'y',
                    ax_min = round_down_to_base(trace_mean-50, 10),
                    ax_max = round_up_to_base(trace_mean+10, 10),
                    pad = None,
                    step = 20,
                    stepminor = 5,
                    limits_n_0 = True,
                    label = 'Current [pA]')


# eventplots
ax = axs[1]

ax.eventplot(positions = x[event_locations],
              orientation = 'horizontal',
              lineoffsets = 2,
              linelengths = 1,
              linewidths = 0.5,
              color = 'k',
              label = 'events')

ax.eventplot(positions = uniform_events,
              orientation = 'horizontal',
              lineoffsets = 0,
              linelengths = 1,
              linewidths = 0.5,
              color = 'gray',
              alpha = 0.5,
              label = 'uniform')

for ax in axs[:2]:
    ax.set_xlim([0-1.8, x_max+1.8])
    ax.spines['bottom'].set_bounds([0, x_max])
    ax.set_xticks(ticks = np.arange(0, x_max+1, 60),
                        labels = [])
    ax.set_xticks(ticks = np.arange(0, x_max+1, 10), minor = True)


remove_spines_n_ticks([axs[0]], axis = 'x')
axs[1].set_xticks(ticks = np.arange(0, x_max+1, 60),
                  labels = np.arange(0, x_max+1, 60, dtype = int))
axs[1].set_xlabel('Time [s]')

axs[1].set_yticks(ticks = [0, 2], labels = ['uniform', 'events'])
axs[1].spines['left'].set_visible(False)


# histogram
ax = axs[2]

ax.stairs(minmax_norm(events_hist), bins, 
          fill = False,
          lw = 1,
          color = colors_dict['primecolor'],
          label = 'IEI')

ax.plot(bins, minmax_norm(theo_pdf),
        color = 'k',
        lw = 1,
        ls = 'dashed',
        label = 'uniform pdf')

ax.plot(bins, minmax_norm(theo_cdf),
        color = 'k',
        lw = 1,
        ls = 'dotted',
        label = 'uniform cdf')

ax.legend(frameon = False, loc = 'lower right', fontsize = 9)

ax.set_xlabel('IEI [s]')
ax.set_xticks(ticks = np.arange(0, 5+0.1, 1))
ax.set_xticks(ticks = np.arange(0, 5+0.1, 0.2), minor = True)
ax.set_xlim([0-0.1, 5+0.1])
ax.spines['bottom'].set_bounds([0, 5])

ax.set_ylabel('Normalized event count')
ax.set_yticks(ticks = np.arange(0, 1+0.01, 0.2))
ax.set_yticks(ticks = np.arange(0, 1+0.01, 0.05), minor = True)
ax.set_ylim([0-0.01, 1+0.01])
ax.spines['left'].set_bounds([0, 1])


# inset axis: full histogram
ax_inset = ax.inset_axes([3/5, 0.55, 1.8/5, 0.35])

ax_inset.stairs(events_hist, bins, 
                fill = False,
                lw = 1,
                color = colors_dict['primecolor'],
                label = 'IEI')

# inset marker
box_xmin   = 0
box_width  = 5 
box_ymin   = 0
box_height = round_up_to_base(np.max(events_hist), 5)

# add rectangle marker
ax_inset.add_patch(Rectangle(xy = (box_xmin, box_ymin), 
                   width = box_width, 
                   height = box_height,
                   fill = False,
                   color = 'r',
                   linestyle = '--',
                   lw = 0.5))

ax_inset.set_ylabel('Count [#]')

ax_inset.set_xlabel('IEI [s]')
ax_inset.set_xticks(ticks = np.arange(0, bins[-1]+0.1, 10))
ax_inset.set_xticks(ticks = np.arange(0, bins[-1]+0.1, 1), minor = True)
ax_inset.set_xlim([0-0.1, bins[-1]+0.1])
ax_inset.spines['bottom'].set_bounds([0, bins[-1]])
[ax_inset.spines[spine].set_visible(False) for spine in ['top', 'right']]

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# display figure
plt.show()

















