# -*- coding: utf-8 -*-
"""
Created on Wed May  7 14:35:33 2025

@author: nesseler
"""

# import packages
from functions.initialize_packages import *
import pickle
import gc

# import directories 
from parameters.directories_win import synaptic_dir, figure_dir

# import functions


# set cell_ID
cell_ID = 'E-304'

# set filename
filename = f'miniMLdetect_{cell_ID}_Erest_3min_ctrl' 

# open a file, where you stored the pickled data
file = open((synaptic_dir + f'/miniML_dtc-Erest-ctrl/' + filename + '.pickle'), 'rb')

# dump information to that file
detection = pickle.load(file)

# close and remove (from memory) the file
file.close()
del file 
gc.collect()

## trace
trace = detection['mini_trace']

## prediction 
prediction = detection['prediction']

## prediction threshold
th = detection['metadata']['miniml_model_threshold']

# filter prediction
# interpolation factor


# extend prediction with nan values
prediction = np.append(prediction, [np.nan]*(trace.shape[0]-prediction.shape[0]))

## event indices
event_peak_idc = detection['event_location_parameters']['event_peak_locations']


SR = 100_000
x = np.arange(0, trace.shape[0] / SR, 1/SR)


# %% plotting

# init plotting
from functions.initialize_plotting import *


fig, axs = plt.subplot_mosaic(mosaic = [[0, 0, 0], 
                                        [1, 1, 1], 
                                        [2, 2, 2], 
                                        [3, 6, 9], 
                                        [4, 7, 10], 
                                        [5, 8, 11], 
                                        [12, 15, 18], 
                                        [13, 16, 19], 
                                        [14, 17, 20]],
                              layout = 'constrained',
                              figsize = get_figure_size(width = 159.2, height = 140),
                              dpi = 300,
                              height_ratios = [4, 1.2, 5, 1, 0.6, 5, 1, 0.6, 5])

# set start axis
ax_start = 0
t = 0
trange = 180
ymin = -45
ymax = 10

# get plotting indices
plt_idc = np.arange((0+t) * SR, (trange+t) * SR, step = 1, dtype = int)
    
# limit event peak indices
plt_peak_idc = [idx for idx in event_peak_idc if (idx > plt_idc[0] and idx < plt_idc[-1])]

# calc y range
ypad = (ymax - ymin) * 0.01

# set axis
ax_pred = axs[ax_start]
ax_event = axs[ax_start+1]
ax_trace = axs[ax_start+2]

# plot prediction
ax_pred.plot(x[plt_idc], prediction[plt_idc],
             color = colors_dict['primecolor'],
             lw = 0.5,
             label = 'prediction')

# plot threshold
ax_pred.hlines(xmin = x[plt_idc][0], xmax = x[plt_idc][-1],
               y = th, 
               color = colors_dict['color3'],
               lw = 0.5,
               ls = 'dashed')

# plot eventplot
ax_event.eventplot(positions = x[plt_peak_idc],
                   orientation = 'horizontal',
                   lineoffsets = 0,
                   linelengths = 1,
                   linewidths = 1,
                   color = colors_dict['color3'],
                   label = 'events')
    
# plot trace
ax_trace.plot(x[plt_idc], trace[plt_idc],
              color='k',
              lw=0.5,
              label='data (filt)',
              zorder = 2)
    
# plot event peak indicators
ax_trace.scatter(x = x[plt_peak_idc],
                 y = trace[plt_peak_idc],
                 marker = 'o',
                 color = 'r',
                 s = 5,
                 lw = 1,
                 zorder = 1)

# zeroline
ax_trace.hlines(xmin = 0, xmax = trace.shape[0] / SR, y = 0, 
                zorder = 0, 
                lw = 0.5, 
                color = 'k', 
                ls = 'dashed', 
                alpha = 0.5)


# y axes
ax_pred.set_ylabel('Probability')
ax_pred.set_ylim([-0.05, 1.05])
ax_pred.spines['left'].set_bounds([0, 1])
ax_pred.set_yticks(ticks = [0, 0.5, 1], labels = [0, '', 1])

# ax_event.set_ylabel('Events')
ax_event.set_ylim([-1, 1])
ax_event.spines['left'].set_bounds([-1, 1])
ax_event.set_yticks(ticks = [0], labels = [])

ax_trace.set_ylabel('Current [pA]')
ax_trace.set_ylim([ymin - ypad, ymax + ypad])
ax_trace.spines['left'].set_bounds([ymin, ymax])
ax_trace.set_yticks(ticks = [ymin, 0, ymax], labels = [ymin, '', ymax])
ax_trace.set_yticks(ticks = np.arange(ymin, ymax+1, 5), minor = True)
    
# x axes
xmax = t+trange
xpad = trange * 0.01

ax_pred.set_xlim([t - xpad, xmax + xpad])
ax_pred.spines['bottom'].set_bounds([t, xmax])
ax_pred.set_xticks(ticks = np.arange(t, xmax+0.0001, 0.1), labels = [])
ax_pred.set_xticks(ticks = np.arange(t, xmax+1, 10), minor = True)

ax_event.set_xlim([t - xpad, xmax + xpad])
remove_spines_n_ticks([ax_event], axis = 'x')
ax_event.set_xticks(ticks = [])

ax_trace.set_xlabel('Time [s]')
ax_trace.set_xlim([t - xpad, xmax + xpad])
ax_trace.spines['bottom'].set_bounds([t, xmax])
ax_trace.set_xticks(ticks = np.arange(t, xmax+1, 30))
ax_trace.set_xticks(ticks = np.arange(t, xmax+1, 5), minor = True)


def plot_pred_event_trace(ax_start = 0, t = 0, trange = 180, ymin = -45, ymax = 10, l = 'A'):
    # set start axis
    # ax_start = 0
    # t = 0
    # trange = 180
    # ymin = -45
    # ymax = 10
    
    # get plotting indices
    plt_idc = np.arange((0+t) * SR, (trange+t) * SR, step = 1, dtype = int)
        
    # limit event peak indices
    plt_peak_idc = [idx for idx in event_peak_idc if (idx > plt_idc[0] and idx < plt_idc[-1])]
    
    # calc y range
    ypad = (ymax - ymin) * 0.01
    
    # set axis
    ax_pred = axs[ax_start]
    ax_event = axs[ax_start+1]
    ax_trace = axs[ax_start+2]
    
    # plot prediction
    ax_pred.plot(x[plt_idc], prediction[plt_idc],
                 color = colors_dict['primecolor'],
                 lw = 0.5,
                 label = 'prediction')
    
    # plot threshold
    ax_pred.hlines(xmin = x[plt_idc][0], xmax = x[plt_idc][-1],
                   y = th, 
                   color = colors_dict['color3'],
                   lw = 0.5,
                   ls = 'dashed')
    
    # plot eventplot
    ax_event.eventplot(positions = x[plt_peak_idc],
                       orientation = 'horizontal',
                       lineoffsets = 0,
                       linelengths = 1,
                       linewidths = 1,
                       color = colors_dict['color3'],
                       label = 'events')
        
    # plot trace
    ax_trace.plot(x[plt_idc], trace[plt_idc],
                  color='k',
                  lw=0.5,
                  label='data (filt)',
                  zorder = 2)
        
    # plot event peak indicators
    ax_trace.scatter(x = x[plt_peak_idc],
                     y = trace[plt_peak_idc],
                     marker = 'o',
                     color = 'r',
                     s = 5,
                     lw = 1,
                     zorder = 1)
    
    # zeroline
    ax_trace.hlines(xmin = 0, xmax = trace.shape[0] / SR, y = 0, 
                    zorder = 0, 
                    lw = 0.5, 
                    color = 'k', 
                    ls = 'dashed', 
                    alpha = 0.5)
    
    ax_pred.set_title(l, loc = 'left')
    
    # y axes
    # ax_pred.set_ylabel('Probability')
    ax_pred.set_ylim([-0.05, 1.05])
    ax_pred.spines['left'].set_bounds([0, 1])
    ax_pred.set_yticks(ticks = [0, 0.5, 1], labels = [0, '', 1])
    
    # ax_event.set_ylabel('Events')
    ax_event.set_ylim([-1, 1])
    ax_event.spines['left'].set_bounds([-1, 1])
    ax_event.set_yticks(ticks = [0], labels = [])
    
    # ax_trace.set_ylabel('Current [pA]')
    ax_trace.set_ylim([ymin - ypad, ymax + ypad])
    ax_trace.spines['left'].set_bounds([ymin, ymax])
    ax_trace.set_yticks(ticks = [ymin, 0, ymax], labels = [ymin, '', ymax])
    ax_trace.set_yticks(ticks = np.arange(ymin, ymax+1, 5), minor = True)
        
    # x axes
    xmax = t+trange
    xpad = trange * 0.01
    
    ax_pred.set_xlim([t - xpad, xmax + xpad])
    ax_pred.spines['bottom'].set_bounds([t, xmax])
    ax_pred.set_xticks(ticks = np.arange(t, xmax+0.0001, 0.1), labels = [])
    # ax_pred.set_xticks(ticks = np.arange(t, xmax+1, 10), minor = True)
    
    ax_event.set_xlim([t - xpad, xmax + xpad])
    remove_spines_n_ticks([ax_event], axis = 'x')
    ax_event.set_xticks(ticks = [])
    
    # ax_trace.set_xlabel('Time [s]')
    ax_trace.set_xlim([t - xpad, xmax + xpad])
    ax_trace.spines['bottom'].set_bounds([t, xmax])
    ax_trace.set_xticks(ticks = np.arange(t, xmax+0.0001, 0.05), labels = np.arange(0, (trange+0.0001)*1e3, 50, dtype = int))
    # ax_trace.set_xticks(ticks = np.arange(t, xmax+1, 10), minor = True)
    
    # add indicator to full trace
    axs[2].arrow(x = t+(trange/2), y = 10,
                 dx = 0, dy = -3,
                 head_length = 4,
                 head_width = 1,
                 length_includes_head = True,
                 facecolor = colors_dict['primecolor'])
    
    axs[2].text(x = t+(trange/2)+0.75, y = 15,
                s = l,
                va = 'top', ha = 'left',
                fontsize = 6)


# full
# 'E-304'
plot_pred_event_trace(ax_start = 3, t = 2.85, trange = 0.1, ymin = -30, ymax = 5, l = 'B')
plot_pred_event_trace(ax_start = 6, t = 26.190, trange = 0.100, ymin = -20, ymax = 5, l = 'C')
plot_pred_event_trace(ax_start = 9, t = 35.600, trange = 0.100, ymin = -40, ymax = 5, l = 'D')
plot_pred_event_trace(ax_start = 12, t = 19.15, trange = 0.3, ymin = -40, ymax = 5, l = 'E')
plot_pred_event_trace(ax_start = 15, t = 62.9, trange = 0.1, ymin = -10, ymax = 10, l = 'F')
plot_pred_event_trace(ax_start = 18, t = 121.7, trange = 0.1, ymin = -10, ymax = 5, l = 'G')

# 'E-305'
# plot_pred_event_trace(ax_start = 3, t = 2.8, trange = 0.2, ymin = -40, ymax = 5, l = 'B')
# plot_pred_event_trace(ax_start = 6, t = 22.4, trange = 0.2, ymin = -10, ymax = 5, l = 'C')
# plot_pred_event_trace(ax_start = 9, t = 19., trange = 0.200, ymin = -20, ymax = 5, l = 'D')
# plot_pred_event_trace(ax_start = 12, t = 35.600, trange = 0.100, ymin = -10, ymax = 5, l = 'E')
# plot_pred_event_trace(ax_start = 15, t = 109.4, trange = 0.3, ymin = -10, ymax = 5, l = 'F')
# plot_pred_event_trace(ax_start = 18, t = 144.2, trange = 0.3, ymin = -15, ymax = 5, l = 'G')


# for axi, l in zip([0, 3, 6, 9, 12, 15, 18], ['A', 'B', 'C', 'D', 'E', 'F', 'G']):
axs[0].set_title('A', loc = 'left')

for axi in [5, 8, 11, 14, 17, 20]:
    axs[axi].set_xlabel('Time [ms]')
    
for axi in [3, 12]:
    axs[axi].set_ylabel('Probability')
    axs[axi+2].set_ylabel('Current [pA]')


# align labels
fig.align_labels()

# remove spines
[axs[axi].spines[spine].set_visible(False) for axi in range(21) for spine in ['top', 'right']]

# display figure
plt.show()
      
# save  
save_figures(fig, 
             f'minimL_validation-{cell_ID}_exp_trace', 
             figure_dir + '/miniML_validation/', 
             darkmode_bool, 
             figure_format='both')
