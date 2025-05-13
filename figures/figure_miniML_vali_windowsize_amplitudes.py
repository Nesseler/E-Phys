# -*- coding: utf-8 -*-
"""
Created on Thu May  8 13:24:49 2025

@author: nesseler
"""


# import packages
from functions.initialize_packages import *

# import directories 
from parameters.directories_win import synaptic_dir, figure_dir

# set cell_ID
cell_IDs = ['E-298', 'E-301', 'E-302', 'E-303', 'E-309', 'E-310', 'E-314']

# cell_IDs = ['E-303']
cell_ID = cell_IDs[0]

# set parameters
winsizes = [36, 96, 114, 276]


# %% dataframes

keys = [f'{cell_ID}_{str(w)}' for cell_ID in cell_IDs for w in winsizes]

amplitudes_hist = pd.DataFrame(columns = keys,
                              index = np.arange(-50, 0+1, 1))

events = dict.fromkeys(keys)
event_peak_locations = dict.fromkeys(keys)
traces = dict.fromkeys(cell_IDs)

# %% load

for cell_ID in tqdm(cell_IDs):
    for winsize in winsizes:
        
        # set filename
        filename = f'miniMLdetect_{cell_ID}_Erest_ctrl_{winsize}_0p5' 
        
        # open a file, where you stored the pickled data
        file = open((synaptic_dir + f'/miniML_dtc-validation/' + filename + '.pickle'), 'rb')
        
        # dump information to that file
        detection = pickle.load(file)
        
        # close and remove (from memory) the file
        file.close()
        del file 
        gc.collect()
        
        # set column
        key = f'{cell_ID}_{str(winsize)}'
        
        ## trace
        traces[cell_ID] = detection['mini_trace']
        events[key] = detection['events']
        event_peak_locations[key] = detection['event_location_parameters']['event_peak_locations']
        
        # amplitude histogram
        amplitudes_hist.loc[:amplitudes_hist.index[-2], key], _ = np.histogram(a = detection['individual_values']['amplitudes'], bins = amplitudes_hist.index.to_list())
    



# %% figure one cell


def plot_trace(ax,
               trace,
               all_peak_idc,
               t = 0,
               trange = 60,
               ymin = -30,
               ymax = 5,
               xstep = 30,
               xstepminor = 10,
               xscale = 'ms',
               label = 'D'):


    # get plotting indices
    x = np.arange(0, trace.shape[0] / SR, 1/SR)
    plt_idc = np.arange((0+t) * SR, (trange+t) * SR, step = 1, dtype = int)
        
    # plot trace
    ax.plot(x[plt_idc], trace[plt_idc],
            color='k',
            lw=0.5,
            label='_nolegend_',
            zorder = 5)
    
    # zeroline
    ax.hlines(xmin = x[plt_idc][0], xmax = x[plt_idc][-1], y = 0, 
                    lw = 0.5, 
                    color = 'k', 
                    ls = 'dashed', 
                    alpha = 0.5,
                    zorder = 0)
    
    for wi, winsize in enumerate(winsizes):
    
        # get event peaks
        event_peak_idc = event_peak_locations[f'{cell_ID}_{str(winsize)}']
        
        # limit event peak indices
        plt_peak_idc = [int(idx) for idx in event_peak_idc if (idx > plt_idc[0] and idx < plt_idc[-1])]
        
        if len(plt_peak_idc) > 0:
            # plot event peak indicators
            ax.scatter(x = x[plt_peak_idc],
                       y = trace[plt_peak_idc],
                       color = winsize_colors[winsize],
                       marker = 'o',
                       s = 20-(wi*6),
                       lw = 1,
                       zorder = wi+1,
                       label = str(winsize) + ' ms')
    
    ypad = (ymax - ymin) * 0.01
    ax.set_ylabel('Current [pA]')
    ax.set_ylim([ymin - ypad, ymax + ypad])
    ax.spines['left'].set_bounds([ymin, ymax])
    ax.set_yticks(ticks = [ymin, 0, ymax], labels = [ymin, '', ymax])
    ax.set_yticks(ticks = np.arange(ymin, ymax+1, 5), minor = True)
    
    tpad = trange * 0.01
    ax.set_xlabel(f'Time [{xscale}]')
    ax.set_xlim([t - tpad, (t+trange) + tpad])
    ax.spines['bottom'].set_bounds([t, (t+trange)])
    
    if xscale == 's':
        ax.set_xticks(ticks = np.arange(t, (t+trange)+0.0001, xstep))
        ax.set_xticks(ticks = np.arange(t, (t+trange)+0.0001, xstepminor), minor = True)
    elif xscale == 'ms':
        ax.set_xticks(ticks = np.arange(t, (t+trange)+0.0001, xstep), 
                      labels = np.arange(0, (trange+0.0001)*1e3, xstep*1e3, dtype = int))
        ax.set_xticks(ticks = np.arange(t, (t+trange)+0.0001, xstepminor), minor = True)

    # add indicator to full trace
    axs[0].arrow(x = t+(trange/2), y = 10,
                 dx = 0, dy = -1,
                 head_length = 1,
                 head_width = 0.5,
                 length_includes_head = True,
                 facecolor = colors_dict['primecolor'])
    
    axs[0].text(x = t+(trange/2), y = 10.2,
                s = label,
                va = 'bottom', ha = 'center',
                fontsize = 6)
    
    # legend
    h, l = ax.get_legend_handles_labels()
    ax.legend(h, ['']*len(h), loc = 'lower right', frameon = False,
              borderaxespad = 0.0, borderpad = 0.0, labelspacing = 0.0,
              markerfirst = False)


# figure one cell


# init plotting
from functions.initialize_plotting import *


winsize_colors = {36: '#003f5c', 96 : '#7a5195', 114 : '#ef5675', 276 : '#ffa600'}

# winsize_colors = {36: '#656565', 96 : '#808782', 114 : '#a6d3a0', 276 : '#d1ffd7'}

cell_ID = 'E-303'
trace = traces[cell_ID]


fig, axs = plt.subplot_mosaic(mosaic = [[0,0,0,0],
                                        [1,1,2,3],
                                        [1,1,4,5],
                                        [6,7,8,9],
                                        [10,11,12,13]],
                              figsize = get_figure_size(width = 159.2, height = 180),
                              dpi = 300,
                              layout = 'constrained')

# # plot trace
ax = axs[0]
ax.set_title('A: Trace', loc = 'left')
t = 0
trange = 180
SR = 100_000
x = np.arange(0, trace.shape[0] / SR, 1/SR)

# get plotting indices
plt_idc = np.arange((0+t) * SR, (trange+t) * SR, step = 1, dtype = int)
    
# plot trace
ax.plot(x[plt_idc], trace[plt_idc],
        color='k',
        lw=0.5,
        label='_nolegend_',
        zorder = 5)

# zeroline
ax.hlines(xmin = x[plt_idc][0], xmax = x[plt_idc][-1], y = 0, 
                lw = 0.5, 
                color = 'k', 
                ls = 'dashed', 
                alpha = 0.5,
                zorder = 0)

for wi, winsize in enumerate(winsizes):

    # get event peaks
    event_peak_idc = event_peak_locations[f'{cell_ID}_{str(winsize)}']
    
    # limit event peak indices
    plt_peak_idc = [int(idx) for idx in event_peak_idc if (idx > plt_idc[0] and idx < plt_idc[-1])]
    
    # plot event peak indicators
    ax.scatter(x = x[plt_peak_idc],
               y = trace[plt_peak_idc],
               color = winsize_colors[winsize],
               marker = 'o',
               s = 8-(wi*2),
               lw = 1,
               zorder = wi+1,
               label = str(winsize) + ' ms')

ax.set_ylabel('Current [pA]')
ax.set_ylim([-30 - 0.4, 8 + 2.5])
ax.spines['left'].set_bounds([-30, 8])
ax.set_yticks(ticks = np.arange(-30, 8, 10))
ax.set_yticks(ticks = np.arange(-30, 8+0.1, 2), minor = True)

ax.set_xlabel('Time [s]')
ax.set_xlim([t - 1.8, (t+trange) + 1.8])
ax.spines['bottom'].set_bounds([t, (t+trange)])
ax.set_xticks(ticks = np.arange(t, (t+trange)+1, 30))
ax.set_xticks(ticks = np.arange(t, (t+trange)+1, 5), minor = True)

ax.legend(loc='lower right', frameon = False, ncols = 4, fontsize = 7,
          title = 'winsize', title_fontsize = 7, handletextpad = 0.1, 
          borderaxespad = 0.1, columnspacing = 0.3)

# plot histogram
ax = axs[1]
ax.set_title('B: Amplitude histogram', loc = 'left')

for winsize in winsizes:
    
    # get data
    values = amplitudes_hist[f'{cell_ID}_{str(winsize)}'].dropna().to_list()
    edges = amplitudes_hist.index.to_list() 
    
    ax.stairs(values, edges,
              lw = 1.5,
              color = winsize_colors[winsize],
              label = str(winsize) + ' ms')
    
ax.legend(loc = 'upper left', frameon = False, fontsize = 7,
          title = 'winsize', title_fontsize = 7)
    
# edit axis
ax.set_ylabel('Event count [#]')
ax.set_ylim([0-2.6, 130+2.6])
ax.set_yticks(ticks = np.arange(0, 130+1, 20))
ax.set_yticks(ticks = np.arange(0, 130+1, 5), minor = True)
ax.spines['left'].set_bounds([0, 130])

ax.set_xlabel('Amplitude [pA]')
ax.set_xticks(ticks = np.arange(-20, 0+1, 10),
              labels = np.arange(-20, 0+1, 10))
ax.set_xticks(ticks = np.arange(-25, 0, 1), minor = True)
ax.set_xlim([-25-.3, 0+.3])
ax.spines['bottom'].set_bounds([-25, 0])
        

# plot average events
for wi, winsize in enumerate(winsizes):

    ax = axs[wi+2]
    
    # get data
    key = f'{cell_ID}_{str(winsize)}'
    event_x = np.arange(0, events[key].shape[1]) * 1/(SR/1e3)
    event_average = np.mean(events[key], axis=0)

    # plot
    ax.plot(event_x, events[key].T,
            c = 'k',
            alpha = 0.1,
            lw = 0.5,
            label = '_nolegend_')

    ax.plot(event_x, event_average,
            c = winsize_colors[winsize], 
            lw = 1.0,
            label = 'average event')
    
    # axis
    if wi in [0, 2]:
        ax.set_ylabel('Current [pA]')
    ax.set_ylim([-30-0.3, 7+0.3])
    ax.set_yticks(ticks = np.arange(-30, 7, 10), labels = ['', -20, '', 0])
    ax.set_yticks(ticks = np.arange(-30, 7, 2), minor = True)
    ax.spines['left'].set_bounds([-30, 7])
    
    xmax = events[key].shape[1]/(SR/1e3)
    xpad = xmax*0.01
    xstepdict = {36: [50, 10], 96 : [50, 25], 114 : [100, 25], 276 : [200, 25]}
    xstep = xstepdict[winsize][0]
    xstepminor = xstepdict[winsize][1]
    
    ax.set_xlim([0-xpad, xmax+xpad])
    ax.spines['bottom'].set_bounds([0, xmax])
    ax.set_xticks(ticks = np.arange(0, xmax+1, xstep))
    ax.set_xticks(ticks = np.arange(0, xmax+1, xstepminor), minor = True)
    
    if wi in [2, 3]:
        ax.set_xlabel('Time [ms]')
    
axs[2].set_title('C: Average events', loc = 'left')


axs[6].set_title('D: Large events', loc = 'left')
plot_trace(ax = axs[6], trace=trace, all_peak_idc=event_peak_idc, t=8.732, trange=0.075, ymin=-20, ymax=5, xstep=0.05, xstepminor=0.025, xscale='ms', label='D')
axs[10].set_title('E:', loc = 'left')
plot_trace(ax = axs[10], trace=trace, all_peak_idc=event_peak_idc, t=26.591, trange=0.075, ymin=-20, ymax=5, xstep=0.05, xstepminor=0.025, xscale='ms', label='E')

axs[7].set_title('F: Medium events', loc = 'left')
plot_trace(ax = axs[7], trace=trace, all_peak_idc=event_peak_idc, t=35.518, trange=0.075, ymin=-15, ymax=5, xstep=0.05, xstepminor=0.025, xscale='ms', label='F')
axs[11].set_title('G:', loc = 'left')
plot_trace(ax = axs[11], trace=trace, all_peak_idc=event_peak_idc, t=58.40, trange=0.075, ymin=-15, ymax=5, xstep=0.05, xstepminor=0.025, xscale='ms', label='G')

axs[8].set_title('H: Small events', loc = 'left')
plot_trace(ax = axs[8], trace=trace, all_peak_idc=event_peak_idc, t=1.205, trange=0.075, ymin=-10, ymax=10, xstep=0.05, xstepminor=0.025, xscale='ms', label='H')
axs[12].set_title('I:', loc = 'left')
plot_trace(ax = axs[12], trace=trace, all_peak_idc=event_peak_idc, t=68.005, trange=0.075, ymin=-10, ymax=10, xstep=0.05, xstepminor=0.025, xscale='ms', label='I')

axs[9].set_title('J: Small events', loc = 'left')
plot_trace(ax = axs[9], trace=trace, all_peak_idc=event_peak_idc, t=105, trange=0.075, ymin=-10, ymax=10, xstep=0.05, xstepminor=0.025, xscale='ms', label='J')
axs[13].set_title('K:', loc = 'left')
plot_trace(ax = axs[13], trace=trace, all_peak_idc=event_peak_idc, t=140.85, trange=0.075, ymin=-10, ymax=10, xstep=0.05, xstepminor=0.025, xscale='ms', label='K')




# align labels
fig.align_labels()

# remove spines
[axs[axi].spines[spine].set_visible(False) for axi in axs for spine in ['top', 'right']]

# display figure
plt.show()


# import directories 
from parameters.directories_win import figure_dir
figure_path = figure_dir + '/miniML_validation'
save_figures(fig, 
             figure_name = f'miniML_validation-window_size_amplitudes-{cell_ID}', 
             save_dir = figure_path,
             figure_format = 'both')


# %% collect histogram data

ampl_hist_pop = pd.DataFrame(columns = winsizes, index = amplitudes_hist.index)

for winsize in winsizes:
    keys = list()
    for cell_ID in cell_IDs:
        keys.append(f'{cell_ID}_{str(winsize)}')
        
    ampl_hist_pop[winsize] = amplitudes_hist.loc[:, keys].mean(axis = 1)

# %%
# # plot histogram

fig, ax = plt.subplots(nrows = 1, ncols = 1,
                       figsize = get_figure_size(width = 100, height = 80),
                       dpi = 300,
                       layout = 'constrained')


for winsize in winsizes:

    # get data
    keys = list()
    for cell_ID in cell_IDs:
        keys.append(f'{cell_ID}_{str(winsize)}')
        
    means = amplitudes_hist.loc[:, keys].mean(axis = 1).dropna().to_list()
    stds = amplitudes_hist.loc[:, keys].std(axis = 1).dropna().to_list()
    edges = amplitudes_hist.index.to_list() 
    
    ax.stairs(means, edges,
              lw = 1.,
              color = winsize_colors[winsize],
              label = str(winsize) + ' ms')
    
    ax.stairs(values = [m+s for m, s in zip(means,stds)], 
              edges = edges,
              baseline = [m-s for m, s in zip(means,stds)],
              fill = True,
              color = winsize_colors[winsize],
              alpha = 0.25,
              label = '_nolegend_')
    
ax.legend(loc = 'upper left', frameon = False, fontsize = 7,
          title = 'winsize', title_fontsize = 7)
    
# edit axis
ax.set_ylabel('Event count [#]')
ax.set_ylim([0-2.6, 130+2.6])
ax.set_yticks(ticks = np.arange(0, 130+1, 20))
ax.set_yticks(ticks = np.arange(0, 130+1, 5), minor = True)
ax.spines['left'].set_bounds([0, 130])

ax.set_xlabel('Amplitude [pA]')
ax.set_xticks(ticks = np.arange(-50, 0+1, 10),
              labels = np.arange(-50, 0+1, 10))
ax.set_xticks(ticks = np.arange(-50, 0, 1), minor = True)
ax.set_xlim([-50-.3, 0+.3])
ax.spines['bottom'].set_bounds([-50, 0])

# align labels
fig.align_labels()

# remove spiness
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# display figure
plt.show()