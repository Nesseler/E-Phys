# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 16:00:28 2025

@author: nesseler
"""


# import standard packages
from functions.initialize_packages import *

# import functions
from functions.functions_import import get_traceIndex_n_file, get_vc_data, get_PSCs_steps
from functions.functions_filter import butter_filter, bessel_filter, cheby_filter

# import parameters
from parameters.PGFs import vc_Erest_parameters

SR = vc_Erest_parameters['SR']
t = vc_Erest_parameters['t']
n_steps = vc_Erest_parameters['n_steps']
PGF = 'vc_Erest_3min'


# %% load data

cell_ID = 'E-301'
condition = 'ctrl'


# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')

# get list of steps to include
idc_steps = get_PSCs_steps(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')

# get data with file path & trace index
i, _, _, _, nsteps_loaded = get_vc_data(file_path, traceIndex, scale='ms')

# concatenate steps
i = i.flatten('C')

i = i[:3000000]
t = t[:3000000]

# # filter
# filter all data with 1kHz cutoff
order = 3
# cutoff = 1000 #Hz

# set dict
i_freqdict = dict.fromkeys([2000, 1000, 750, 500, 250])

# filter
for freq in i_freqdict.keys():
    i_freqdict[freq] = bessel_filter(i, order=order, cutoff=freq, sampling_rate=SR)


# %% plotting

# single event index
start = 4.365 #s
winsize = 0.020 #s
event_idc = np.arange(start*SR, (start+winsize)*SR, dtype = int)
# event_idc = np.arange(0, 30*SR, dtype = int)

# init plotting
from functions.initialize_plotting import *

# specify color map
cmap_str = 'cividis_r'

# min max normalize time for color-code
norm = mtl.colors.Normalize(0, 2000)

# create mappable colormap object for colorbar
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

fig = plt.figure(dpi = 600,
                 figsize = get_figure_size(width = 130, height = 100))

plt.plot(t[event_idc], i[event_idc], 
         lw = 0.75,
         label = 'unfiltered',
         color = 'k')

plt.text(start, y = 0, s = 'unfiltered', color = colors_dict['primecolor'], va = 'top')

for freq_i, freq in enumerate(i_freqdict.keys()):
    plt.plot(t[event_idc], i_freqdict[freq][event_idc] + ((freq_i+1)*10), 
             lw = 0.75,
             label = freq,
             color = cmap.to_rgba(freq))
    
    plt.text(start, y = ((freq_i+1)*10), s = f'{freq} Hz', color = cmap.to_rgba(freq), va = 'top')
    
# plt.legend(loc = 'lower right',
#            frameon = False,
#            fontsize = 9,
#            handletextpad = 0.4) 

plt.vlines(x = [4.3696, 4.3705, 4.3751], ymin = -30, ymax = 50, lw = 0.5, color = 'k', alpha = 0.5)
# plt.vlines(x = 4.3705, ymin = -30, ymax = 50, lw = 0.5, color = 'k', alpha = 0.5)

ydict = {'ax_min' : -10,
         'ax_max' : 0,
         'pad' : None,
         'step' : 10,
         'stepminor' : 10,
         'label' : ''}

plt.gca().spines['left'].set_bounds([-10, 0])
plt.gca().set_ylim([-30, 50])  
plt.gca().set_yticks(ticks = np.arange(-10, 0.1, 10),
                     labels = np.arange(-10, 0.1, 10))
plt.gca().set_ylabel('Current [pA]')
plt.gca().yaxis.set_label_coords(-0.098, 0.31)

# remove spines
[plt.gca().spines[spine].set_visible(False) for spine in ['top', 'right']]

# x axis
plt.gca().set_xlabel('Time [ms]')
plt.gca().set_xticks(ticks = np.arange(start, start+winsize+0.001, 0.005),
                     labels = np.arange(0, (winsize+0.001)*1e3, 5, dtype = int))
plt.gca().spines['bottom'].set_bounds([start, start+winsize])
plt.gca().set_xlim([start-(winsize/100), start+winsize+(winsize/100)])
    
plt.show()

# create saving path and save
from parameters.directories_win import figure_dir
path_fig = join(figure_dir, 'filter_choice')
save_figures(fig, f'filter_frequency-single_event', path_fig, darkmode_bool, figure_format='both')


# %% full trace figure

lw = 0.75

# initialize figure
fig, axs = plt.subplot_mosaic(mosaic = 'aa;bc;de',
                              figsize = get_figure_size(width = 159.2, height = 225.154),
                              layout = 'constrained',
                              height_ratios = [2, 1, 1],
                              dpi = 600)

# set title
fig.suptitle(f'Filter frequency choice - Bessel low pass filter',
             fontsize = 12)


ax = axs['a']

ax.set_title('A: Full trace', loc = 'left', fontsize = 9)

ax.plot(t, i, lw = lw, label = 'unfiltered', color = colors_dict['primecolor'])
ax.text(x = 0.1, y = 2, s = 'unfiltered', color = colors_dict['primecolor'], va = 'bottom')

for freq_i, freq in enumerate(i_freqdict.keys()):
    ax.plot(t, i_freqdict[freq] + ((freq_i+1)*30), 
             lw = 0.75,
             label = freq,
             color = cmap.to_rgba(freq))
    
    ax.text(x = 0.1, y = (5+(freq_i+1)*29.5), s = f'{freq} Hz', color = cmap.to_rgba(freq), va = 'top')

# x
ax.spines['bottom'].set_bounds([0, 5])
ax.set_xlim([-0.3, 30.3])  
ax.set_xticks(ticks = np.arange(0, 5.1, 5))
ax.set_xlabel('Time [s]')
ax.xaxis.set_label_coords(0.09, -0.06)

# y
ax.spines['left'].set_bounds([-30, 0])
ax.set_ylim([-35, 155])  
ax.set_yticks(ticks = np.arange(-30, 0.1, 30))
ax.set_ylabel('Current [pA]')
ax.yaxis.set_label_coords(-0.07, 0.105)

# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]


# single events

def create_single_event_plot(axs_l, start = 4.365, winsize = 0.020, title = ''):
    
    # single event index
    # start = 4.365 #s
    # winsize = 0.020 #s
    event_idc = np.arange(start*SR, (start+winsize)*SR, dtype = int)
    
    # add indicator to panel a
    axs['a'].arrow(x = start, y = -31,
                   dx = 0, dy = 1,
                   head_length = 1,
                   head_width = 0.1,
                   length_includes_head = True)
    axs['a'].text(x = start+0.3, y = -30.5,
                  s = axs_l,
                  va = 'center', ha = 'center',
                  fontsize = 7)
    
    # set axis
    ax = axs[axs_l]
    
    # axis title
    ax.set_title(axs_l.capitalize() + ': ' + title, loc = 'left', fontsize = 9)

    ax.plot(t[event_idc], i[event_idc], 
             lw = 0.75,
             label = 'unfiltered',
             color = 'k')

    ax.text(x = start, y = 1.5, s = 'unfiltered', color = colors_dict['primecolor'], va = 'top')
    
    for freq_i, freq in enumerate(i_freqdict.keys()):
        ax.plot(t[event_idc], i_freqdict[freq][event_idc] + ((freq_i+1)*10), 
                 lw = 0.75,
                 label = freq,
                 color = cmap.to_rgba(freq))
        
        ax.text(x = start, y = (1+(freq_i+1)*10), s = f'{freq} Hz', color = cmap.to_rgba(freq), va = 'top')
        
    # x
    ax.spines['bottom'].set_bounds([start, start+winsize])
    ax.set_xlim([start-(winsize/100), start+winsize+(winsize/100)])  
    ax.set_xticks(ticks = np.arange(start, start+winsize+0.001, 0.005),
                  labels = np.arange(0, (winsize+0.001)*1e3, 5, dtype = int))
    ax.set_xlabel('Time [ms]')

    # y
    ax.spines['left'].set_bounds([-30, 0])
    ax.set_ylim([-35, 50])  
    ax.set_yticks(ticks = np.arange(-30, 0.1, 30))
    ax.set_ylabel('Current [pA]')
    ax.yaxis.set_label_coords(-0.15, 0.24)
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
        
        
create_single_event_plot(axs_l = 'b', start = 4.360, winsize = 0.030, title = 'Fast event')
create_single_event_plot(axs_l = 'c', start = 7.513, winsize = 0.030, title = 'Background') #10.253
create_single_event_plot(axs_l = 'd', start = 10.885, winsize = 0.030, title = 'Slow event') #12.723
create_single_event_plot(axs_l = 'e', start = 17.869, winsize = 0.030, title = 'Noise')

# align labels
fig.align_labels()

# create saving path and save
from parameters.directories_win import figure_dir
path_fig = join(figure_dir, 'filter_choice')
save_figures(fig, f'filter_frequency', path_fig, darkmode_bool, figure_format='both')
