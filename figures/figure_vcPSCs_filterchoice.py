# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 13:35:43 2025

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

# # filter
# filter all data with 1kHz cutoff
order = 3
cutoff = 1000 #Hz

i_butter = butter_filter(i, order=order, cutoff=cutoff, sampling_rate=SR)
i_bessel = bessel_filter(i, order=order, cutoff=cutoff, sampling_rate=SR)
i_cheby = cheby_filter(i, order=order, rp = 2, cutoff=cutoff, sampling_rate=SR)

# # replace first values with nans to eliminate filter artifact
# i[:100] = np.nan

# # add nan values if protocol was stopped before complete recording
# if n_steps != nsteps_loaded:
#     # expected datapoints
#     n_datapoints = vc_Erest_parameters['dur_steps'] * n_steps * SR
    
#     # generate filler array
#     filler = np.full(shape = (n_datapoints - len(i)), 
#                      fill_value = np.nan)
    
#     # add filler
#     i = np.concatenate([i, filler])

# # split back into steps
# i = np.array_split(i, n_steps)

# # initialize i
# i_new = np.full(shape = (n_steps, vc_Erest_parameters['dur_steps'] * SR), 
#                 fill_value = np.nan)

# # include only listed steps
# for i_step in range(n_steps):
    
#     if i_step in idc_steps:
#         i_new[i_step] = i[i_step]

# # concatenate steps
# i = np.array(i_new).flatten('C')

# # write to dict
# i_s[condition] = i

# %% plotting

# single event index
start = 4.365 #s
winsize = 0.025 #s
event_idc = np.arange(start*SR, (start+winsize)*SR, dtype = int)

# init plotting
from functions.initialize_plotting import *

# line plot dict
plt_dict = {'lw' : 0.5}

# x
xdict = {'ax_min' : 0,
         'ax_max' : 30,
         'pad' : None,
         'step' : 10,
         'stepminor' : 5,
         'label' : ''}

xdict_zoom = {'ax_min' : start,
              'ax_max' : start + winsize,
              'pad' : None,
              'step' : winsize,
              'stepminor' : winsize,
              'label' : ''}

# y
ydict = {'ax_min' : -40,
         'ax_max' : 15,
         'pad' : None,
         'step' : 40,
         'stepminor' : 10,
         'label' : ''}

ydict_zoom = {'ax_min' : -30,
              'ax_max' : 0,
              'pad' : None,
              'step' : 30,
              'stepminor' : 2,
              'label' : ''}

# create rectangle
def create_rect():
    return Rectangle(xy = (start, -35), 
                           width = winsize, 
                           height = 40,
                           fill = False,
                           color = colors_dict['color3'],
                           linestyle = '--',
                           lw = 0.5)

# initialize figure
fig, axs = plt.subplots(nrows = 4,
                        ncols = 2,
                        figsize = get_figure_size(width = 159.2, height = 160),
                        width_ratios = [2, 1.5],
                        layout = 'constrained',
                        dpi = 600,
                        sharey = 'col',
                        sharex = 'col')

# set title
fig.suptitle(f'Filter choice - low pass\ncutoff_freq: {cutoff} Hz - filter_order: {order}',
             fontsize = 12)

axs = axs.flatten()

# original trace
# full trace
ax = axs[0]

ax.plot(t, i,
        color = colors_dict['primecolor'],
        **plt_dict)

# add rectangle marker
ax.add_patch(create_rect())

# single event
ax = axs[1]

ax.plot(t[event_idc], i[event_idc],
        color = colors_dict['primecolor'],
        label = 'unfiltered',
        **plt_dict)

ax.legend(loc = 'lower right',
          frameon = False,
          fontsize = 9,
          handletextpad = 0.4)


def add_traces(ax_i, t, i, i_filtered, zoom_idc, label, plt_dict):
    axs[ax_i].plot(t, i, alpha = 0.5, color = 'gray', **plt_dict)
    axs[ax_i].plot(t, i_filtered, color = colors_dict['primecolor'], **plt_dict)
    axs[ax_i].add_patch(create_rect())
    axs[ax_i+1].plot(t[event_idc], i[event_idc], color = 'gray', alpha = 0.5, label = 'unfiltered', **plt_dict)
    axs[ax_i+1].plot(t[event_idc], i_filtered[event_idc], color = colors_dict['primecolor'], label = label, **plt_dict)
    axs[ax_i+1].legend(loc = 'lower right', frameon = False, fontsize = 9, handletextpad = 0.4)
    
    
# butter filtered trace
add_traces(2, t, i, i_butter, event_idc, 'butterworth', plt_dict)

# bessel filtered trace
add_traces(4, t, i, i_bessel, event_idc, 'bessel', plt_dict)

# chebyshev filtered trace
add_traces(6, t, i, i_cheby, event_idc, 'chebyshev', plt_dict)

# set sup labels
fig.supylabel('Current [pA]', fontsize = 9)
axs[-2].set_xlabel('Time [s]')
axs[-1].set_xlabel('Time [ms]')

# axis titles
axis_titles = [r'$\mathregular{A_{i}}$: Original trace',
               r'$\mathregular{A_{ii}}$',
               r'$\mathregular{B_{i}}$: Butterworth filter',
               r'$\mathregular{B_{ii}}$',
               r'$\mathregular{C_{i}}$: Bessel filter',
               r'$\mathregular{C_{ii}}$',
               r'$\mathregular{D_{i}}$: Chebyshev filter',
               r'$\mathregular{D_{ii}}$']

for ax_i, ax in enumerate(axs):
    ax.set_title(axis_titles[ax_i],
                 loc = 'left',
                 fontsize = 9)


for ax in axs[::2]:
    apply_axis_settings(ax, axis = 'x', **xdict)
    apply_axis_settings(ax, axis = 'y', **ydict)
    
for ax in axs[1::2]:
    apply_axis_settings(ax, axis = 'x', **xdict_zoom)
    apply_axis_settings(ax, axis = 'y', **ydict_zoom)
    
    ax.set_xticks(ticks = np.arange(start, start+winsize+0.001, 0.005),
                  labels = np.arange(0, (winsize+0.001)*1e3, 5, dtype = int))

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# align labels
fig.align_labels()

# create saving path and save
from parameters.directories_win import figure_dir
path_fig = join(figure_dir, 'filter_choice')
save_figures(fig, f'filter_type', path_fig, darkmode_bool, figure_format='both')
