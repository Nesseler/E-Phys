# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 19:03:13 2024

@author: nesseler
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 10:36:47 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import numpy as np


from parameters.directories_win import table_file, raw_data_dir
from parameters.PGFs import cc_th1Ap_parameters

import pandas as pd
from functions.functions_useful import calc_time_series, butter_filter
import os

from functions.functions_ccIF import get_IF_data
from functions.functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes

import scipy as sc
from parameters.parameters import min_peak_distance, min_peak_prominence

from matplotlib import animation
from matplotlib.animation import FuncAnimation 
from matplotlib.collections import LineCollection


# %%

table = pd.read_excel(table_file,
                      sheet_name="PGFs",
                      index_col='cell_ID')


# frequencies = ['1Hz', '30Hz', '75Hz']

# # loop to create string to include all frequencies in query
# query_str = ''

# for idx, frequency in enumerate(frequencies):
#     PGF = 'cc_APs_' + frequency
    
#     if idx > 0:
#         query_str = query_str + ' and '
        
#     query_str = query_str + f'{PGF}.notnull()'
    

# limit lookup table
lookup_table = table

# %%

# test cell E-092
cell_ID = 'E-092'

    
# PGF to load
PGF = 'cc_th1AP'
# PGF_parameters = cc_APs_parameters[frequency]

# lookup_table = table.query(f'{PGF}.notnull()')


# get indices of current cell with the dataframe containing all indices    
group_idx = int(lookup_table.at[cell_ID, 'group'])-1
series_idx = int(lookup_table.at[cell_ID, f'{PGF}'])-1

# construct traceIndex with indices
traceIndex = [group_idx, series_idx, 0, 0]

# call on data file with indices from dataframe above
current_file = lookup_table.at[cell_ID, 'file']

data_file_path = os.path.join(raw_data_dir, current_file + '.dat')

data_file_path_str = fr"{data_file_path}"

# get IF data form file
i, v, t, SR, n_steps = get_IF_data(data_file_path_str, traceIndex, 'ms')

# sampling rate in ms
SR_ms = SR / 1e3

# concatenate individual steps
n_points = int(np.shape(i)[0] * np.shape(i)[1])

i_concat = i.flatten() #.reshape([n_points], order='F')

v_concat = v.flatten() #.reshape([n_points], order='F')

t_ms = calc_time_series(v_concat, SR)
t_s = calc_time_series(v_concat, SR, scale = 's')

# plt.plot(v_concat[0:15000])
# plt.show()


# filter voltage (to vf)
vf = butter_filter(v_concat, 
                   order = 3,
                   cutoff = 1e3,
                   sampling_rate = SR)


# %% find peaks

idx_peaks, dict_peak = sc.signal.find_peaks(vf, 
                                            prominence = min_peak_prominence, 
                                            distance = min_peak_distance * (SR_ms))



# %% create current array

i_hold = 0 # pA

i_start = cc_th1Ap_parameters['i_start']
i_delta = cc_th1Ap_parameters['i_delta']

i_steps= np.arange(i_start, i_start + (i_delta * n_steps), i_delta)

i = [None] * n_steps

for idx, i_stim in enumerate(i_steps):
    i_pre = np.full(int((SR_ms * cc_th1Ap_parameters['t_pre'])), i_hold)
    i_stim = np.full(int((SR_ms * cc_th1Ap_parameters['t_stim'])), i_stim)
    i_post = np.full(int((SR_ms * cc_th1Ap_parameters['t_post'])-1), i_hold)
    
    i_step = np.concatenate((i_pre, i_stim, i_post))
    i[idx] = i_step



# %% split concatenate arrays back to steps wise 
# needs to occurs after filtering because of the filtering artifact

v = [None] * n_steps
t = [None] * n_steps
peaks = [None] * n_steps

step_dur = cc_th1Ap_parameters['t_pre'] + cc_th1Ap_parameters['t_stim'] + cc_th1Ap_parameters['t_post']
step_points = step_dur * SR_ms

for idx in np.arange(0, n_steps, 1):
    start_idx = int(step_points * idx)
    stop_idx = int(start_idx + step_points - 1)
    
    v[idx] = vf[start_idx:stop_idx]
    t[idx] = t_ms[start_idx:stop_idx]
    peaks[idx] = [(idx_peak / SR_ms) - (step_dur * idx) for idx_peak in idx_peaks if idx_peak > start_idx and idx_peak < stop_idx]


# %%

def return_segments(x, ys):
    '''
    Function returns segments from single x array and multiple y array to use
    with matplotlib.collections.LineCollection.
    Parameters:
        x : Single array of common x-coordinates.
        ys : 2D array / list of y-coordinates to be plotted on common x-coordinates.
    Returns:
        segs : Segmented lines ready to use with 'LineCollection(segs)'
    '''
    n_ys = len(ys)
    n_x = len(x)
    segs = np.zeros((n_ys, n_x, 2))
    segs[:, :, 1] = ys
    segs[:, :, 0] = x
    return segs


# %%

darkmode_bool = True

color_dict, _ = get_colors(darkmode_bool)

fig, ax = plt.subplots(2,1,
                       gridspec_kw={'height_ratios': [1,4]}) 


ax[0].plot(t[0], i[9], 'm')

segs = return_segments(t[0], v[:9])

line_segments = LineCollection(segs, colors = 'grey')

ax[1].add_collection(line_segments)

ax[1].plot(t[0], v[9])
# ax[1].eventplot(t_peaks_s, color = 'r', lineoffsets=30, linelengths=5)
ax[1].set_ylim([-100, 50])

plt.show()


# %% where is the threshold

idx_u_th = next(p for p, idx_peak in enumerate(peaks) if len(idx_peak) > 0)

i_u_th = i_steps[idx_u_th]


# %% animation

# see link
# https://spikesandbursts.wordpress.com/2024/01/04/patch-clamp-data-analysis-animate-time-series/

# Initialise figure
fig_ani, ax_ani = plt.subplots(2, 1, 
                               figsize = get_figure_size(width=246.502),
                               sharex = 'col',
                               gridspec_kw={'height_ratios': [1,4]})

# set font sizes
set_font_sizes()

# plot parameters
## initialise a line plot
linesegments_v = LineCollection([], colors = 'grey', lw = 2)
linecollection = ax_ani[1].add_collection(linesegments_v)

line_v, = ax_ani[1].plot([],[], lw = 2, color = color_dict['color1'])

linesegments_i = LineCollection([], colors = 'grey', lw = 2)
linecollection = ax_ani[0].add_collection(linesegments_i)

line_i, = ax_ani[0].plot([],[], lw = 2, color = color_dict['color2'])

## initialise eventmarker
events, = ax_ani[1].eventplot([], color = 'r', lineoffsets=50, linelengths=10, linewidth = 3)

## initialise text for threshold label
text = ax_ani[0].text(x = 375,
                      y = 100,
                      s = '', 
                      fontsize = 14,
                      verticalalignment='center')

# axis settings
ax_ani[0].set_ylabel('Current [pA]')
ax_ani[0].set_ylim([-50, 200])
ax_ani[0].set_yticks(ticks = np.arange(-50, 200+1, 50), 
                     labels = [None, 0, None, 100, None, 200])
ax_ani[0].set_yticks(np.arange(-50, 200+1, 10), minor = True)

v_range = [-100, 75]

ax_ani[1].set_xlabel('Time [ms]')
ax_ani[1].set_xlim([0, 500])

ax_ani[1].set_ylabel('Voltage [mV]')
ax_ani[1].set_ylim(v_range)

ax_ani[1].set_yticks(np.arange(v_range[0], v_range[1] + 1, 50))
ax_ani[1].set_yticks(np.arange(v_range[0], v_range[1] + 1, 25), minor = True)

[ax.grid(False) for ax in ax_ani]

# animation function 
interval = 1   # in datapoints

def animate(frame):
    # end frame for each update in the animation
    # end_frame = (frame + 1) * interval
    print(frame)
    
    if frame > 0:
        segs_v = return_segments(t[0], v[:frame])
        linesegments_v.set_segments(segs_v)
        
        segs_i = return_segments(t[0], i[:frame])
        linesegments_i.set_segments(segs_i)
    
    # update line plot data
    line_v.set_data(t[0], v[frame])
    line_i.set_data(t[0], i[frame])
    
    # eventplot
    # end_frame_s = end_frame / SR
    events.set_positions(peaks[frame])
    
    if frame >= idx_u_th:
        text.set_text(f'Threshold between\n{i_u_th - 10} pA and {i_u_th} pA!')
    else:
        text.set_text('')
    
    return line_v, #events

# create the animation
frames = len(t) // interval

anim = animation.FuncAnimation(fig_ani, animate,
                                frames = frames,
                                interval = 5,
                                blit = True, 
                                repeat = False)



print(f'Animation frames: {frames}')

fig_ani.tight_layout()
plt.show()

# %% animation settings

# Animation parameters
duration_s = 10  # Select the length of the video here
 
anim_fps = frames/duration_s
 
extra_args = [
    '-vcodec', 'libx264',  # Video codec
    ]  # Set the output resolution, same ratio as figure size
 

#%% MP4

# Create the video with FFMpegWriter
writer = animation.FFMpegWriter(fps=anim_fps, bitrate=3000, 
                                codec="h264",  extra_args=extra_args)
 
# Save path
anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/1APth_video.mp4', writer = writer)

# anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/video.mp4', writer=writer)
 
# Print the fps and duration of the final video
print('Video_duration (s):', frames/anim_fps)
print('Video_fps:', anim_fps)


# %% GIF


writergif = animation.PillowWriter(fps=anim_fps, bitrate=2000)
 
anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/1APth_video.gif', writer=writergif)

print('Done!')


























