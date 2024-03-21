# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 19:03:13 2024

@author: nesseler
"""

import os
os.chdir('C:/Users/nesseler/E-Phys')

import matplotlib.pyplot as plt
import numpy as np


from parameters.directories_win import table_file, raw_data_dir
from parameters.PGFs import cc_IF_parameters

import pandas as pd
from functions.functions_useful import calc_time_series, butter_filter, round_to_base

from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_plotting import get_colors, save_figures, get_figure_size, set_font_sizes
from functions.functions_constructors import construct_current_array

import scipy as sc
from parameters.parameters import min_peak_prominence_ccIF, min_peak_distance_ccIF, min_max_peak_width_ccIF

from matplotlib import animation
from matplotlib.animation import FuncAnimation 
from matplotlib.collections import LineCollection


# %%

table = pd.read_excel(table_file,
                      sheet_name="PGFs",
                      index_col='cell_ID')

# test cell E-092
cell_ID = 'E-122'

    
# PGF to load
PGF = 'cc_IF'
# PGF_parameters = cc_APs_parameters[frequency]




# %%


# get indices of current cell with the dataframe containing all indices    
# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)

# get IF data form file
i, v, t, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')

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

# find peaks
idc_peaks, dict_peak = sc.signal.find_peaks(vf, 
                                            prominence = min_peak_prominence_ccIF, 
                                            distance = min_peak_distance_ccIF * (SR_ms),
                                            width = np.multiply(min_max_peak_width_ccIF, SR_ms))



# %% create current array

# get hold current as table
i_hold = pd.read_excel(table_file, sheet_name="V_or_I_hold", index_col='cell_ID').at[cell_ID, PGF]

# calculate current steps relative to I_hold
## rounded to nearest 5
i_hold_rounded = round_to_base(i_hold, 5)

# get current arrays and list of input current relative to i_hold
i, i_input = construct_current_array(i_hold = i_hold_rounded,
                                     n_steps = n_steps,
                                     parameters_dict = cc_IF_parameters,
                                     SR_ms = SR_ms)




# %% split concatenate arrays back to steps wise 
# needs to occurs after filtering because of the filtering artifact

v = [None] * n_steps
t = [None] * n_steps
peaks = [None] * n_steps

step_dur = cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim'] + cc_IF_parameters['t_post']
step_points = step_dur * SR_ms

for idx in np.arange(0, n_steps, 1):
    start_idx = int(step_points * idx)
    stop_idx = int(start_idx + step_points - 1)
    
    v[idx] = vf[start_idx:stop_idx]
    t[idx] = t_ms[start_idx:stop_idx]
    peaks[idx] = [(idx_peak / SR_ms) - (step_dur * idx) for idx_peak in idc_peaks if idx_peak > start_idx and idx_peak < stop_idx]


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
                      height_ratios = [1,4]) 


ax[0].plot(t[0], i[9], 'm')

segs = return_segments(t[0], v[:9])

line_segments = LineCollection(segs, colors = 'grey')

ax[1].add_collection(line_segments)

ax[1].plot(t[0], v[9])
# ax[1].eventplot(t_peaks_s, color = 'r', lineoffsets=30, linelengths=5)
ax[1].set_ylim([-100, 50])

plt.show()


# %% where is the threshold

# idx_u_th = next(p for p, idx_peak in enumerate(peaks) if len(idx_peak) > 0)

# i_u_th = i_steps[idx_u_th]


# %% animation

# see link
# https://spikesandbursts.wordpress.com/2024/01/04/patch-clamp-data-analysis-animate-time-series/

# Initialise figure
fig_ani, ax_ani = plt.subplots(2, 1, 
                               figsize = get_figure_size(width=245.252),
                               sharex = 'col',
                               height_ratios = [1,4])

# set font sizes
set_font_sizes()

# plot parameters
## initialise a line plot
linesegments_v = LineCollection([], colors = 'grey', lw = 2)
linecollection = ax_ani[1].add_collection(linesegments_v)

line_v, = ax_ani[1].plot([],[], lw = 2, color = color_dict['primecolor'])

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
ax_ani[0].set_ylim([-100, 300])
ax_ani[0].set_yticks(ticks = np.arange(-50, 200+1, 50), 
                     labels = [None, 0, None, 100, None, 200])
ax_ani[0].set_yticks(np.arange(-50, 200+1, 10), minor = True)

v_range = [-150, 75]

ax_ani[1].set_xlabel('Time [ms]')
ax_ani[1].set_xlim([0, 1500])

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
    
    # if frame >= idx_u_th:
    #     text.set_text(f'Threshold between\n{i_u_th - 10} pA and {i_u_th} pA!')
    # else:
    #     text.set_text('')
    
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
duration_s = n_steps * 0.5  # Select the length of the video here
 
anim_fps = frames/duration_s
 
extra_args = [
    '-vcodec', 'libx264',  # Video codec
    ]  # Set the output resolution, same ratio as figure size

anim_fps_str = str(round(anim_fps, 2))


#%% MP4

# Create the video with FFMpegWriter
writer = animation.FFMpegWriter(fps=anim_fps, bitrate=3000, 
                                codec="h264",  extra_args=extra_args)


 
# Save path
anim.save(f'C:/Users/nesseler/Desktop/local E-Phys/figures/{cell_ID}_{PGF}_animation-{anim_fps_str}.mp4', writer = writer)

# anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/video.mp4', writer=writer)
 
# Print the fps and duration of the final video
print('Video_duration (s):', frames/anim_fps)
print('Video_fps:', anim_fps)


# %% GIF


# writergif = animation.PillowWriter(fps=anim_fps, bitrate=2000)
 
# anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/1APth_video.gif', writer=writergif)

# print('Done!')


























