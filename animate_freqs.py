# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 10:36:47 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import directories_win as directories
from PGFs import cc_APs_parameters
import pandas as pd
from useful_functions import calc_time_series, butter_filter
import os
from cc_IF_functions import get_IF_data
from plotting_functions import get_colors, save_figures, get_figure_size, set_font_sizes
import scipy as sc
import parameters

from matplotlib import animation
from matplotlib.animation import FuncAnimation 


# %%

table = pd.read_excel(directories.table_dir + 'InVitro_Database.xlsx',
                      sheet_name="PGFs",
                      index_col='cell_ID')


frequencies = ['1Hz', '30Hz', '75Hz']

# loop to create string to include all frequencies in query
query_str = ''

for idx, frequency in enumerate(frequencies):
    PGF = 'cc_APs_' + frequency
    
    if idx > 0:
        query_str = query_str + ' and '
        
    query_str = query_str + f'{PGF}.notnull()'
    

# limit lookup table
lookup_table = table.query(query_str)

# %%

# test cell E-092
cell_ID = 'E-092'
frequency = '10Hz'



    
# PGF to load
PGF = 'cc_APs_' + frequency
PGF_parameters = cc_APs_parameters[frequency]

# lookup_table = table.query(f'{PGF}.notnull()')


# get indices of current cell with the dataframe containing all indices    
group_idx = int(lookup_table.at[cell_ID, 'group'])-1
series_idx = int(lookup_table.at[cell_ID, f'{PGF}'])-1

# construct traceIndex with indices
traceIndex = [group_idx, series_idx, 0, 0]

# call on data file with indices from dataframe above
current_file = lookup_table.at[cell_ID, 'file']

data_file_path = os.path.join(directories.raw_data_dir, current_file + '.dat')

data_file_path_str = fr"{data_file_path}"

# get IF data form file
i, v, t, SR, n_steps = get_IF_data(data_file_path_str, traceIndex, 'ms')

# sampling rate in ms
SR_ms = SR / 1e3

# concatenate individual steps
n_points = int(np.shape(i)[0] * np.shape(i)[1])

i_concat = i.flatten() #.reshape([n_points], order='F')

v_concat = v.flatten() #.reshape([n_points], order='F')

t = calc_time_series(v_concat, SR)
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
                                            prominence = parameters.min_peak_prominence, 
                                            distance = parameters.min_peak_distance * (SR_ms))

t_peaks_s = np.divide(idx_peaks, SR)

# %%

darkmode_bool = True

color_dict = get_colors(darkmode_bool)

plt.plot(t_s, vf)
plt.eventplot(t_peaks_s, color = 'r', lineoffsets=30, linelengths=5)
plt.ylim([-100, 50])
plt.show() 


# %% animation

# see link
# https://spikesandbursts.wordpress.com/2024/01/04/patch-clamp-data-analysis-animate-time-series/

# Initialise figure
fig_ani, ax_ani = plt.subplots(figsize = get_figure_size())

# set font sizes
set_font_sizes()

# plot parameters
## initialise a line plot
line, = ax_ani.plot([],[], lw = 1, color = color_dict['color1'])

## initialise eventmarker
events, = ax_ani.eventplot([], color = 'r', lineoffsets=40, linelengths=5)

# axis settings
ax_ani.set_xlabel('Time [s]')
ax_ani.set_xlim([0, 10])

ax_ani.set_ylabel('Voltage [mV]')
ax_ani.set_ylim([-100, 50])


# animation function 
interval = 10000   # in datapoints

def animate(frame):
    # end frame for each update in the animation
    end_frame = (frame + 1) * interval
    
    # update line plot data
    line.set_data(t_s[:end_frame], vf[:end_frame])
    
    # eventplot
    end_frame_s = end_frame / SR
    events.set_positions(t_peaks_s[t_peaks_s < end_frame_s])
    
    return line, events

# create the animation
frames = len(t_s) // interval

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
anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/video.mp4', writer = writer)

# anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/video.mp4', writer=writer)
 
# Print the fps and duration of the final video
print('Video_duration (s):', frames/anim_fps)
print('Video_fps:', anim_fps)


# %% GIF


writergif = animation.PillowWriter(fps=anim_fps, bitrate=2000)
 
anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/video.gif', writer=writergif)

print('Done!')


























