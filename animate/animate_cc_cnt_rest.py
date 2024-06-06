# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 15:35:53 2024

@author: nesseler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from os.path import join, exists
from os import mkdir

from parameters.directories_win import table_dir, quant_data_dir, vplot_dir, table_file, figure_dir
from parameters.PGFs import cc_cntrest_parameters


from functions.functions_import import get_traceIndex_n_file
from functions.functions_ccIF import get_IF_data
from functions.functions_plotting import get_colors, save_figures, plot_t_vs_v, get_figure_size, set_font_sizes
from functions.functions_useful import calc_time_series, butter_filter


from matplotlib import animation
from matplotlib.animation import FuncAnimation 


table = pd.read_excel(table_file,
                      sheet_name="PGFs",
                      index_col='cell_ID')

# loop to create string to include all frequencies in query
PGF = 'cc_rest'  

# limit lookup table
lookup_table = table.query(f'{PGF}.notnull()')

# cell IDs 
cell_IDs = lookup_table.query('cc_cnt_rest.notnull()').index.to_list()

autocorr_dir = join(quant_data_dir, 'cnt_rest', 'autocorrelations')

# autocorrelation dataframe will be filled for all cells
autocorr_df = pd.DataFrame(columns=cell_IDs)


# plotting specifications
darkmode_bool = True

vplot_bool = True

save_bool = True

autocorr_bool = False

colors_dict, _ = get_colors(darkmode_bool)

cell_ID = 'E-176'


# vplot_dir_cell = join(vplot_dir, PGF, cell_ID)
# if not exists(vplot_dir_cell):
#     mkdir(vplot_dir_cell)

# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)

# get IF data form file
_, v, _, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')

# sampling rate in ms
SR_ms = SR / 1e3

# concatenate individual steps
n_points = int(np.shape(v)[0] * np.shape(v)[1])

# reconstruct i array from parameters
i = np.multiply(np.ones(n_points), cc_cntrest_parameters['i_hold'])

v_concat = v.flatten() 

print(f'Started: {cell_ID}')

t_ms = calc_time_series(v_concat, SR)
t_s = calc_time_series(v_concat, SR, scale = 's')
t_total = len(t_s) / SR


# filter voltage (to vf)
vf = butter_filter(v_concat, 
                    order = 3,
                    cutoff = 1e3,
                    sampling_rate = SR)


# plot_t_vs_v(t_s[:5000000], vf[:5000000])

t_total = 30 #s

points_total = int(t_total * SR)

t = t_s[:points_total]
v = vf[:points_total]
i = np.zeros(len(t))


# %% animation

# see link
# https://spikesandbursts.wordpress.com/2024/01/04/patch-clamp-data-analysis-animate-time-series/

# darkmode_bool = True

# colors_dict, _ = get_colors(darkmode_bool)

# Initialise figure
fig_ani, ax_ani = plt.subplots(2, 1, 
                               figsize = get_figure_size(width=80.917, height=56.209),
                               sharex = 'col',
                               height_ratios = [1, 5])

# set font sizes
set_font_sizes()

# plot parameters
line_v, = ax_ani[1].plot([], [], color=colors_dict['primecolor'], lw=1)
line_i, = ax_ani[0].plot([], [], color='r', lw=1)

# axis settings
# ax_ani[0].set_ylabel('Current [pA]')
ax_ani[0].set_ylabel('')
ax_ani[0].set_ylim([-50, 100])
ax_ani[0].set_yticks(ticks = np.arange(-50, 100+1, 50), 
                      labels = [None, 0, None, None])
ax_ani[0].set_yticks(np.arange(-50, 100+1, 25), minor = True)

v_range = [-100, 40]

# ax_ani[1].set_xlabel('Time [s]')
ax_ani[1].set_xlabel('')
ax_ani[1].set_xlim([0, t_total])

# ax_ani[1].set_ylabel('Voltage [mV]')
ax_ani[1].set_ylabel('')
ax_ani[1].set_ylim(v_range)

ax_ani[1].set_yticks(ticks = np.arange(v_range[0], v_range[1] + 1, 50), 
                     labels = [None, None, 0])
ax_ani[1].set_yticks(np.arange(v_range[0], v_range[1] + 1, 25), minor = True)

ax_ani[1].set_xticks(np.arange(0, 30 + 1, 30))
ax_ani[1].set_xticks(np.arange(0, 30 + 1, 1), minor = True)

[ax.grid(False) for ax in ax_ani]
[ax.spines[spine].set_visible(False) for ax in ax_ani for spine in ['top', 'right']]


# animation function 
interval = int(SR / 10)  # in datapoints

def animate(frame):
    # end frame for each update in the animation
    # end_frame = (frame + 1) * interval
    print(frame)
    
    end_frame = interval * (frame + 1)
    
    if frame > 0:
        line_i.set_data(t[:end_frame], i[:end_frame])
        line_v.set_data(t[:end_frame], v[:end_frame])
    
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

# animate(100)

# %% animation settings

# Animation parameters
duration_s = t_total  # Select the length of the video here
# duration_s = 10
 
anim_fps = frames/duration_s
 
extra_args = [
    '-vcodec', 'libx264',  # Video codec
    ]  # Set the output resolution, same ratio as figure size
 

#%% MP4

# Create the video with FFMpegWriter
writer = animation.FFMpegWriter(fps=anim_fps, bitrate=3000, 
                                codec="h264",  extra_args=extra_args)
 
# Save path
anim.save(f'{figure_dir}/{cell_ID}-{PGF}-animation-realtime.mp4', writer = writer)

# anim.save('C:/Users/nesseler/Desktop/local E-Phys/figures/video.mp4', writer=writer)
 
# Print the fps and duration of the final video
print('Video_duration (s):', frames/anim_fps)
print('Video_fps:', anim_fps)

# %% gif
# writergif = animation.PillowWriter(fps=anim_fps, bitrate=2000)
 
# anim.save(f'C:/Users/nesseler/Desktop/local E-Phys/figures/{cell_ID}-{PGF}-animation-realtime.gif', writer=writergif)

# print('Done!')









    
    
    