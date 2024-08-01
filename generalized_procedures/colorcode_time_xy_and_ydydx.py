#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 13:31:48 2024

@author: moritznesseler
"""

import numpy as np
import matplotlib as mtl
import matplotlib.pyplot as plt


plt.style.use('default')
plt.style.use('dark_background')


# x (time)
x = np.arange(0, 2 * np.pi * 5, 0.01)

# y (amplitude)
# calc sine wave
sine_pure = np.sin(x)

# calculate a linear decline to multiply on the pure sine wave
decline = np.linspace(1, 0.3, num = len(x))

# multiply the pure sine wave and its linear decline
y = np.multiply(sine_pure, decline)
 
#calculate the first derivate dy/dx
dy = np.diff(y)
dx = np.diff(x)
dydx = dy / dx

# pad with first value as nan
dydx = np.pad(dydx,
              pad_width = (1, 0),
              mode = 'constant',
              constant_values = (np.nan,))

### time color coding ###

# specify color map
cmap_str = 'viridis'

# min max normalize time for color-code
norm = mtl.colors.Normalize(x.min(), x.max())

# create mappable colormap object for colorbar
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)


## create first line collection for xy
# create individual points as a combination of their x and y component
points_xy = np.array([x, y]).T.reshape(-1, 1, 2)

# concatenate these points to form the segments needed for plotting
segments_xy = np.concatenate([points_xy[:-1], points_xy[1:]], axis=1)

# create a line collection
# solution for plotting lots of individual lines
lc_xy = mtl.collections.LineCollection(segments_xy, cmap=cmap_str, norm=norm)

# Set the values used for colormapping
# in this case use x (time) component
lc_xy.set_array(x)


## create second linecollection of ydydx
# create individual points as a combination of their x and y component
points_ydydx = np.array([y, dydx]).T.reshape(-1, 1, 2)

# concatenate these points to form the segments needed for plotting
segments_ydydx = np.concatenate([points_ydydx[:-1], points_ydydx[1:]], axis=1)

# create a line collection
# solution for plotting lots of individual lines
lc_ydydx = mtl.collections.LineCollection(segments_ydydx, cmap=cmap_str, norm=norm)

# Set the values used for colormapping
# in this case use x (time) component
lc_ydydx.set_array(x)



### plotting ###
fig, axs = plt.subplots(nrows = 1, ncols = 2,
                        layout = 'constrained')

# colorbar
fig.colorbar(cmap, ax = axs[-1], orientation = 'vertical')

# plot amplitude vs time
lines_xy = axs[0].add_collection(lc_xy)

axs[0].set_xlim([0, x.max()])
axs[0].set_ylim([-1, 1])

# plot amp v dydx (equivalent to phase plane plot)
lines_ydydx = axs[1].add_collection(lc_ydydx)

axs[1].set_xlim([-1, 1])
axs[1].set_xticks([-1, 0, 1])
axs[1].set_ylim([-1, 1])


# format axis
for ax in axs:
    ax.set_yticks([-1, 0, 1])