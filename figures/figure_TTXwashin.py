# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 13:39:23 2025

@author: nesseler
"""

# import standard packages
from functions.initialize_packages import *

# initialize plotting
from functions.initialize_plotting import *

# get directories
from parameters.directories_win import quant_data_dir, figure_dir

# load data
peakCurrents = pd.read_excel(join(quant_data_dir, 'TTX_washin', 'TTX_washin-peakCurrents.xlsx'), 
                             index_col = 'step_idc')

# get cell IDs
from functions.get_cell_IDs import get_cell_IDs_one_protocol
cell_IDs = get_cell_IDs_one_protocol('vc_TTX_washin', sheet_name = 'PGFs_Syn')
cell_IDs_leak = get_cell_IDs_one_protocol('vc_TTX_washin_leak', sheet_name = 'PGFs_Syn')

# %%

# min max normalize currents
mm_peakCurrents = (peakCurrents - peakCurrents.min()) / (peakCurrents.max() - peakCurrents.min())

plt_cell_IDs = cell_IDs + cell_IDs_leak
# plt_cell_IDs = cell_IDs_leak

#  plot graph
fig, ax = plt.subplots(nrows = 1,
                       ncols = 1,
                       figsize = get_figure_size(width = 318.67/2),
                       dpi = 300,
                       layout = 'constrained')

# create x dimension
x = np.arange(0, 120) * 5 + 2.5

for cell_ID in plt_cell_IDs:
    
    
    if cell_ID in cell_IDs_leak:
        linestyle = 'dashed'
    else:
        linestyle = 'solid'

    # plot
    ax.plot(x, mm_peakCurrents[cell_ID],
            c = 'grey',
            lw = 1,
            ls = linestyle)
    
# calc mean
mean_mm_peakcurrent = mm_peakCurrents[plt_cell_IDs].mean(axis = 1)
std_mm_peakcurrent = mm_peakCurrents[plt_cell_IDs].std(axis = 1)


# plot mean and std
ax.plot(x, mean_mm_peakcurrent,
        c = colors_dict['primecolor'],
        lw = 2)

# y axis
apply_axis_settings(ax = ax, 
                    axis = 'y', 
                    ax_min = 0, 
                    ax_max = 1, 
                    pad = 0.05, 
                    step = 0.2, 
                    stepminor = 0.1, 
                    label = 'Current [pA]' )

# x axis
apply_axis_settings(ax = ax, 
                    axis = 'x', 
                    ax_min = 0, 
                    ax_max = 600, 
                    pad = 30, 
                    step = 60, 
                    stepminor = 30, 
                    label = 'Time [s]' )
# despine
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# show figure
plt.show()

# set directory for figure
fig_dir = join(figure_dir, 'temp_figs')

# save figure
save_figures(fig, 'figure-TTXwashin-minmax', 
              save_dir = fig_dir,
              darkmode_bool = darkmode_bool,
              figure_format = 'both')