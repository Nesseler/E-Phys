# -*- coding: utf-8 -*-
"""
Created on Wed May 21 16:38:22 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, synaptic_dir

# set parameter
cell_IDs = ['E-298', 'E-301', 'E-302', 'E-303', 'E-309', 'E-310', 'E-314']
ehold = 'Erest'
treatment = 'ctrl'
winsize = 114
th = 0.5

traces = dict.fromkeys(cell_IDs)
peak_idx = dict.fromkeys(cell_IDs)

same_y = False

# %% load

for cell_ID in tqdm(cell_IDs):
    
    # set filename
    filename = f'miniMLdetect_{cell_ID}_{ehold}_{treatment}_{winsize}_{str(th).replace(".", "p")}' 
    
    # open a file, where you stored the pickled data
    file = open((synaptic_dir + f'/miniML_dtc-validation/' + filename + '.pickle'), 'rb')
    
    # dump information to that file
    detection = pickle.load(file)
    
    # close and remove (from memory) the file
    file.close()
    del file 
    gc.collect()

    # write to dict
    peak_idx[cell_ID] = detection['event_location_parameters']['event_peak_locations']
    traces[cell_ID] = detection['mini_trace']


# %% figure

# init plotting
from functions.initialize_plotting import *

# init figure
fig, axs = plt.subplots(nrows = len(cell_IDs), ncols = 1,
                        figsize = get_figure_size(width = 159.2, height = 246.2),
                        dpi = 300,
                        layout = 'constrained',
                        sharex = True)

# set figure title
fig.suptitle('Validation traces')

for ci, cell_ID in enumerate(cell_IDs):
    
    # set titles
    axs[ci].set_title(f'{cell_ID}')
    
    # calc x
    x_full = np.arange(0, traces[cell_ID].shape[0] / 100_000, step=1/100_000)

    # trace
    axs[ci].plot(x_full, traces[cell_ID],
                 color='k',
                 lw=0.5,
                 label='data (filt)',
                 zorder = 0)

    ### event peak indicators
    axs[ci].scatter(x = x_full[peak_idx[cell_ID]],
                    y = traces[cell_ID][peak_idx[cell_ID]],
                    marker = 'o',
                    color = 'r',
                    s = 2,
                    lw = 1,
                    zorder = 1)

# axes
for ax in axs:
    apply_axis_settings(ax, axis = 'x', ax_min=0, ax_max=180, pad=None, step=30, stepminor=5, label='')
axs[-1].set_xlabel('Time [s]')

# y
if not same_y:
    for ci, cell_ID in enumerate(cell_IDs):
        if cell_ID == 'E-298':
            apply_axis_settings(axs[ci], axis='y', ax_min=-40, ax_max=10, pad=None, step=20, stepminor=5, label='')
        elif cell_ID == 'E-301':
            apply_axis_settings(axs[ci], axis='y', ax_min=-60, ax_max=0, pad=None, step=20, stepminor=5, label='')
        elif cell_ID == 'E-302':
            apply_axis_settings(axs[ci], axis='y', ax_min=-40, ax_max=10, pad=None, step=20, stepminor=5, label='') 
        elif cell_ID == 'E-303':
            apply_axis_settings(axs[ci], axis='y', ax_min=-30, ax_max=10, pad=None, step=10, stepminor=5, label='')
        elif cell_ID == 'E-309':
            apply_axis_settings(axs[ci], axis='y', ax_min=-40, ax_max=15, pad=None, step=20, stepminor=5, label='') 
        elif cell_ID == 'E-310':
            apply_axis_settings(axs[ci], axis='y', ax_min=-20, ax_max=15, pad=None, step=20, stepminor=5, label='') 
        elif cell_ID == 'E-314':
            apply_axis_settings(axs[ci], axis='y', ax_min=-40, ax_max=5, pad=None, step=20, stepminor=5, label='') 
            
elif same_y:
    for ax in axs:
        apply_axis_settings(ax, axis='y', ax_min=-60, ax_max=10, pad=None, step=20, stepminor=5, label='') 

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = f'figure-miniML_validation-alltraces', 
             save_dir = figure_dir + '/miniML_validation',
             figure_format = 'both')
