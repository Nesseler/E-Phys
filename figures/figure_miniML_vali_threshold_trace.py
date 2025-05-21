# -*- coding: utf-8 -*-
"""
Created on Wed May 21 15:09:32 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, synaptic_dir

# set parameter
cell_ID = 'E-303'
ehold = 'Erest'
treatment = 'ctrl'
winsize = 114

ths = [0.5, 0.75, 0.9]
event_idc = dict.fromkeys(ths)


# %% load

for i, th in enumerate(ths):
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
    event_idc[th] = detection['event_location_parameters']['event_peak_locations']
    
    # get trace and prediction
    if i == 0:
        trace = detection['mini_trace']
        prediction = detection['prediction']
        
        # filter prediction
        filtered_prediction = sc.ndimage.maximum_filter1d(prediction, size=int(5*detection['metadata']['interpolation_factor']), origin=-2)

        # fill with nan values
        prediction = np.append(filtered_prediction, [np.nan]*(trace.shape[0]-filtered_prediction.shape[0]))


# %%

t = 5
trange = 10
SR = 100_000

# init plotting
from functions.initialize_plotting import *
# th_colors = {0.5 : '#292f56', 0.75 : '#00a3a4', 0.9 : '#acfa70'}
th_colors = {0.5 : '#0F4F69', 0.75 : '#2AB8A2', 0.9 : '#61BA8D'}

# init figure
fig, axs = plt.subplots(nrows = 3, ncols = 1,
                        figsize = get_figure_size(width = 159.2, height = 78.73),
                        dpi = 300,
                        layout = 'constrained',
                        height_ratios = [3, 1, 5],
                        sharex = True)

# set figure title
fig.suptitle(f'{cell_ID} - Prediction - Trace',
            fontsize = 9)

# set plotting 
x_full = np.arange(0, trace.shape[0] / SR, step=1/SR)
plt_idc = np.arange((0+t) * SR, (trange+t) * SR, step = 1, dtype = int)
plt_peak_idc = {th : [idx for idx in peak_idx if (idx > plt_idc[0] and idx < plt_idc[-1])] for th, peak_idx in event_idc.items()}

### prediction
axs[0].plot(x_full[plt_idc], prediction[plt_idc],
            color='k',
            lw=0.75,
            label='filtered',
            zorder = 0)

### trace
axs[2].plot(x_full[plt_idc], trace[plt_idc],
            color='k',
            lw=0.5,
            label='data (filt)',
            zorder = 0)

for i_th, th in enumerate(ths):
    ### threshold line
    axs[0].hlines(xmin = x_full[plt_idc][0], xmax = x_full[plt_idc][-1],
                  y = th, 
                  color = th_colors[th],
                  lw = 0.75,
                  ls = 'dashed',
                  zorder = 1)
    
    ### threshold text
    axs[0].text(x = t,
                y = th,
                s = th,
                color = th_colors[th],
                fontsize = 6,
                ha = 'left', 
                va = 'bottom')

    ### eventplot
    axs[1].eventplot(positions = x_full[plt_peak_idc[th]],
                     orientation = 'horizontal',
                     lineoffsets = i_th,
                     linelengths = 0.8,
                     linewidths = 1,
                     color = th_colors[th],
                     label = 'events')

    ### event peak indicators
    axs[2].scatter(x = x_full[plt_peak_idc[th]],
                   y = trace[plt_peak_idc[th]],
                   marker = 'o',
                   color = th_colors[th],
                   s = [13, 7, 2][i_th],
                   lw = 1,
                   zorder = i_th+1)


# axes
for ax in axs:
    apply_axis_settings(ax, axis = 'x', ax_min=t, ax_max=(t+trange), pad=None, step=5, stepminor=1, label='')
axs[2].set_xlabel('Time [s]')  

# y
axs[0].set_ylim([-0.05, 1.05])
axs[0].set_ylabel('Probability')
axs[0].spines['left'].set_bounds([0, 1])
axs[1].set_ylim([-1, 3])
axs[1].set_ylabel('Events')
axs[1].set_yticks(ticks = [0, 1, 2], labels = [])
axs[1].spines['left'].set_bounds([0, 2])
apply_axis_settings(axs[2], axis = 'y', ax_min=-20, ax_max=5, pad=None, step=10, stepminor=2, label='Current [pA]')

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = f'figure-miniML_validation-score_threshold-{t}_{trange}', 
             save_dir = figure_dir + '/miniML_validation',
             figure_format = 'both')