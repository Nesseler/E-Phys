# -*- coding: utf-8 -*-
"""
Created on Tue May 13 17:45:35 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, synaptic_dir
from parameters.parameters import PSC_bins

# import functions
from functions.functions_pscs import map_detection_to_df

# set validation cell_IDs
cell_IDs = ['E-298', 'E-301', 'E-302', 'E-303', 'E-309', 'E-310', 'E-314']

winsize = 114
th = 0.75


# %% load

# create dataframe
allevents  = pd.DataFrame()

for cell_ID in tqdm(cell_IDs):
    
    # set filename
    filename = f'miniMLdetect_{cell_ID}_Erest_ctrl_{winsize}_{str(th).replace(".", "p")}' 
    
    # open a file, where you stored the pickled data
    file = open((synaptic_dir + f'/miniML_dtc-validation/' + filename + '.pickle'), 'rb')
    
    # dump information to that file
    detection = pickle.load(file)
    
    # close and remove (from memory) the file
    file.close()
    del file 
    gc.collect()

    # write to dataframe
    allevents = map_detection_to_df(detection, df = allevents, cell_ID = cell_ID, include_events=True)

# drop nan valued rows
allevents.dropna(axis = 'index', how = 'any', inplace = True)


# %% calc histogram

score_bins = PSC_bins['score']
scores_hist, _ = np.histogram(allevents['score'], bins = score_bins)


# %% figure

# init plotting
from functions.initialize_plotting import *

fig, axs = plt.subplot_mosaic(mosaic = [[0,1,2,3],[0,4,5,6],[0,7,8,9]],
                              figsize = get_figure_size(width = 159.2, height = 80),
                              width_ratios = [4,1,1,1],
                              dpi = 600,
                              layout = 'constrained')

# set histogram axis
ax = axs[0]
ax.set_title('A: Event score histogram')

ax.stairs(values = scores_hist, 
          edges = score_bins,
          baseline = [0]*len(scores_hist-1),
          fill = False,
          lw = 0.75,
          color = colors_dict['primecolor'],
          label = 'score')

for i_sbin, sbin in enumerate(score_bins):
    if sbin >= 0.7 and sbin < 1.:
        ax.text(x = sbin+0.005,
                y = scores_hist[i_sbin] + 20,
                s = scores_hist[i_sbin],
                fontsize = 6,
                rotation = 90,
                ha = 'center', va = 'bottom')

# x
apply_axis_settings(ax, axis = 'x', ax_min=0.70, ax_max=1.0, pad=None, step=0.1, stepminor=0.02, label='Event score')

# y
apply_axis_settings(ax, axis = 'y', ax_min=0, ax_max=1800, pad=None, step=500, stepminor=100, label='Event count [#]')

# full range inset axis
ax_inset = ax.inset_axes([0.15, 0.65, 0.30, 0.30])

# plot in inset
ax_inset.stairs(values = scores_hist, 
                edges = score_bins,
                fill = False,
                lw = 0.5,
                color = colors_dict['primecolor'],
                label = 'score')

# add rectangle marker
ax_inset.add_patch(Rectangle(xy = (0.7, 0-40), 
                   width = 1-0.7+0.02, 
                   height = 1800+40,
                   fill = False,
                   color = 'r',
                   linestyle = '--',
                   lw = 0.5))

# x
apply_axis_settings(ax_inset, axis = 'x', ax_min=0, ax_max=1.0, pad=0.05, step=1, stepminor=0.1, label='')

# y
apply_axis_settings(ax_inset, axis = 'y', ax_min=0, ax_max=1800, pad=80, step=500, stepminor=100, label='', ticklabels=['', '', '', ''])


# event plots

ax = axs[1]








# align labels
fig.align_labels()

# display figure
plt.show()


# %%

event_t = np.arange(0, len(allevents.at[0, 'event']) / (100_000/1e3), 1/(100_000/1e3))

# calculate points added with formulas from miniML
add_points = int(19*600/3)
after = 144 + add_points
before = add_points

event_i = 1179
event_trace = allevents.at[event_i, 'event']
peak_idx = before + (allevents.at[event_i, 'peak_idx'] - allevents.at[event_i, 'event_idx'])
onset_idx = before + (allevents.at[event_i, 'onset_idx'] - allevents.at[event_i, 'event_idx'])
half_decayidx = int(before + (allevents.at[event_i, 'halfdecay_idx'] - allevents.at[event_i, 'event_idx']))

plt.plot(event_t, event_trace)

plt.scatter(event_t[peak_idx], event_trace[peak_idx], c = 'r')
plt.scatter(event_t[onset_idx], event_trace[onset_idx], c = 'g')
plt.scatter(event_t[half_decayidx], event_trace[half_decayidx], c = 'y')


# %%

# grouped_allevents = {cell_ID : group for cell_ID, group in allevents.groupby(by = 'cell_ID')}

# cell_allevents = grouped_allevents['E-303']

# # for i in cell_allevents.query('score > 0.80 and score < 0.99').index:
# for i in cell_allevents.query('score >= 0.70 and score < 0.80').index:
    
#     if cell_allevents.at[i, 'amplitude'] <= 0:
#     # if cell_allevents.at[i, 'amplitude'] >= -15 and cell_allevents.at[i, 'amplitude'] <= -5:
#     # if cell_allevents.at[i, 'amplitude'] >= -5:
        
#         plt.plot(cell_allevents.at[i, 'event'])
#         plt.title(f'{allevents.at[i, "cell_ID"]} event i : {i}')
#         plt.show()


# %% figure - scatter ampl

fig, axs = plt.subplots(nrows = 2, ncols = 2,
                        figsize = get_figure_size(width = 120, height = 100),
                        dpi = 600,
                        layout = 'constrained',
                        height_ratios = [1, 7],
                        width_ratios = [7, 1],
                        sharex = 'col',
                        sharey = 'row')

# flatten axis array
axs = axs.flatten()

# scatter 
ax = axs[2]

# scatter plot
ax.scatter(x = allevents['amplitude'],
           y = allevents['score'],
           marker = '.',
           s = 3,
           alpha = 0.25,
           color = 'k')

# half violin 1
plot_half_violin(allevents['amplitude'], 
                 ax = axs[0],
                 v_orientation = 'horizontal',
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1,
                 v_kde_cutoff = np.nan,
                 v_abs_cutoff = [-50, 0])

# half violin 2
plot_half_violin(allevents['score'], ax = axs[3],
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1,
                 v_kde_cutoff = np.nan,
                 v_abs_cutoff = [0.5, 1.01])

# # x
apply_axis_settings(axs[0], axis = 'x', ax_min=-50, ax_max=0, pad=None, step=10, stepminor=2, label='')
apply_axis_settings(axs[2], axis = 'x', ax_min=-50, ax_max=0, pad=None, step=10, stepminor=2, label='Amplitude [mV]')

# y
apply_axis_settings(axs[2], axis = 'y', ax_min=0.5, ax_max=1., pad=0.02, step=0.1, stepminor=0.02, label='Event score')
apply_axis_settings(axs[3], axis = 'y', ax_min=0.5, ax_max=1., pad=0.02, step=0.1, stepminor=0.02, label='')
    
# remove spines
remove_spines_ticks_labels([axs[0]], axis = 'y')
remove_spines_ticks_labels([axs[3]], axis = 'x')

# remove axis
plt.delaxes(axs[1])

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = 'figure-miniML_validation-ampl_score-scatter', 
             save_dir = figure_dir + '/miniML_validation',
             figure_format = 'both')


# %% figure - scatter ampl & charge

fig, axs = plt.subplots(nrows = 2, ncols = 3,
                        figsize = get_figure_size(width = 159.2, height = 85),
                        dpi = 600,
                        layout = 'constrained',
                        height_ratios = [1, 7],
                        width_ratios = [7, 7, 1],
                        sharex = 'col',
                        sharey = 'row')

# flatten axis array
axs = axs.flatten()

# scatter 1
ax = axs[3]

# scatter plot 1
ax.scatter(x = allevents['amplitude'],
           y = allevents['score'],
           marker = '.',
           s = 3,
           alpha = 0.25,
           color = 'k')

# half violin x 1
plot_half_violin(allevents['amplitude'], 
                 ax = axs[0],
                 v_orientation = 'horizontal',
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1,
                 v_kde_cutoff = np.nan,
                 v_abs_cutoff = [-50, 0])

# x 1
apply_axis_settings(axs[0], axis = 'x', ax_min=-50, ax_max=0, pad=None, step=10, stepminor=2, label='')
apply_axis_settings(axs[3], axis = 'x', ax_min=-50, ax_max=0, pad=None, step=10, stepminor=2, label='Amplitude [mV]')

# scatter 2
ax = axs[4]

# scatter plot 2
ax.scatter(x = allevents['charge'],
           y = allevents['score'],
           marker = '.',
           s = 3,
           alpha = 0.25,
           color = 'k')

# half violin x 2
plot_half_violin(allevents['charge'], 
                 ax = axs[1],
                 v_orientation = 'horizontal',
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1,
                 v_kde_cutoff = np.nan,
                 v_abs_cutoff = [-0.8, 0])

# x 2
apply_axis_settings(axs[1], axis = 'x', ax_min=-0.8, ax_max=0.0, pad=None, step=0.2, stepminor=0.02, label='')
apply_axis_settings(axs[4], axis = 'x', ax_min=-0.8, ax_max=0.0, pad=None, step=0.2, stepminor=0.02, label='Charge [pC]')

# half violin y
plot_half_violin(allevents['score'], ax = axs[5],
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1,
                 v_kde_cutoff = np.nan,
                 v_abs_cutoff = [0.5, 1.01])


# # y
apply_axis_settings(axs[3], axis = 'y', ax_min=0.5, ax_max=1., pad=0.02, step=0.1, stepminor=0.02, label='Event score')
apply_axis_settings(axs[4], axis = 'y', ax_min=0.5, ax_max=1., pad=0.02, step=0.1, stepminor=0.02, label='')
apply_axis_settings(axs[5], axis = 'y', ax_min=0.5, ax_max=1., pad=0.02, step=0.1, stepminor=0.02, label='')   
 
# # remove spines
remove_spines_ticks_labels([axs[0]], axis = 'y')
remove_spines_ticks_labels([axs[1]], axis = 'y')
remove_spines_ticks_labels([axs[5]], axis = 'x')

# remove axis
plt.delaxes(axs[2])

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = 'figure-miniML_validation-ampl_charge_score-scatter', 
             save_dir = figure_dir + '/miniML_validation',
             figure_format = 'both')

