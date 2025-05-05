# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 14:44:52 2025

@author: nesseler
"""

# import packages
from functions.initialize_packages import *
import pickle
import gc

# import directories 
from parameters.directories_win import synaptic_dir
miniML_path = synaptic_dir + '/miniML_dtc-validation'

# import functions
from functions.functions_import import get_onlyfiles_list

# get parameter from filenames
onlyfiles = get_onlyfiles_list(miniML_path)

# init sets
cell_IDs = set()
winsizes = set()
ths = set()

# iterate through files list
for filename in onlyfiles:
    
    # remove filename extension
    filename = splitext(filename)[0]
    
    # split filename string
    _, cell_ID, hold, treat, winsize_ms, th = filename.split('_')

    # convert to usable numbers
    th = float(th.replace('p', '.'))
    winsize_ms = int(winsize_ms)

    # add to set if not already included
    cell_IDs.add(cell_ID)
    winsizes.add(winsize_ms)
    ths.add(th)
    
# convert to list
cell_IDs = list(cell_IDs)
winsizes = list(winsizes)
ths = list(ths)

# sort lists
for ls in [cell_IDs, winsizes, ths]:
    ls.sort()

# init dataframes
n_events = pd.DataFrame(columns = cell_IDs, index = winsizes)
avg_score = pd.DataFrame(columns = cell_IDs, index = winsizes)

for df in [n_events, avg_score]:
    df.index.name = 'winsize'

 
# %% get data

# th = ths[0]

# for cell_ID in tqdm(cell_IDs):
    
#     for winsize in winsizes:
        
#         # set filename
#         filename = f'miniMLdetect_{cell_ID}_Erest_ctrl_{str(winsize)}_{str(th).replace(".", "p")}' 
        
#         # open a file, where you stored the pickled data
#         file = open((miniML_path + '/' +filename + '.pickle'), 'rb')
        
#         # dump information to that file
#         detection = pickle.load(file)
        
#         # close and remove (from memory) the file
#         file.close()
#         del(file)
#         gc.collect()

#         # get data
#         n_events.at[winsize, cell_ID] = detection['event_location_parameters']['event_peak_locations'].shape[0]
#         avg_score.at[winsize, cell_ID] = np.mean(detection['event_location_parameters']['event_scores'])

# n_events.to_excel(synaptic_dir + '/miniML_validation' + '/n_events.xlsx', index_label = 'winsize')
# avg_score.to_excel(synaptic_dir + '/miniML_validation' + '/avg_score.xlsx', index_label = 'winsize')

# %% read xlsx

n_events = pd.read_excel(synaptic_dir + '/miniML_validation' + '/n_events.xlsx', index_col = 'winsize')
avg_score = pd.read_excel(synaptic_dir + '/miniML_validation' + '/avg_score.xlsx', index_col = 'winsize')

# %% figure

# init plotting
from functions.initialize_plotting import *

# # init figure
# fig, axs = plt.subplots(nrows = 2, ncols = 1,
#                         sharex = True,
#                         dpi = 300)

# axs[0].plot(winsizes, n_events,
#             lw = 0.75,
#             color = 'gray',
#             alpha = 0.7,
#             marker = '.',
#             ms = 2,
#             label = 'recordings',
#             zorder = 0)

# axs[0].errorbar(x = winsizes, 
#                 y = n_events.mean(axis = 1),
#                 lw = 1,
#                 color = 'k',
#                 ls = 'dashed',
#                 elinewidth = 0.75,
#                 capsize = 1,
#                 capthick = 0.75,
#                 label = 'mean',
#                 zorder = 2)

# axs[0].fill_between(x = winsizes,
#                     y1 = (n_events.mean(axis = 1) - n_events.std(axis = 1)).to_list(),
#                     y2 = (n_events.mean(axis = 1) + n_events.std(axis = 1)).to_list(),
#                     facecolor = 'k',
#                     alpha = 0.1,
#                     zorder = 1)
    

# axs[1].plot(winsizes, avg_score,
#             lw = 0.75,
#             color = 'gray',
#             alpha = 0.7,
#             marker = '.',
#             ms = 2,
#             label = 'recordings')

# axs[1].errorbar(x = winsizes, 
#                 y = avg_score.mean(axis = 1),
#                 lw = 1,
#                 color = 'k',
#                 ls = 'dashed',
#                 elinewidth = 0.75,
#                 capsize = 1,
#                 capthick = 0.75,
#                 label = 'mean')

# axs[1].fill_between(x = winsizes,
#                     y1 = (avg_score.mean(axis = 1) - avg_score.std(axis = 1)).to_list(),
#                     y2 = (avg_score.mean(axis = 1) + avg_score.std(axis = 1)).to_list(),
#                     facecolor = 'k',
#                     alpha = 0.1,
#                     zorder = 1)

# # axs[0].legend(frameon = False,
# #               fontsize = 7,
# #               loc = 'upper right')

# axs[0].set_ylim([0-12, 600+12])
# axs[0].spines['left'].set_bounds([0, 600])
# axs[0].set_ylabel('Number of\ndetected events [#]')
    
# # axs[1].legend(frameon = False,
# #               fontsize = 7,
# #               loc = 'lower right')

# axs[1].set_ylim([0.70, 1])
# axs[1].set_ylabel('Average\nevent score')

# # x
# axs[1].set_xlabel('Sliding window size [ms]')
# axs[1].set_xlim([0, 300])
# axs[1].set_xticks(ticks = np.arange(0, 300+0.1, 60))
# axs[1].set_xticks(ticks = np.arange(0, 300+0.1, 6), minor = True)

# # remove spines
# [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# fig.align_labels()

# plt.show()


# %% figure 2 - normed n_events

# n_events_norm = (n_events-n_events.min()) / (n_events.max()-n_events.min())
n_events_norm = (n_events-0) / (n_events.max()-0)
avg_score_norm = (avg_score-0) / (avg_score.max()-0)

# init figure
fig, axs = plt.subplot_mosaic(mosaic = 'AB;AC;DE;DF',
                              dpi = 300,
                              width_ratios = [3, 1],
                              figsize = get_figure_size(width = 159.2, height = 159.2), # powerpoint width = 318.67, height = 160.5 # A4 width = 159.2, height = 159.2
                              layout = 'constrained')

fig.suptitle('miniML validation - window size')

for ax_key, data in zip(['A', 'B', 'C', 'D', 'E', 'F'], [n_events_norm, n_events_norm, n_events, avg_score_norm, avg_score_norm, avg_score_norm]):
    # set axis
    ax = axs[ax_key]
    
    # plot
    ax.plot(winsizes, data,
            lw = 0.75, 
            color = 'gray', 
            alpha = 0.5,
            label = 'recordings', 
            zorder = 0)
    
    ax.errorbar(x = winsizes, 
                y = data.mean(axis = 1),
                lw = 1,
                color = 'k',
                ls = 'solid',
                elinewidth = 0.75,
                capsize = 1,
                capthick = 0.75,
                label = 'mean',
                zorder = 2,
                marker = '.',
                ms = 2)
    
    ax.fill_between(x = winsizes,
                    y1 = (data.mean(axis = 1) - data.std(axis = 1)).to_list(),
                    y2 = (data.mean(axis = 1) + data.std(axis = 1)).to_list(),
                    facecolor = 'k',
                    alpha = 0.1,
                    zorder = 1)
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    

for ax_key, y, h, s in zip(['A', 'D'], [0.8, 0.94], [0.2, 0.14], ['B', 'E']):
    axs[ax_key].set_title(ax_key, loc='left')
    
    # x axis
    axs[ax_key].spines['bottom'].set_bounds([0, 300])
    axs[ax_key].set_xlim([0-3, 300+3])
    axs[ax_key].set_xticks(ticks = np.arange(0, 300+0.1, 60), labels = [])
    axs[ax_key].set_xticks(ticks = np.arange(0, 300+0.1, 6), minor = True)
    
    # add zoom in marker
    axs[ax_key].add_patch(Rectangle(xy = (84, y), 
                            width = 48, 
                            height = h,
                            fill = False,
                            color = colors_dict['primecolor'],
                            linestyle = '--',
                            lw = 0.5,
                            alpha = 0.5))
    
    # add text
    axs[ax_key].text(x = 132-1, y = y, s = s, fontsize = 7, ha = 'right', va = 'bottom')

for ax_key in ['B', 'E']:
    axs[ax_key].set_title(ax_key, loc='left')
    # x axis
    axs[ax_key].spines['bottom'].set_bounds([84, 132])
    axs[ax_key].set_xlim([84-2, 132+2])
    axs[ax_key].set_xticks(ticks = np.arange(84, 132+0.1, 12))
    axs[ax_key].set_xticks(ticks = np.arange(84, 132+0.1, 6), minor = True)  
    axs[ax_key].set_xlabel('Window size [ms]')

for ax_key in ['C', 'F']:
    axs[ax_key].set_title(ax_key, loc='left')
    # x axis
    axs[ax_key].spines['bottom'].set_bounds([0, 300])
    axs[ax_key].set_xlim([0-3, 300+3])
    axs[ax_key].set_xticks(ticks = np.arange(0, 300+0.1, 60), labels = np.arange(0, 300+0.1, 60, dtype = int))
    axs[ax_key].set_xlabel('Window size [ms]')

    
# A: number of events  
axs['A'].set_ylim([-0.01, 1.01])
axs['A'].spines['left'].set_bounds([0, 1])
axs['A'].set_ylabel('Normalized\nnumber of detected events')
axs['A'].set_yticks(ticks = np.arange(0, 1.+0.001, 0.2))
axs['A'].set_yticks(ticks = np.arange(0, 1.+0.001, 0.05), minor = True)

# B: zoom in
axs['B'].set_ylim([0.8-0.01, 1.01])
axs['B'].spines['left'].set_bounds([0.8, 1])
axs['B'].set_ylabel('Normalized number\nof detected events')
axs['B'].set_yticks(ticks = np.arange(0.8, 1.+0.001, 0.1))
axs['B'].set_yticks(ticks = np.arange(0.8, 1.+0.001, 0.05), minor = True)

# C: abs number 
axs['C'].set_ylim([0-5, 600+5])
axs['C'].spines['left'].set_bounds([0, 600])
axs['C'].set_ylabel('Number of\ndetected events [#]')
axs['C'].set_yticks(ticks = np.arange(0, 600.+1, 200))
axs['C'].set_yticks(ticks = np.arange(0, 600.+1, 50), minor = True)

# D: average event score
axs['D'].set_xticks(ticks = np.arange(0, 300+0.1, 60), labels = np.arange(0, 300+0.1, 60, dtype = int))
axs['D'].set_xlabel('Window size [ms]')

axs['D'].set_ylim([0.0-0.003, 1+0.003])
axs['D'].spines['left'].set_bounds([0.0, 1])
axs['D'].set_ylabel('Normalized\naverage event score')
axs['D'].set_yticks(ticks = np.arange(0.0, 1.+0.001, 0.2))
axs['D'].set_yticks(ticks = np.arange(0.0, 1.+0.001, 0.05), minor = True)

# E: zoom in
axs['E'].set_ylim([0.94-0.003, 1+0.003])
axs['E'].spines['left'].set_bounds([0.94, 1])
axs['E'].set_ylabel('Normalized\naverage event score')
axs['E'].set_yticks(ticks = np.arange(0.94, 1.+0.001, 0.06))
axs['E'].set_yticks(ticks = np.arange(0.94, 1.+0.001, 0.02), minor = True)

# F: abs average event score swarm
axs['F'].set_ylim([0.7-0.003, 1+0.003])
axs['F'].spines['left'].set_bounds([0.7, 1])
axs['F'].set_ylabel('Average\nevent score')
axs['F'].set_yticks(ticks = np.arange(0.80, 1.+0.001, 0.2))
axs['F'].set_yticks(ticks = np.arange(0.70, 1.+0.001, 0.05), minor = True)

# labels
fig.align_labels()

# display figure
plt.show()

# import directories 
from parameters.directories_win import figure_dir
figure_path = figure_dir + '/miniML_validation'
save_figures(fig, 
             figure_name = 'miniML_validation-window_size', 
             save_dir = figure_path,
             figure_format = 'both')

