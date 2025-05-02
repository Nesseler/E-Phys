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
miniML_path = synaptic_dir + '/miniML_validation'

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

 
# # %% get data

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

# n_events.to_excel(synaptic_dir + '/n_events.xlsx', index_label = 'winsize')
# avg_score.to_excel(synaptic_dir + '/avg_score.xlsx', index_label = 'winsize')

# %% read xlsx

n_events = pd.read_excel(synaptic_dir + '/n_events.xlsx', index_col = 'winsize')
avg_score = pd.read_excel(synaptic_dir + '/avg_score.xlsx', index_col = 'winsize')

# %% figure

# init plotting
from functions.initialize_plotting import *

# init figure
fig, axs = plt.subplots(nrows = 2, ncols = 1,
                        sharex = True,
                        dpi = 300)

axs[0].plot(winsizes, n_events,
            lw = 0.75,
            color = 'gray',
            alpha = 0.7,
            marker = '.',
            ms = 2,
            label = 'recordings',
            zorder = 0)

axs[0].errorbar(x = winsizes, 
                y = n_events.mean(axis = 1),
                lw = 1,
                color = 'k',
                ls = 'dashed',
                elinewidth = 0.75,
                capsize = 1,
                capthick = 0.75,
                label = 'mean',
                zorder = 2)

axs[0].fill_between(x = winsizes,
                    y1 = (n_events.mean(axis = 1) - n_events.std(axis = 1)).to_list(),
                    y2 = (n_events.mean(axis = 1) + n_events.std(axis = 1)).to_list(),
                    facecolor = 'k',
                    alpha = 0.1,
                    zorder = 1)
    

axs[1].plot(winsizes, avg_score,
            lw = 0.75,
            color = 'gray',
            alpha = 0.7,
            marker = '.',
            ms = 2,
            label = 'recordings')

axs[1].errorbar(x = winsizes, 
                y = avg_score.mean(axis = 1),
                lw = 1,
                color = 'k',
                ls = 'dashed',
                elinewidth = 0.75,
                capsize = 1,
                capthick = 0.75,
                label = 'mean')

axs[1].fill_between(x = winsizes,
                    y1 = (avg_score.mean(axis = 1) - avg_score.std(axis = 1)).to_list(),
                    y2 = (avg_score.mean(axis = 1) + avg_score.std(axis = 1)).to_list(),
                    facecolor = 'k',
                    alpha = 0.1,
                    zorder = 1)

# axs[0].legend(frameon = False,
#               fontsize = 7,
#               loc = 'upper right')

axs[0].set_ylim([0-12, 600+12])
axs[0].spines['left'].set_bounds([0, 600])
axs[0].set_ylabel('Number of\ndetected events [#]')
    
# axs[1].legend(frameon = False,
#               fontsize = 7,
#               loc = 'lower right')

axs[1].set_ylim([0.70, 1])
axs[1].set_ylabel('Average\nevent score')

# x
axs[1].set_xlabel('Sliding window size [ms]')
axs[1].set_xlim([0, 300])
axs[1].set_xticks(ticks = np.arange(0, 300+0.1, 60))
axs[1].set_xticks(ticks = np.arange(0, 300+0.1, 6), minor = True)

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

fig.align_labels()

plt.show()


# %% figure 2 - normed n_events

n_events_norm = n_events / n_events.max()

# init figure
fig, axs = plt.subplots(nrows = 2, ncols = 1,
                        sharex = True,
                        dpi = 300)

axs[0].plot(winsizes, n_events_norm,
            lw = 0.75,
            color = 'gray',
            alpha = 0.7,
            label = 'recordings',
            zorder = 0)

axs[0].errorbar(x = winsizes, 
                y = n_events_norm.mean(axis = 1),
                lw = 1,
                color = 'k',
                ls = 'dashed',
                elinewidth = 0.75,
                capsize = 1,
                capthick = 0.75,
                label = 'mean',
                zorder = 2)

axs[0].fill_between(x = winsizes,
                    y1 = (n_events_norm.mean(axis = 1) - n_events_norm.std(axis = 1)).to_list(),
                    y2 = (n_events_norm.mean(axis = 1) + n_events_norm.std(axis = 1)).to_list(),
                    facecolor = 'k',
                    alpha = 0.1,
                    zorder = 1)
    

axs[1].plot(winsizes, avg_score,
            lw = 0.75,
            color = 'gray',
            alpha = 0.7,
            label = 'recordings')

axs[1].errorbar(x = winsizes, 
                y = avg_score.mean(axis = 1),
                lw = 1,
                color = 'k',
                ls = 'dashed',
                elinewidth = 0.75,
                capsize = 1,
                capthick = 0.75,
                label = 'mean')

axs[1].fill_between(x = winsizes,
                    y1 = (avg_score.mean(axis = 1) - avg_score.std(axis = 1)).to_list(),
                    y2 = (avg_score.mean(axis = 1) + avg_score.std(axis = 1)).to_list(),
                    facecolor = 'k',
                    alpha = 0.1,
                    zorder = 1)

# axs[0].legend(frameon = False,
#               fontsize = 7,
#               loc = 'upper right')

axs[0].set_ylim([-0.05, 1.1])
axs[0].spines['left'].set_bounds([0, 1])
axs[0].set_ylabel('Normalized number of\ndetected events')
    
# axs[1].legend(frameon = False,
#               fontsize = 7,
#               loc = 'lower right')

axs[1].set_ylim([0.70, 1])
axs[1].set_ylabel('Average\nevent score')

# x
axs[1].set_xlabel('Sliding window size [ms]')
axs[1].set_xlim([0, 300])
axs[1].set_xticks(ticks = np.arange(0, 300+0.1, 60))
axs[1].set_xticks(ticks = np.arange(0, 300+0.1, 6), minor = True)

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

fig.align_labels()

plt.show()
