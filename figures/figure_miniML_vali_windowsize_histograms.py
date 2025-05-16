# -*- coding: utf-8 -*-
"""
Created on Fri May 16 16:01:33 2025

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
cell_ID = cell_IDs[3]

th = 0.5
winsizes = [96, 114]

# %% load

# create dataframe
allevents  = pd.DataFrame()

# set previouse length
pre_len = 0

for winsize in winsizes:
    for cell_ID in cell_IDs:
        
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
        allevents = map_detection_to_df(detection, df = allevents, cell_ID = cell_ID)
    
        # get number of events
        n_events = detection['individual_values']['amplitudes'].shape[0]
        
        # get current length of allevents
        allevents_len = allevents.shape[0]
        
        # indexing
        allevents_idx = np.arange(pre_len, allevents_len, dtype = int)
    
        # extend with windows size
        allevents.loc[allevents_idx, 'winsize'] = [winsize] * int(allevents_idx.shape[0])
        
        # update previouse length
        pre_len = allevents_len

# drop nan valued rows
allevents.dropna(axis = 'index', how = 'any', inplace = True)


# %% calc histogram

score_bins = PSC_bins['scores']
score_hist = pd.DataFrame(index = score_bins[:-1], columns = winsizes)
ampl_bins = PSC_bins['amplitudes']
ampl_hist = pd.DataFrame(index = ampl_bins[:-1], columns = winsizes)

for winsize in winsizes:
    
    # get occurances
    scores_hist_win, _ = np.histogram(allevents[allevents['winsize'] == winsize]['score'], bins = score_bins)
    ampl_hist_win, _ = np.histogram(allevents[allevents['winsize'] == winsize]['amplitude'], bins = ampl_bins)
        
    # write to dataframe
    score_hist.loc[:, winsize] = scores_hist_win
    ampl_hist.loc[:, winsize] = ampl_hist_win
    

# %% figure

# init plotting
from functions.initialize_plotting import *
winsize_colors = {36: '#003f5c', 96 : '#7a5195', 114 : '#ef5675', 276 : '#ffa600'}

fig, ax = plt.subplots(nrows = 1, ncols = 1,
                       figsize = get_figure_size(width = 79.6, height = 70),
                       dpi = 300,
                       layout = 'constrained')


for winsize, txt_o in zip([96, 114], [0, 1]):
    
    # get data
    values = score_hist[winsize].to_list()
    edges = PSC_bins['scores']
    
    ax.stairs(values, edges,
              lw = 1.,
              color = winsize_colors[winsize],
              label = str(winsize) + ' ms')
    
for i_sbin, sbin in enumerate(score_bins):
    if sbin >= 0.7 and sbin < 1.:
        
        v_w1 = score_hist.at[sbin, winsizes[0]]
        v_w2 = score_hist.at[sbin, winsizes[1]]
        
        xs = [sbin+0.005]*2
        
        if v_w1 < v_w2:
            ys = [y+60 for y in [v_w1, v_w2]]
            ts = [v_w1, v_w2]
            cs = [winsize_colors[winsizes[0]], winsize_colors[winsizes[1]]]
        else:
            ys = [y+60 for y in [v_w2, v_w1]]
            ts = [v_w2, v_w1]
            cs = [winsize_colors[winsizes[1]], winsize_colors[winsizes[0]]]
        
        for wi, winsize in enumerate(winsizes):
            ax.text(x = xs[wi],
                    y = ys[wi] + (wi*100), #values[i_sbin] + 60 + (txt_o * 110),
                    s = ts[wi],
                    color = cs[wi],
                    fontsize = 6,
                    rotation = 90,
                    ha = 'center', va = 'bottom')

    
ax.legend(loc = 'upper left', frameon = False, fontsize = 7,
          title = 'winsize', title_fontsize = 7)
    
# x
apply_axis_settings(ax, axis = 'x', ax_min=0.7, ax_max=1.0, pad=None, step=0.1, stepminor=0.02, label='Event score')

# y
apply_axis_settings(ax, axis = 'y', ax_min=0, ax_max=1800, pad=None, step=500, stepminor=100, label='Event count [#]')


# align labels
fig.align_labels()

# remove spiness
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# display figure
plt.show()


# %%

fig, ax = plt.subplots(nrows = 1, ncols = 1,
                       figsize = get_figure_size(width = 79.6, height = 70),
                       dpi = 300,
                       layout = 'constrained')


for winsize in winsizes:
    
    # get data
    values = ampl_hist[winsize].to_list()
    edges = PSC_bins['amplitudes']
    
    ax.stairs(values, edges,
              lw = 1.,
              color = winsize_colors[winsize],
              label = str(winsize) + ' ms')
    


    
ax.legend(loc = 'upper left', frameon = False, fontsize = 7,
          title = 'winsize', title_fontsize = 7)
    
# x
apply_axis_settings(ax, axis = 'x', ax_min=-50, ax_max=0, pad=None, step=10, stepminor=5, label='Event amplitude [pA]')

# y
apply_axis_settings(ax, axis = 'y', ax_min=0, ax_max=600, pad=None, step=100, stepminor=20, label='Event count [#]')


# align labels
fig.align_labels()

# remove spiness
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# display figure
plt.show()