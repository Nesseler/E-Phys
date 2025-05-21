# -*- coding: utf-8 -*-
"""
Created on Fri May 16 16:01:33 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, synaptic_dir
from parameters.parameters import PSC_bins

# set threshold
th = 0.5

create = False
load = not create

labels = True

# %% load or create

if create: 
    from miniML_dtc.miniML_validation_create_histdfs import *
    
# load excel sheet
elif load:
    
    # create dict for dataframes
    hist_dfs = dict.fromkeys(['score', 'amplitude', 'halfdecaytime', 'risetime'])
    
    for k in hist_dfs.keys():
        hist_dfs[k] = pd.read_excel(synaptic_dir + '/miniML_validation' + f'/{k}_hist_{str(th).replace(".", "p")}.xlsx', 
                                   index_col = f'{k}_bins')


# %% figure

# set window sizes for plotting
winsizes = [96, 114] #[36, 96, 114, 276]

# init plotting
from functions.initialize_plotting import *
winsize_colors = {36: '#003f5c', 96 : '#7a5195', 114 : '#ef5675', 276 : '#ffa600'}

# iterate through measures
for measure in hist_dfs.keys():
    
    # get data
    hist = hist_dfs[measure]
    bins = PSC_bins[measure]

    fig, ax = plt.subplots(nrows = 1, ncols = 1,
                           figsize = get_figure_size(width = 79.6, height = 70),
                           dpi = 300,
                           layout = 'constrained')
    
    
    for txt_o, winsize in enumerate(winsizes):
        
        # get data
        values = hist[winsize].to_list()
        edges = bins
        
        ax.stairs(values, edges,
                  lw = 1.,
                  color = winsize_colors[winsize],
                  label = str(winsize) + ' ms')
        
    if labels:
        
        # get necessary values
        bin_steps = np.diff(bins)[0]
        
        # dict for text labels (limit 0, limit 1, offset in y, spacing_factor)
        labelrange = {'score' : [0.7, 1., 60, 10], 
                      'amplitude' : [-30, 0, 20, 3.4],
                      'halfdecaytime' : [0, 20, 8, 1.25],
                      'risetime' : []}
    
    if len(labelrange[measure]) == 4:
        for i_sbin, sbin in enumerate(bins):
            if sbin >= labelrange[measure][0] and sbin < labelrange[measure][1]:
                v_w1 = hist.at[sbin, winsizes[0]]
                v_w2 = hist.at[sbin, winsizes[1]]
                
                xs = [sbin+(bin_steps/1.8)]*2
                
                if v_w1 < v_w2:
                    ys = [y+labelrange[measure][2] for y in [v_w2, v_w2]]
                    ts = [v_w1, v_w2]
                    cs = [winsize_colors[winsizes[0]], winsize_colors[winsizes[1]]]
                else:
                    ys = [y+labelrange[measure][2] for y in [v_w1, v_w1]]
                    ts = [v_w2, v_w1]
                    cs = [winsize_colors[winsizes[1]], winsize_colors[winsizes[0]]]
                
                for wi, winsize in enumerate(winsizes):
                    
                    if ys[wi] < 100:
                        spacing_f = 10 * labelrange[measure][3]
                    elif ys[wi] > 100 and ys[wi] < 1000:
                        spacing_f = 12 * labelrange[measure][3]
                    elif ys[wi] > 1000:
                        spacing_f = 14 * labelrange[measure][3]
                    
                    ax.text(x = xs[wi],
                            y = ys[wi] + (wi*spacing_f), #values[i_sbin] + 60 + (txt_o * 110),
                            s = ts[wi],
                            color = cs[wi],
                            fontsize = 5,
                            rotation = 90,
                            ha = 'center', va = 'bottom')
    
    if measure == 'score':
        # score
        apply_axis_settings(ax, axis = 'x', ax_min=0.7, ax_max=1.0, pad=None, step=0.1, stepminor=0.02, label='Event score')
        apply_axis_settings(ax, axis = 'y', ax_min=0, ax_max=1900, pad=None, step=500, stepminor=100, label='Event count [#]')
        legendloc = 'upper left'
        
    elif measure == 'amplitude':
        apply_axis_settings(ax, axis = 'x', ax_min=-30, ax_max=0, pad=None, step=10, stepminor=1, label='Event amplitude [pA]')
        apply_axis_settings(ax, axis = 'y', ax_min=0, ax_max=600, pad=None, step=100, stepminor=20, label='Event count [#]')
        legendloc = 'upper left'

    elif measure == 'halfdecaytime':
        apply_axis_settings(ax, axis = 'x', ax_min=0, ax_max=20, pad=None, step=5, stepminor=1, label='Half decay time [ms]')
        apply_axis_settings(ax, axis = 'y', ax_min=0, ax_max=250, pad=None, step=50, stepminor=5, label='Event count [#]')
        legendloc = 'upper right'
        
    elif measure == 'risetime':
        apply_axis_settings(ax, axis = 'x', ax_min=0, ax_max=15, pad=None, step=5, stepminor=1, label='Rise time [ms]')
        apply_axis_settings(ax, axis = 'y', ax_min=0, ax_max=300, pad=None, step=50, stepminor=5, label='Event count [#]')
        legendloc = 'upper right'
        
    # legend
    ax.legend(loc = legendloc, frameon = False, fontsize = 7,
              title = 'winsize', title_fontsize = 7)
    
    # align labels
    fig.align_labels()
    
    # display figure
    plt.show()
    
    # save figure 
    save_figures(fig, 
                 figure_name = f'figure-miniML_validation-histo-{measure}', 
                 save_dir = figure_dir + '/miniML_validation',
                 figure_format = 'both')


