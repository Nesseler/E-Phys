# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 19:02:06 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, vplot_dir, synaptic_dir, quant_data_dir

# custom functions
from functions.functions_import import get_vc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_filter import butter_filter

import pickle
import gc

plot = False

PGF = 'vc-Erest-3min'

conditions = ['ctrl', 'adaEk']

SR = 100000

t = np.arange(0, (6*30), 1/SR)

# get all cell_IDs
cell_IDs = get_cell_IDs_one_protocol('vc-Erest-3min-adaEk', 'PGFs_Syn')


# %% set dicts and dataframes

treatments = ['ctrl', 'adaEk']

amplitudes_hist = pd.DataFrame(index = np.arange(-50, 0+1, 1),
                               columns = [f'{cell_ID}-{t}' for cell_ID in cell_IDs for t in treatments])

risetimes_hist = pd.DataFrame(index = np.arange(0, 20+0.2, 0.2),
                              columns = [f'{cell_ID}-{t}' for cell_ID in cell_IDs for t in treatments])

halfdecay_hist = pd.DataFrame(index = np.arange(0, 20+0.5, 0.5),
                              columns = [f'{cell_ID}-{t}' for cell_ID in cell_IDs for t in treatments])

IEI_hist = pd.DataFrame(index = np.arange(0, 30+0.05, 0.05),
                        columns = [f'{cell_ID}-{t}' for cell_ID in cell_IDs for t in treatments])

for df in [amplitudes_hist, risetimes_hist, halfdecay_hist, IEI_hist]:
    df.index.name = 'bins'
    

# %% load event detection

# init plotting
from functions.initialize_plotting import *

# set colors
adaEk_color = {'ctrl' : colors_dict['primecolor'],
               'adaEk' : 'goldenrod'}



for cell_ID in cell_IDs:

    # set dicts
    traces = dict.fromkeys(['ctrl', 'adaEk'])
    amplitudes = dict.fromkeys(['ctrl', 'adaEk'])
    risetimes = dict.fromkeys(['ctrl', 'adaEk'])
    events = dict.fromkeys(['ctrl', 'adaEk'])
    event_locations = dict.fromkeys(['ctrl', 'adaEk'])
    
    
    for treatment in ['ctrl', 'adaEk']:
    
        # set filename
        filename = f'miniMLdetect_{cell_ID}_Erest_3min_{treatment}' 
        
        # open a file, where you stored the pickled data
        file = open((synaptic_dir + f'/miniML_dtc-Erest-{treatment}/' + filename + '.pickle'), 'rb')
        
        # dump information to that file
        detection = pickle.load(file)
        
        # close and remove (from memory) the file
        file.close()
        del file 
        gc.collect()
        
        # write to dicts
        traces[treatment] = detection['mini_trace']
        events[treatment] = detection['events']
        amplitudes[treatment] = detection['individual_values']['amplitudes']
        risetimes[treatment] = detection['individual_values']['risetimes']
        event_locations[treatment] = detection['event_location_parameters']['event_peak_locations']
        
        # calc occurrances
        amplitudes_hist.loc[:amplitudes_hist.index[-2], f'{cell_ID}-{treatment}'], _ = np.histogram(a = detection['individual_values']['amplitudes'], bins = amplitudes_hist.index.to_list())
        risetimes_hist.loc[:risetimes_hist.index[-2], f'{cell_ID}-{treatment}'], _ = np.histogram(a = detection['individual_values']['risetimes'] *1e3, bins = risetimes_hist.index.to_list())
        IEI_hist.loc[:IEI_hist.index[-2], f'{cell_ID}-{treatment}'], _ = np.histogram(a = np.diff(detection['event_location_parameters']['event_peak_locations'] / SR, prepend = 0), bins = IEI_hist.index.to_list())
        halfdecay_hist.loc[:halfdecay_hist.index[-2], f'{cell_ID}-{treatment}'], _ = np.histogram(a = detection['individual_values']['half_decaytimes'] *1e3, bins = halfdecay_hist.index.to_list())

    
    if plot:
        
        # init figure
        fig, axs = plt.subplots(nrows = 5,
                                ncols = 2,  
                                dpi = 300, 
                                layout = 'constrained',
                                figsize = get_figure_size(width = 159.2, height = 170),
                                height_ratios = [3, 0.5, 3, 4.5, 4.5],
                                width_ratios =  [1, 1])
        
        axs = axs.flatten()
        
        # figure title
        fig.suptitle(f'{cell_ID} - PSCs - adapted Ek')
        
        axs_titles = {0 : 'A: PSCs ctrl', 1 : 'B: PSCs adaEk',
                      2 : '', 3 : '',
                      4 : 'C: Events ctrl', 5: 'D: Events adaEk',
                      6 : 'E: Ampltiude', 7 : 'F: IEI',
                      8 : 'G: Risetime', 9 : 'H: Half decay time'}
        
        # traces
        for ti, treatment in enumerate(treatments):
            
            # traces
            ax = axs[ti]
            
            # title
            ax.set_title(axs_titles[ti], loc = 'left', fontsize = 9)
            
            # set x axis
            x = np.arange(0, traces[treatment].shape[0] / SR, 1/SR)
            
            # plot trace
            ax.plot(x, traces[treatment],
                    color= adaEk_color[treatment],
                    lw=0.5,
                    label= 'trace', 
                    zorder = 1)
            
            ax.scatter(x[event_locations[treatment]], traces[treatment][event_locations[treatment]],
                       color = 'r',
                       alpha = 0.5,
                       s = 3,
                       lw=0.5,
                       label='events',
                       zorder = 0)
            
            ax.legend(frameon = False,
                      fontsize = 9,
                      loc = 'lower right',
                      labelspacing = 0.1,
                      borderpad = 0.1,
                      borderaxespad = 0.1)
            
            # y axis
            ax.set_ylim([-40-0.4, 10+0.4])
            ax.spines['left'].set_bounds([-40, 10])
        
                
            # Eventplots
            ax = axs[ti+2]
        
            # plot eventplot
            ax.eventplot(positions = x[event_locations[treatment]],
                          orientation = 'horizontal',
                          lineoffsets = 0,
                          linelengths = 1,
                          linewidths = 0.5,
                          color = adaEk_color[treatment],
                          label = 'events')
        
        axs[0].set_ylabel('Current [pA]')
        axs[0].set_yticks(ticks = np.arange(-40, 10+1, 20))
        axs[0].set_yticks(ticks = np.arange(-40, 10+1, 5), minor = True)
        axs[1].set_yticks(ticks = [])
        
        axs[2].set_ylabel('Events')
        
        for axi in range(4):
            axs[axi].set_xlim([0-1.8, 180+1.8])
            axs[axi].spines['bottom'].set_bounds([0, 180])
            axs[axi].set_xticks(ticks = np.arange(0, 180+1, 60),
                                labels = [])
            axs[axi].set_xticks(ticks = np.arange(0, 180+1, 10), minor = True)
            
        
        remove_spines_n_ticks([axs[0], axs[1]], axis = 'x')
        remove_spines_n_ticks([axs[1], axs[2], axs[3]], axis = 'y')
        
        for axi in [2, 3]:
            axs[axi].set_yticks(ticks = [])
            axs[axi].set_xticks(ticks = np.arange(0, 180+1, 60),
                                labels = np.arange(0, 180+1, 60))
            axs[axi].set_xlabel('Time [ms]')
        
        
        # average events
        
        # traces
        for ti, treatment in enumerate(treatments):
            
            # axis
            axi = ti+4
            ax = axs[axi]
            ax.set_title(axs_titles[axi], loc = 'left', fontsize = 9)
        
            # set data
            event_x = np.arange(0, events[treatment].shape[1]) * 1/(SR/1e3)
            event_average = np.mean( events[treatment], axis=0)
            
            # plot
            ax.plot(event_x, events[treatment].T,
                    c = 'k',
                    alpha = 0.3,
                    lw = 1,
                    label = '_nolegend_')
            
            ax.plot(event_x, event_average,
                    c = 'r', 
                    lw = 1.5,
                    label = 'average event')
            
            ax.legend(frameon = False,
                       fontsize = 9,
                       loc = 'lower right')
            
            # edit axis
            ax.set_ylabel('Current [pA]')
            ax.set_ylim([-41, 10])
            ax.set_yticks(ticks = np.arange(-40, 10+1, 20))
            ax.set_yticks(ticks = np.arange(-40, 10+1, 5), minor = True)
            ax.spines['left'].set_bounds([-40, 10])
            
            ax.set_xlabel('Time [ms]')
            ax.set_xticks(ticks = np.arange(0, round(event_x[-1])+1, 50))
            ax.set_xticks(ticks = np.arange(0, round(event_x[-1])+1, 10), minor = True)
            ax.set_xlim([-5, round(event_x[-1])+5])
            ax.spines['bottom'].set_bounds([0, round(event_x[-1])])
        
        
        # set amplitude histogram axis
        ax = axs[6]
        
        # title
        ax.set_title(axs_titles[6], loc = 'left', fontsize = 9)
        
        ax.stairs(amplitudes_hist[f'{cell_ID}-ctrl'].dropna().to_list(), amplitudes_hist.index.to_list(), 
                  fill = False,
                  lw = 1,
                  alpha = 0.7,
                  color = colors_dict['primecolor'],
                  label = 'ctrl')
        
        # plot treatment
        ax.stairs(amplitudes_hist[f'{cell_ID}-adaEk'].dropna().to_list(), amplitudes_hist.index.to_list(), 
                  fill = False,
                  lw = 1,
                  color = 'goldenrod',
                  label = 'adaEk')
            
        ax.legend(title = 'Treatment', 
                  frameon = False,
                  loc = 'upper left',
                  fontsize = 9)
        
        # edit axis
        ax.set_ylabel('Event count [#]')
        ax.set_ylim([0-2.6, 130+2.6])
        ax.set_yticks(ticks = np.arange(0, 130+1, 20))
        ax.set_yticks(ticks = np.arange(0, 130+1, 5), minor = True)
        ax.spines['left'].set_bounds([0, 130])
        
        ax.set_xlabel('Amplitude [pA]')
        ax.set_xticks(ticks = np.arange(-50, 0+1, 10),
                      labels = np.arange(-50, 0+1, 10))
        ax.set_xticks(ticks = np.arange(-50, 0, 1), minor = True)
        ax.set_xlim([-51, 1])
        ax.spines['bottom'].set_bounds([-50, 0])
        
        
        # set IEI histogram axis
        ax = axs[7]
        
        # title
        ax.set_title(axs_titles[7], loc = 'left', fontsize = 9)
        
        ax.stairs(values = IEI_hist[f'{cell_ID}-adaEk'].dropna().to_list(), 
                  edges = IEI_hist.index.to_list(), 
                  fill = False,
                  lw = 1,
                  alpha = 0.7,
                  color = colors_dict['primecolor'],
                  label = 'ctrl')
        
        # plot treatment
        ax.stairs(values = IEI_hist[f'{cell_ID}-ctrl'].dropna().to_list(), 
                  edges = IEI_hist.index.to_list(), fill = False,
                  lw = 1,
                  color = 'goldenrod',
                  label = 'adaEk')
            
        ax.legend(title = 'Treatment', 
                  frameon = False,
                  loc = 'upper right',
                  fontsize = 9)
        
        # edit axis
        ax.set_ylabel('Event count [#]')
        ax.set_ylim([0-2, 80+2])
        ax.spines['left'].set_bounds([0, 80])
        ax.set_yticks(ticks = np.arange(0, 80+1, 20))
        ax.set_yticks(ticks = np.arange(0, 80+1, 5), minor = True)
        
        ax.set_xlabel('IEI [ms]')
        ax.set_xticks(ticks = np.arange(0, 5+0.1, 1))
        ax.set_xticks(ticks = np.arange(0, 5+0.1, 0.2), minor = True)
        ax.set_xlim([0-0.05, 5+0.05])
        ax.spines['bottom'].set_bounds([0, 5])
        
        
        # set risetimes histogram axis
        ax = axs[8]
        
        # title
        ax.set_title(axs_titles[8], loc = 'left', fontsize = 9)
        
        ax.stairs(values =risetimes_hist[f'{cell_ID}-adaEk'].dropna().to_list(), 
                  edges = risetimes_hist.index.to_list(), 
                  fill = False,
                  lw = 1,
                  alpha = 0.7,
                  color = colors_dict['primecolor'],
                  label = 'ctrl')
        
        # plot treatment
        ax.stairs(values = risetimes_hist[f'{cell_ID}-ctrl'].dropna().to_list(), 
                  edges = risetimes_hist.index.to_list(), fill = False,
                  lw = 1,
                  color = 'goldenrod',
                  label = 'adaEk')
            
        ax.legend(title = 'Treatment', 
                  frameon = False,
                  loc = 'upper right',
                  fontsize = 9)
        
        # edit axis
        ax.set_ylabel('Event count [#]')
        ax.set_ylim([0-2, 80+2])
        ax.spines['left'].set_bounds([0, 80])
        ax.set_yticks(ticks = np.arange(0, 80+1, 20))
        ax.set_yticks(ticks = np.arange(0, 80+1, 5), minor = True)
        
        ax.set_xlabel('Risetime [ms]')
        ax.set_xticks(ticks = np.arange(0, 20+1, 2))
        ax.set_xticks(ticks = np.arange(0, 20+0.1, 0.4), minor = True)
        ax.set_xlim([0-0.2, 20+0.2])
        ax.spines['bottom'].set_bounds([0, 20])
        
        
        # set risetimes histogram axis
        ax = axs[9]
        
        # title
        ax.set_title(axs_titles[9], loc = 'left', fontsize = 9)
        
        ax.stairs(values = halfdecay_hist[f'{cell_ID}-adaEk'].dropna().to_list(), 
                  edges = halfdecay_hist.index.to_list(), 
                  fill = False,
                  lw = 1,
                  alpha = 0.7,
                  color = colors_dict['primecolor'],
                  label = 'ctrl')
        
        # plot treatment
        ax.stairs(values = halfdecay_hist[f'{cell_ID}-ctrl'].dropna().to_list(), 
                  edges = halfdecay_hist.index.to_list(), fill = False,
                  lw = 1,
                  color = 'goldenrod',
                  label = 'adaEk')
            
        ax.legend(title = 'Treatment', 
                  frameon = False,
                  loc = 'upper right',
                  fontsize = 9)
        
        # edit axis
        ax.set_ylabel('Event count [#]')
        ax.set_ylim([0-1, 60+1])
        ax.spines['left'].set_bounds([0, 60])
        ax.set_yticks(ticks = np.arange(0, 60+1, 20))
        ax.set_yticks(ticks = np.arange(0, 60+1, 5), minor = True)
        
        ax.set_xlabel('Half decay times [ms]')
        ax.set_xticks(ticks = np.arange(0, 20+1, 5))
        ax.set_xticks(ticks = np.arange(0, 20+0.1, 1), minor = True)
        ax.set_xlim([0-0.2, 20+0.2])
        ax.spines['bottom'].set_bounds([0, 20])
        
        # align labels
        fig.align_labels()
        
        # remove spines
        [ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]
        
        # display figure
        plt.show()
        
        # create saving path and save
        vplots_path_fig = join(vplot_dir, 'vc_adaEk-PSCs')
        save_figures(fig, f'{cell_ID}-adaEK_PSCs', vplots_path_fig, darkmode_bool, figure_format='png')


# %% save to dataframe

path = quant_data_dir + '/vc-PSCs-adaEk/'

for k, df in {'amplitudes' : amplitudes_hist, 'risetimes' : risetimes_hist, 'halfdecaytimes' : halfdecay_hist, 'IEI' : IEI_hist}.items():
    df.to_excel(path + 'PSCs-adaEk-' + k + '-histogram' + '.xlsx',
                index_label='bins')
    


