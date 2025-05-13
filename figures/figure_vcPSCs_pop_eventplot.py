# -*- coding: utf-8 -*-
"""
Created on Fri May  9 13:41:22 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, synaptic_dir

# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol

ehold = 'Erest'
treatment = 'ctrl'
PGF = f'vc-{ehold}-3min'

# get all cell_IDs
cell_IDs = get_cell_IDs_one_protocol(f'vc-{ehold}-3min-{treatment}', 'PGFs_Syn')


# %% set dataframe / dicts

event_peak_locs = dict.fromkeys(cell_IDs)
amplitudes = dict.fromkeys(cell_IDs)
halfdecays = dict.fromkeys(cell_IDs)
charges = dict.fromkeys(cell_IDs)
risetimes = dict.fromkeys(cell_IDs)

# %% load

for cell_ID in tqdm(cell_IDs):
    # set filename
    filename = f'miniMLdetect_{cell_ID}_{ehold}_3min_{treatment}' 
    
    # open a file, where you stored the pickled data
    file = open((synaptic_dir + f'/miniML_dtc-{ehold}-{treatment}/' + filename + '.pickle'), 'rb')
    
    # dump information to that file
    detection = pickle.load(file)
    
    # close and remove (from memory) the file
    file.close()
    del file 
    gc.collect()
    
    # write to dicts
    amplitudes[cell_ID] = detection['individual_values']['amplitudes']
    charges[cell_ID] = detection['individual_values']['charges']
    risetimes[cell_ID] = detection['individual_values']['risetimes']*1e3
    halfdecays[cell_ID] = detection['individual_values']['half_decaytimes'] *1e3
    event_peak_locs[cell_ID] = detection['event_location_parameters']['event_peak_locations'] / 100_000


# %% processing

n_events = pd.DataFrame(index = cell_IDs)
for cell_ID in cell_IDs:
    n_events.at[cell_ID, 'n_events'] = event_peak_locs[cell_ID].shape[0]

cell_IDs = n_events.sort_values(by = ['n_events']).index.to_list()

# 
from parameters.parameters import PSC_bins
ampl_bins = PSC_bins['amplitudes']

ampl_hist = pd.DataFrame(columns = cell_IDs,
                         index = ampl_bins)

# all ampltiudes histogram
# all_amplitudes = np.empty(1)
for cell_ID in cell_IDs:
    ampl_hist.loc[ampl_bins[:-1], cell_ID], _ = np.histogram(amplitudes[cell_ID], bins = ampl_bins)
    # all_amplitudes = np.append(all_amplitudes, amplitudes[cell_ID])

# all_amp_hist, _ = np.histogram(all_amplitudes, bins = amp_bins)


# map all event measurements into one dataframe
event_measures = pd.DataFrame(columns = ['amplitudes', 'charges', 'risetimes', 'half_decaytimes', 'cell_ID'])

for cell_ID in cell_IDs:
    
    cell_eventmeasures = pd.DataFrame.from_dict({'amplitudes' : amplitudes[cell_ID],
                                                 'charges' : charges[cell_ID],
                                                 'risetimes' : risetimes[cell_ID],
                                                 'half_decaytimes' : halfdecays[cell_ID],
                                                 'cell_ID' : [cell_ID] * len(amplitudes[cell_ID])})

    # append to full dataframe
    event_measures = pd.concat([event_measures, cell_eventmeasures], ignore_index=True)
    
# drop events with nan values
event_measures.dropna(axis = 'rows', inplace = True)
    

# %% plot functions

def vcPSCs_xaxis(ax):
    ax.set_xlim([0-1.8, 180+1.8])
    ax.spines['bottom'].set_bounds([0, 180])
    ax.set_xticks(ticks = np.arange(0, 180+1, 30, dtype = int))
    ax.set_xticks(ticks = np.arange(0, 180+1, 5, dtype = int), minor = True)
    ax.set_xlabel('Time [ms]')
    
    
def vcPSCs_ampl_histo_axis(ax):
    ax.set_ylabel('Event count [#]')
    ax.set_ylim([0-2.6, 140+2.6])
    ax.set_yticks(ticks = np.arange(0, 140+1, 20))
    ax.set_yticks(ticks = np.arange(0, 140+1, 5), minor = True)
    ax.spines['left'].set_bounds([0, 140])
    
    
    from parameters.parameters import PSC_bins
    xmin = PSC_bins['amplitudes'][0]
    xpad = abs(xmin)*0.01
    
    ax.set_xlabel('Amplitude [pA]')
    ax.set_xlim([xmin-xpad, 0+xpad])
    ax.spines['bottom'].set_bounds([xmin, 0])   
    ax.set_xticks(ticks = np.arange(xmin, 0+1, 10),
                  labels = np.arange(xmin, 0+1, 10))
    ax.set_xticks(ticks = np.arange(xmin, 0, 1), minor = True)


# %% figure - eventplot

# init plotting
from functions.initialize_plotting import *


fig, ax = plt.subplots(nrows = 1, ncols = 1,
                       figsize = get_figure_size(width = 159.2, height = 80),
                       dpi = 600,
                       layout = 'constrained')

# color coding #
cmap_str = 'inferno_r'
norm = mtl.colors.Normalize(-40, 0)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap=cmap_str)

# title
fig.suptitle(f'Eventplot - {ehold}')

for ci, cell_ID in enumerate(cell_IDs):
    p = event_peak_locs[cell_ID]

    lw = halfdecays[cell_ID]/30
    
    # fixed linelength, fixed color
    ax.eventplot(positions = p, 
                 lineoffsets = ci,
                 linelengths = 0.8,
                 lw = 0.25, 
                 color = 'k')
    
    # linelength proportional to amplitude
    # ax.eventplot(positions = p, 
    #              lineoffsets = [[ci-(amp/2) for amp in np.abs(amplitudes[cell_ID]/60)]],
    #              linelengths = [np.abs(amplitudes[cell_ID])/60],
    #              lw = 0.25, 
    #              color = 'k')
    
    # linelength proportional to charge
    # ax.eventplot(positions = p, 
    #              lineoffsets = [[ci-(ch/2) for ch in np.abs(charges[cell_ID]*2)]],
    #              linelengths = [np.abs(charges[cell_ID]*2)],
    #              lw = 0.25, 
    #              color = 'k')
    
    
# y 
ax.set_ylim([0-0.4-(len(cell_IDs)*0.01), (len(cell_IDs)-0.6)*1.01])
ax.spines['left'].set_bounds([0, len(cell_IDs)-1])
ax.set_ylabel('cells')
ax.set_yticks(ticks = np.arange(0, len(cell_IDs)),
              labels = cell_IDs)

# x
vcPSCs_xaxis(ax)

# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = f'figure-vcPSCs-population-{ehold}-eventplot', 
             save_dir = figure_dir + '/vcPSCs-population',
             figure_format = 'both')


# %% figure - mean amplitude histogram 

fig, ax = plt.subplots(nrows = 1, ncols = 1,
                       figsize = get_figure_size(width = 159.2, height = 80),
                       dpi = 600,
                       layout = 'constrained')

# title
fig.suptitle('Ampltiudes histogram')

mean = ampl_hist.loc[ampl_bins[:-1], :].mean(axis = 1)
std = ampl_hist.loc[ampl_bins[:-1], :].std(axis = 1)

ax.stairs(values = mean.to_list(), 
          edges = ampl_bins,
          fill = False,
          lw = 0.75,
          color = colors_dict['primecolor'],
          label = 'mean')

ax.stairs(values = (mean+std).to_list(), 
          edges = ampl_bins,
          baseline = (mean-std).to_list(),
          fill = True,
          lw = 1,
          alpha = 0.25,
          color = colors_dict['primecolor'],
          label = 'std')

# edit axis
# y
ax.set_ylabel('Event count [#]')
ax.set_ylim([0-1, 100+1])
ax.set_yticks(ticks = np.arange(0, 100+1, 20))
ax.set_yticks(ticks = np.arange(0, 100+1, 5), minor = True)
ax.spines['left'].set_bounds([0, 100])

# x
ax.set_xlabel('Amplitude [pA]')
ax.set_xlim([-70-0.7, 0+0.7])
ax.spines['bottom'].set_bounds([-70, 0])   
ax.set_xticks(ticks = np.arange(-70, 0+1, 10),
              labels = np.arange(-70, 0+1, 10))
ax.set_xticks(ticks = np.arange(-70, 0, 1), minor = True)

# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = f'figure-vcPSCs-population-{ehold}-ampli_hist', 
             save_dir = figure_dir + '/vcPSCs-population',
             figure_format = 'both')


# %% figure - individual amplitude histograms ridge

fig, ax = plt.subplots(nrows = 1, ncols = 1,
                       figsize = get_figure_size(width = 100, height = 120),
                       dpi = 600,
                       layout = 'constrained')

# sort cell_IDs after max occurance
cell_IDs = ampl_hist.max(axis = 0).sort_values().index.to_list()


for ci, cell_ID in enumerate(cell_IDs):
        
    ax.stairs(values = (ampl_hist.loc[ampl_bins[:-1], cell_ID] + (ci*45)).to_list() ,
              edges = ampl_bins+ci/2,
              baseline = (ci*45),
              fill = False,
              lw = 0.5, 
              color = colors_dict['primecolor'])

# y
ax.set_ylabel('Event count [#]')
ax.set_ylim([0-10, 1200])
ax.set_yticks(ticks = np.arange(0, 140+1, 100))
ax.set_yticks(ticks = np.arange(0, 140+1, 10), minor = True)
ax.spines['left'].set_bounds([0, 140])

# x
ax.set_xlabel('Amplitude [pA]')
ax.set_xlim([-70-0.7, len(cell_IDs)/2+0.7])
ax.spines['bottom'].set_bounds([-70, 0])   
ax.set_xticks(ticks = np.arange(-70, 0+1, 10),
              labels = np.arange(-70, 0+1, 10))
ax.set_xticks(ticks = np.arange(-70, 0, 1), minor = True)
    
# remove spines
[ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = f'figure-vcPSCs-population-{ehold}-ampli_hist-indi_ridge', 
             save_dir = figure_dir + '/vcPSCs-population',
             figure_format = 'both')


# %% figure - scatter events 

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
ax.scatter(x = event_measures['amplitudes'],
           y = event_measures['half_decaytimes'],
           marker = '.',
           s = 3,
           alpha = 0.25,
           color = 'k')

# half violin 1
plot_half_violin(event_measures['amplitudes'], 
                 ax = axs[0],
                 v_orientation = 'horizontal',
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1)

# half violin 2
plot_half_violin(event_measures['half_decaytimes'], ax = axs[3],
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1)

# x
apply_axis_settings(axs[0], axis = 'x', ax_min=-70, ax_max=0, pad=None, step=10, stepminor=2, label='')
apply_axis_settings(axs[2], axis = 'x', ax_min=-70, ax_max=0, pad=None, step=10, stepminor=2, label='Amplitude [mV]')

# y
apply_axis_settings(axs[2], axis = 'y', ax_min=0, ax_max=20, pad=None, step=5, stepminor=1, label='Half decay time [ms]')
apply_axis_settings(axs[3], axis = 'y', ax_min=0, ax_max=20, pad=None, step=5, stepminor=1, label='')
    
# remove spines
remove_spines_ticks_labels([axs[0]], axis = 'y')
remove_spines_ticks_labels([axs[3]], axis = 'x')

# remove axis
plt.delaxes(axs[1])

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = f'figure-vcPSCs-population-{ehold}-amp_decay-scatter', 
             save_dir = figure_dir + '/vcPSCs-population',
             figure_format = 'both')


# %% figure - scatter events 

from parameters.parameters import PSC_bins

fig, axs = plt.subplots(nrows = 2, ncols = 2,
                        figsize = get_figure_size(width = 130, height = 100),
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

# 2d histogram 
h2d = ax.hist2d(x = event_measures['amplitudes'],
                 y = event_measures['half_decaytimes'],
                 bins = [PSC_bins['amplitudes'], PSC_bins['halfdecay_times']], 
                 cmax = 50,
                 cmap = 'cividis')

plt.colorbar(h2d[3], ax=ax,
             aspect = 30,
             location='left',
             label = 'Number of events per bin [#]',
             drawedges = False)


# half violin 1
plot_half_violin(event_measures['amplitudes'], 
                 ax = axs[0],
                 v_orientation = 'horizontal',
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1)

# half violin 2
plot_half_violin(event_measures['half_decaytimes'], ax = axs[3],
                 v_direction = 1,
                 v_offset = 0,
                 v_color = colors_dict['primecolor'],
                 v_width = 1)

# x
apply_axis_settings(axs[0], axis = 'x', ax_min=-50, ax_max=0, pad=None, step=10, stepminor=2, label='')
apply_axis_settings(axs[2], axis = 'x', ax_min=-50, ax_max=0, pad=None, step=10, stepminor=2, label='Amplitude [mV]')

# y
apply_axis_settings(axs[2], axis = 'y', ax_min=0, ax_max=20, pad=None, step=5, stepminor=1, label='Half decay time [ms]')
apply_axis_settings(axs[3], axis = 'y', ax_min=0, ax_max=20, pad=None, step=5, stepminor=1, label='')
    
# remove spines
remove_spines_ticks_labels([axs[0]], axis = 'y')
remove_spines_ticks_labels([axs[3]], axis = 'x')

# remove axis
plt.delaxes(axs[1])

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# align labels
fig.align_labels()

# display figure
plt.show()

# save figure 
save_figures(fig, 
             figure_name = f'figure-vcPSCs-population-{ehold}-amp_decay-heat', 
             save_dir = figure_dir + '/vcPSCs-population',
             figure_format = 'both')