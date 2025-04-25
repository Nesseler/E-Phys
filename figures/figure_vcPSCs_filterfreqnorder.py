# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 14:37:18 2025

@author: nesseler
"""

# import standard packages
from functions.initialize_packages import *

# import functions
from functions.functions_import import get_traceIndex_n_file, get_vc_data, get_PSCs_steps
from functions.functions_filter import bessel_filter

# import parameters
from parameters.PGFs import vc_Erest_parameters

SR = vc_Erest_parameters['SR']
t = vc_Erest_parameters['t']
n_steps = vc_Erest_parameters['n_steps']
PGF = 'vc_Erest_3min'


# %% load & filter data

cell_ID = 'E-301'
condition = 'ctrl'

# get the traceIndex and the file path string for data import functions
traceIndex, file_path = get_traceIndex_n_file(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')

# get list of steps to include
idc_steps = get_PSCs_steps(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')

# get data with file path & trace index
i, _, _, _, nsteps_loaded = get_vc_data(file_path, traceIndex, scale='ms')

# concatenate steps
i = i.flatten('C')

# limit dataset
i = i[:3000000]
t = t[:3000000]

# filter
orders = [5, 4, 3, 2, 1, 0]
freqs = [2000, 1000, 750, 500, 250]

# single event index
start = 4.365 #s
winsize = 0.020 #s
event_idc = np.arange(start*SR, (start+winsize)*SR, dtype = int)

# limit time
t_event = np.arange(0, winsize * 1e3, 1/(SR/1e3))
i_unfiltered = i[event_idc]

# set dict
i_dict = dict.fromkeys(freqs)

# filter freqs
for freq in freqs:
    
    # set 2nd dict
    i_dict[freq] = dict.fromkeys(orders)
    
    for order in orders:
        i_filtered = bessel_filter(i, order = order, cutoff = freq, sampling_rate = SR)
        
        # write to dict
        i_dict[freq][order] = i_filtered[event_idc]
        
# add unfiltered 0th order
# i_dict['unfiltered'][0] = i_unfiltered

# %% plotting

# init plotting
from functions.initialize_plotting import *

# init figure
fig, axs = plt.subplots(nrows = len(orders),
                        ncols = len(freqs),
                        dpi = 600,
                        layout = 'constrained',
                        sharex = True,
                        sharey = True)

# loop through frequencies and orders
for freq_i, freq in enumerate(freqs):
    for order_i, order in enumerate(orders):
        axs[order_i][freq_i].plot(t_event, i_dict[freq][order],
                                  lw = 0.5,
                                  color = colors_dict['primecolor'],
                                  zorder = 2)


# edit axis
for freq_i, freq in enumerate(freqs):
    for order_i, order in enumerate(orders):
        
        # set axis
        ax = axs[order_i][freq_i]
        
        # plot title
        ax.set_title(f'{freq} Hz - Order {order}',
                     loc = 'left',
                     size = 7,
                     y = 0.86, x = 0.02,
                     va = 'top')
        
        # plot zero line
        ax.hlines(xmin = 0, xmax = 20, y = 0,
                  lw = 0.5,
                  ls = 'dashed',
                  color = 'gray',
                  zorder = 1)
        
        # plot event onset
        ax.vlines(ymin = -30, ymax = 10, x = 5.5,
                  lw = 0.5,
                  ls = 'dashed',
                  color = 'gray',
                  zorder = 1)
        
        # remove spines    
        [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
        
        # x
        ax.set_xlim([0, 20])
        ax.set_xticks(ticks = np.arange(0, 20+1, 10))
        ax.set_xticks(ticks = np.arange(0, 20+0.5, 5), minor = True)
        
        # y
        ax.set_ylim([-30, 10])
        ax.set_yticks(ticks = np.arange(-30, 10+1, 30))
        ax.set_yticks(ticks = np.arange(-30, 10+0.5, 5), minor = True)

# figure title
fig.suptitle('Filter choice - Bessel filter')

# axis labels
fig.supxlabel('Time [ms]', fontsize = 9)
fig.supylabel('Current [pA]', fontsize = 9)

# align labels
fig.align_labels()

# create saving path and save
from parameters.directories_win import figure_dir
path_fig = join(figure_dir, 'filter_choice')
save_figures(fig, f'filter-single_event-freqs_and_orders', 
             path_fig, 
             darkmode_bool, 
             figure_format='both')


plt.show()



# %% 

# Applying FFT
fft_result = np.fft.fft(i)
freq = np.fft.fftfreq(t.shape[-1], d=1/SR)


# Plotting the spectrum
plt.plot(freq, np.abs(fft_result))
# plt.title('FFT of a Musical Note')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.ylim([0, 0.25*1e7])
plt.xlim([0, 1000])
plt.show()

# Apply Welch's method to estimate the PSD
frequencies, psd_values = sc.signal.welch(i, SR, nperseg=4096)


# Plotting the estimated PSD
plt.semilogy(frequencies, psd_values)
plt.title('Power Spectral Density (PSD) Estimate using Welch\'s Method')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (V^2/Hz)')
plt.xlim([0, 10000])
plt.tight_layout()
plt.show()
