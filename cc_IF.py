# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 19:20:26 2023

@author: nesseler
"""

import pandas as pd
import scipy as sc
import matplotlib as mtl
import matplotlib.pyplot as plt
import numpy as np
from useful_functions import calc_time_series, butter_filter
import time

test = sc.io.loadmat('C:/Users/nesseler/Desktop/local E-Phys/AMY-20230905-E004-cc_IF_06.mat')


potential_df = pd.DataFrame()
current_df = pd.DataFrame()
time_in_ms_df = pd.DataFrame()

for i in range(45):
    current_df[i] = pd.DataFrame(np.transpose(test['Trace_2_6_' + str(i+1) + '_2'])[1])
    potential_df[i] = pd.DataFrame(np.transpose(test['Trace_2_6_' + str(i+1) + '_1'])[1])
    time_in_ms_df[i] = pd.DataFrame(np.transpose(test['Trace_2_6_' + str(i+1) + '_1'])[0])


# %%


def calc_dvdt(v,t):
    #calculate the first derivate dv/dt
    dv = np.diff(v)
    dt = np.diff(t)
    dvdt = dv / dt
    
    return dvdt



def phase_plane_plot(data, axis=None, plot_dict={'c':'k'}, sampling_rate=20e3):
    """
    Calculate and generate phase plane plot for single action potential
    Time is dealt in ms
    
    Parameters:
        data : One-dimensional array with voltage in mV
        axes : Axis to plot phase plane plot on
        sampling_rate : Sampling rate in Hz. Default is 20 kHz.
    Returns:
        v
        dvdt
    """
    v = data
    t = calc_time_series(v, sampling_rate, scale='ms')

    dvdt = calc_dvdt(v,t)
    
    if axis != None:
        axis.plot(v[1:], dvdt, **plot_dict)
        axis.set_ylabel('Rate of membrane potential change\n[mV/ms]')
        axis.set_xlabel('Membrane potential [mV]')

    return v, dvdt



def plot_voltage_v_time(data, axis, plot_dict, v_range = [-100, 20], sampling_rate = 20e3, scale='ms'):
    '''
    Plots the voltage on the given axis with calculated time series.
    Parameters:
        data : One-dimensional array with voltage in mV
        axes : Axis to plot phase plane plot on
        plot_dict : Plotting dictionary that is passed to the matplotlib plot
                    function to specify the plotted line.
        v_range : Range of voltage values covered. Used for the y axis limits.
                  Default is -100 mV to + 20 mV.
        sampling_rate : Sampling rate in Hz. Default is 20 kHz.
        scale : {'ms', 's'} Time scale in which the time series is being 
                calculated. 'ms' Milliseconds is default.
    '''
    #calculate time series
    t = calc_time_series(data, sampling_rate, scale)

    axis.plot(t, data, **plot_dict)

    axis.set_xlabel(f'Time [{scale}]')
    axis.set_xlim([t[0], t[-1]])

    axis.set_ylim(v_range)
    axis.set_ylabel('Membrane potential [mV]')



test_step = potential_df[21] * 1e3
test_AP = test_step[5000:5401]


phaseplane, ppaxs = plt.subplots(1,2, layout = 'constrained')

plot_voltage_v_time(test_AP, ppaxs[0], {'c':'k', 'label':'unfiltered'})

phase_plane_plot(test_AP, ppaxs[1], {'c':'k', 'label':'unfiltered'})





test_AP_filtered = butter_filter(test_AP, order=1, cutoff=1e3, sampling_rate=20e3)

plot_voltage_v_time(test_AP_filtered, ppaxs[0], {'c':'r', 'label':'filtered'})

phase_plane_plot(test_AP_filtered[20:], ppaxs[1], {'c':'r', 'label':'filtered'})

ppaxs[0].legend(loc = 'upper right')
    

# %%

# filter orders figure


test_step = potential_df[21] * 1e3
test_AP = test_step[5000:5301]


phaseplane, ppaxs = plt.subplots(1,2, layout = 'constrained')

#plot_voltage_v_time(test_AP, ppaxs[0], {'c':'k', 'label':'unfiltered'})

#phase_plane_plot(test_AP, ppaxs[1], {'c':'k', 'label':'unfiltered'})


filter_orders = 5

cmap = plt.get_cmap('plasma', filter_orders+1)

for idx, order in enumerate(np.linspace(0, filter_orders, filter_orders+1)):
      
    test_AP_filtered = butter_filter(test_AP, order=order, cutoff=1e3, sampling_rate=20e3)
    
    plot_voltage_v_time(test_AP_filtered, ppaxs[0], {'c':cmap(idx), 'label':f'filtered ({int(order)} order)'})
    
    phase_plane_plot(test_AP_filtered[20:], ppaxs[1], {'c':cmap(idx), 'label':f'filtered ({int(order)} order)'})


# Normalizer
norm = mtl.colors.Normalize(vmin=0, vmax=filter_orders)
  
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
  

cbar = phaseplane.colorbar(sm, ticks=np.linspace(0.4, filter_orders-0.4, 5+1), label = "Butterworth filter order")
cbar.ax.set_yticklabels(np.linspace(0, filter_orders, 5+1,dtype=np.int16))


# TODO
## new figure with insets at threshold, peak & repolarisation




# %%


def get_colorcode(x, y, data_fc, norm=None, cmap='seismic', plot_dict={'c':'k'}):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    return_bool = False

    if norm is None:
        # Create a continuous norm to map from data points to colors
        norm = plt.Normalize(data_fc.min(), data_fc.max())
        return_bool = True
        
        
    lc = mtl.collections.LineCollection(segments, cmap=cmap, norm=norm)
    # Set the values used for colormapping
    lc.set_array(data_fc)
    
    return (norm, lc) if return_bool else lc


test_step = potential_df[14] * 1e3

test_AP = test_step

test_AP_filtered = butter_filter(test_AP, order=1, cutoff=1e3, sampling_rate=20e3)

v = test_AP_filtered[5150:5401]



def plot_vt_n_dvdtv_colorcoded(v, v_range = [-100, 20], sampling_rate = 20e3, scale='ms', cmap='seismic'):
    t = calc_time_series(v, sampling_rate, scale)
    v, dvdt = phase_plane_plot(v,sampling_rate=sampling_rate)
    data_fc = dvdt
    
    cc_phaseplane, cc_ppaxs = plt.subplots(1,2, layout = 'constrained')
    
    norm = mtl.colors.CenteredNorm(0, data_fc.max())

    lc = get_colorcode(t, v, data_fc, norm = norm, cmap=cmap)
    
    line = cc_ppaxs[0].add_collection(lc)
    cc_phaseplane.colorbar(line, ax=cc_ppaxs[1])
    
    lc = get_colorcode(v[1:], data_fc, data_fc, norm = norm, cmap=cmap)
    
    line = cc_ppaxs[1].add_collection(lc)
    
    cc_ppaxs[0].grid(False)
    cc_ppaxs[0].set_xlim(t.min(), t.max())
    cc_ppaxs[0].set_xlabel(f'Time [{scale}]')
    cc_ppaxs[0].set_ylim(v_range)
    cc_ppaxs[0].set_ylabel('Membrane potential [mV]')
    
    cc_ppaxs[1].set_ylim([-100, 250])
    cc_ppaxs[1].set_xlim(v_range)
    cc_ppaxs[1].set_ylabel('Rate of membrane potential change\n[mV/ms]')
    cc_ppaxs[1].set_xlabel('Membrane potential [mV]')
    
    cc_ppaxs[1].grid(False)
    

plot_vt_n_dvdtv_colorcoded(v, cmap = 'viridis')


# %%


test_step = potential_df[11] * 1e3

test_AP = test_step

test_AP_filtered = butter_filter(test_AP, order=3, cutoff=1e3, sampling_rate=20e3)



v = test_AP_filtered
t = calc_time_series(v, 20e3, 'ms')
dvdt = calc_dvdt(v, t)
d2vdt = calc_dvdt(dvdt, t[:-1])


v_derivatives, dv_axs = plt.subplots(3,1, layout = 'constrained', sharex=True)
v_derivatives.set_constrained_layout_pads(wspace=0.0, w_pad=0.0)



norm = mtl.colors.CenteredNorm(0, dvdt.max())
lc = get_colorcode(t, v, dvdt, norm = norm, cmap='seismic')

line_v = dv_axs[0].add_collection(lc)
v_derivatives.colorbar(line_v, ax=dv_axs[0])

#plot_voltage_v_time(v, dv_axs[0], {'c':'k', 'label':'v', 'marker':'.', 'linewidth':0})


#plot_voltage_v_time(dvdt[:-1], dv_axs[1], {'c':'r', 'label':'dv'}, v_range=[-150, 250])
lc = get_colorcode(t[:-1], dvdt, dvdt, norm = norm, cmap='seismic')
dv_axs[1].set_ylim([-150, 250])



line = dv_axs[1].add_collection(lc)
v_derivatives.colorbar(line, ax=dv_axs[1])



#plot_voltage_v_time(d2vdt[:-2], dv_axs[2], {'c':'r', 'label':'d2v'}, v_range=[-500, 500])

lc = get_colorcode(t[:-2], d2vdt, dvdt, norm = norm, cmap='seismic')
dv_axs[2].set_ylim([-500, 500])

line = dv_axs[2].add_collection(lc)
v_derivatives.colorbar(line, ax=dv_axs[2])

#
#dv_axs[2].set_xlim([6.5, 8.5])
#dv_axs[2].set_xlim([5, 13])

std = np.std(v)
ipeaks = sc.signal.find_peaks(v, prominence = 20)


SR = 20e3
peak_x = ipeaks[0] / (SR / 1e3)

#dv_axs[0].vlines(peak_x, -100, 20, 'r')
#dv_axs[1].vlines(peak_x, -150, 250, 'r')
#dv_axs[2].vlines(peak_x, -500, 500, 'r')


dvdtpeaks = sc.signal.find_peaks(dvdt, threshold = 50)

#thresholding of the dvdt
#https://stackoverflow.com/questions/62745970/identify-when-time-series-passes-through-threshold-both-in-graph-and-table

threshold = 25

dv_axs[1].hlines(threshold, t.min(), t.max(), 'r', ':')

dvdt_above_th = np.where(dvdt >= threshold, 1, 0)
dvdt_change = np.diff(dvdt_above_th)

#+1 for backwards derivative
th_crossing_idx = np.where(dvdt_change == 1)[0]+1

peak_x = th_crossing_idx / (SR / 1e3)
dv_axs[0].vlines(peak_x, -100, 20, 'grey')
dv_axs[1].vlines(peak_x, -150, 250, 'grey')



v_at_threshold = v[th_crossing_idx[0]]
t_at_threshold = th_crossing_idx[0] / (SR / 1e3)

print('Threshold',v_at_threshold)
v_at_peak = v[ipeaks[0]]
v_amplitude = v_at_peak-v_at_threshold
print('Amplitude', v_amplitude)


t_peak = (ipeaks[0] / (SR / 1e3))[0]
dv_axs[0].vlines(t_peak, v_at_threshold, v_at_peak, 'k')


# AP Afterhyperpolarisation (AHP)
peak_idx = int(ipeaks[0])

t_post_peak = 5  # in ms
v_post_peak_idx = int(peak_idx + (t_post_peak * SR/1e3))

dv_axs[0].fill_between([t_peak, t_peak+t_post_peak], [-100, -100], [20, 20], color = 'lightgrey')

v_post_peak = v[peak_idx:v_post_peak_idx]


v_AHP = np.min(v_post_peak)

v_AHP_idx = np.argmin(v_post_peak) + peak_idx

t_AHP = v_AHP_idx / (SR/1e3)

dv_axs[0].vlines(t_AHP, -100, 20, 'grey')

v_AHP_amplitude = v_AHP-v_at_threshold

print('AP afterhyperpolarisation', v_AHP_amplitude)

print('AP time to afterhyperpolarisation', t_AHP-t_peak, 'ms')

dv_axs[0].vlines(t_AHP, v_at_threshold, v_AHP, 'lightblue')




# AP full width at half maximum (FWHM)

HM = v_amplitude / 2

v_at_HM = v_at_threshold + HM




#limit timeframe to look for the FWHM
t_pre = 2
t_post = 5
pre_idx = int(peak_idx - (t_pre * (SR/1e3)))
post_idx = int(peak_idx + (t_post * (SR/1e3)))
v_AP = v[pre_idx:post_idx]


v_AP_interp = np.interp

v_AP_above_HM = np.where(v_AP >= v_at_HM, 1, 0)
v_AP_change = np.diff(v_AP_above_HM)

v_AP_change_idx = np.where(v_AP_change != 0)[0]
v_AP_change = v_AP[v_AP_change_idx]
t_AP_change = np.divide(v_AP_change_idx,(SR/1e3))
t_AP_change = np.add(t_AP_change, t_peak - t_pre)

FWHM = np.diff(t_AP_change)[0]

print('FWHM', FWHM, 'ms')

dv_axs[0].scatter(t_AP_change, v_AP_change)

dv_axs[0].hlines(v_at_HM, t_AP_change[0], t_AP_change[1])

dv_axs[0].hlines(v_at_threshold, t_peak-t_pre, t_peak+t_post, 'grey', ':')



# AP time to peak
t_toPeak = np.subtract(t_peak, t_at_threshold)
print('AP time to peak', t_toPeak, 'ms')


# AP rise time
v_pre_peak = v[pre_idx:peak_idx+1]
v_20perc = v_at_threshold + (v_amplitude * 0.2)
v_80perc = v_at_threshold + (v_amplitude * 0.8)

v_rise = np.where((v_pre_peak > v_20perc) & (v_pre_peak < v_80perc))[0]

t_rise = len(v_rise) / (SR / 1e3)

print('AP rise time', t_rise, 'ms')

# TODO
    # Single AP 
        # AP threshold
        # AP amplitude
        # AP FWHM
        # AP time to peak
        # AP rise time
        # AP afterhyperpolarisation
        # AP time to afterhyperpolarisation
        
    # Spike train
        # number of APs
        # times of APs
        # ISIs
        # idx of APs
        # adaptation factor (amplitude, FWHM, ISI)
            #adaptation index: ratio of first to last ISI
        # 



# %%

def get_AP_parameters(v, idx_peaks, SR=20e3, dvdt_threshold=25, t_pre=2, t_post=5):
    """
    Function calculates all parameters associated with an action potential (AP).
    Parameters:
        v : One-dimensional array with voltage in mV.
        peak_idx : One-dimensional array of peak indices.
        SR : Sampling rate in Hz. Default is 20 kHz.
        dvdt_threshold : Threshold in first derivative to calculate threshold
            crossing of the AP (in ms/mV). Default is 25 ms/mV.
        t_pre : Time before peak to investigate in ms. Default is 2 ms.
        t_post : Time after peak to investigate in ms. Default is 5 ms.
    Returns:
        AP_parameters: Pandas Dataframe of all parameters for the provided peaks.
            v_peaks    
            t_peaks
            v_threshold
            t_threshold
            idx_threshold
            v_amplitude
            t_toPeak
            v_AHP
            t_AHP
            idx_AHP
            v_AHP_amplitude
            t_to_AHP
            FWHM        
    """
    
    t = calc_time_series(v, SR, 'ms')
    
    dvdt = calc_dvdt(v, t)
    
    d2vdt = calc_dvdt(dvdt, t[:-1])
    
    t_peaks = np.divide(idx_peaks, (SR/1e3))
    v_peaks = v[idx_peaks]
    
    ### AP THRESHOLD & AMPLITUDE
    # thresholding of the dvdt
    # https://stackoverflow.com/questions/62745970/identify-when-time-series-passes-through-threshold-both-in-graph-and-table
    # creates array like v with ones and zeros for below and above threshold
    dvdt_above_th = np.where(dvdt >= dvdt_threshold, 1, 0)
    
    # caculates the different from element to element so that only the changes remain
    # with the sign indicating the directionality
    dvdt_change = np.diff(dvdt_above_th)
    
    # AP threshold parameters
    idx_threshold = np.where(dvdt_change == 1)[0] + 1 #+1 for backwards derivative
    v_threshold = v[idx_threshold]
    t_threshold = np.divide(idx_threshold, (SR / 1e3))
    
    # AP amplitude
    v_amplitude = v_peaks-v_threshold
    
    # AP time to peak
    t_toPeak = np.subtract(t_peaks, t_threshold)
    
    
    ### AP AFTERHYPERPOLRISATION (AHP)
    v_AHP = np.zeros_like(idx_peaks, dtype=float)
    t_AHP = np.zeros_like(idx_peaks, dtype=float)
    idx_AHP  = np.zeros_like(idx_peaks, dtype=float)
    v_AHP_amplitude = np.zeros_like(idx_peaks, dtype=float)
    t_to_AHP = np.zeros_like(idx_peaks, dtype=float)
    
    for idx, i_peak in enumerate(idx_peaks):
        #limit v array to time post one peak, since looking for local min
        i_post = int(i_peak + (t_post * SR/1e3))
    
        v_post = v[i_peak:i_post]
        
        #voltage minimum
        v_AHP[idx] = np.min(v_post)
    
        #index
        idx_AHP[idx] = np.argmin(v_post) + i_peak
        
        #time
        t_AHP[idx] = idx_AHP[idx] / (SR/1e3)
    
        #AHP amplitude
        v_AHP_amplitude[idx] = v_AHP[idx] - v_threshold[idx]
    
        #time to afterhyperpolarisation
        t_to_AHP = t_AHP[idx] - t_peaks[idx]
    
    
    ### AP FULL WIDTH AT HALF MAXIMUM (FWHM)
    FWHM = np.zeros_like(idx_peaks, dtype=float)
    
    HM = v_amplitude / 2
    v_HM = v_threshold + HM
    
    for idx, i_peak in enumerate(idx_peaks):
        #limit timeframe to look for the FWHM
        pre_idx = int(i_peak - (t_pre * (SR/1e3)))
        post_idx = int(i_peak + (t_post * (SR/1e3)))
        v_AP = v[pre_idx:post_idx+1]
    
        #use thresholding to find data points above half maximum
        v_AP_above_HM = np.where(v_AP >= v_HM[idx], 1, 0)
        v_change = np.diff(v_AP_above_HM)
    
        idx_change = np.where(v_change != 0)[0]
    
        # v_change = v_AP[idx_change]
        # t_change = np.divide(idx_change, (SR/1e3))
        
        FWHM[idx] = np.diff(idx_change) / (SR/1e3)
    
    
    ### AP RISE TIME
    t_rise = np.zeros_like(idx_peaks, dtype=float)
    
    # rise time is calculated between 20 % and 80 % of the voltage amplitude.
    v_20perc = v_threshold + (v_amplitude * 0.2)
    v_80perc = v_threshold + (v_amplitude * 0.8)
    
    for idx, i_peak in enumerate(idx_peaks):
        #limit v array to time before peak
        pre_idx = int(i_peak - (t_pre * (SR/1e3)))
        v_pre = v[pre_idx:i_peak+1]
        
        
        v_rise = np.where((v_pre > v_20perc[idx]) & (v_pre < v_80perc[idx]))[0]
        
        t_rise[idx] = len(v_rise) / (SR / 1e3)
        
    
    APs_dataframe = pd.DataFrame({'v_peaks' : v_peaks,
                                 't_peaks' : t_peaks,
                                 'v_threshold' : v_threshold,
                                 't_threshold' : t_threshold,
                                 'idx_threshold' : idx_threshold,
                                 'v_amplitude' : v_amplitude,
                                 't_toPeak' : t_toPeak,
                                 'v_AHP' : v_AHP,
                                 't_AHP' : t_AHP,
                                 'idx_AHP' : idx_AHP,
                                 'v_AHP_amplitude' : v_AHP_amplitude,
                                 't_to_AHP' : t_to_AHP,
                                 't_rise' : t_rise,
                                 'FWHM' : FWHM})

    return APs_dataframe
    
    
    
  
test_step = potential_df[11] * 1e3
test_AP = test_step
test_AP_filtered = butter_filter(test_AP, order=3, cutoff=1e3, sampling_rate=20e3)

v = test_AP_filtered
idx_peaks, peak_dict = sc.signal.find_peaks(v, prominence = 20)


APs_parameters = get_AP_parameters(v, idx_peaks)

plt.scatter(x=np.arange(len(idx_peaks)), y=APs_parameters['FWHM'])

import seaborn as sbn

def plot_AP_parameter_v_index(AP_parameter, e_range=None, label=''):
    fig_AP, AP_axs = plt.subplots(1,1, layout = 'constrained')
    
    #sbn.violinplot(ax=AP_axs, data=AP_parameter)
    sbn.swarmplot(ax=AP_axs, data=AP_parameter)
    
    if e_range is not None:
        AP_axs.set_ylim(e_range)
        
    AP_axs.set_ylabel(label)
    

plot_AP_parameter_v_index(APs_parameters['v_amplitude'])#, [0,1], 'AP rise time [ms]')
    
    
    
    
    
    
# %%

def get_step_parameters(v, i, SR=20e3):
    """
    Function calculates all parameters associated with a step of an IF measurement.
    Parameters:
        v : One-dimensional array with voltage in mV.
        i : Current applied at that step in pA.
        SR : Sampling rate in Hz. Default is 20 kHz.
    Returns:
        step_parameters: Dictionary of all parameters for the provided step.
            n_peaks
            t_peaks
            idx_peaks
            ISI
            instant_freqs
            ISI_adapt
            v_amplitude_adapt
    """

    # Spike train
        # number of APs
        # times of APs
        # ISIs
        # idx of APs
        # Adaptation
        
        #factor (amplitude, FWHM, ISI)
            #adaptation index: ratio of first to last ISI
        # 
    print('test')

test_step = potential_df[26] * 1e3

test_step_filtered = butter_filter(test_step, order=3, cutoff=1e3, sampling_rate=20e3)



v = test_step_filtered
i = (np.mean(current_df[i].iloc[15]) - np.mean(current_df[i].iloc[15]))*1e12

t = calc_time_series(v, 20e3, 'ms')

plt.plot(t, v)

###PEAK IDX
idx_peaks, peak_dict = sc.signal.find_peaks(v, prominence = 40)

AP_df = get_AP_parameters(v, idx_peaks, dvdt_threshold=25)

plt.eventplot(AP_df['t_peaks'], lineoffsets=20, linelengths=5, color = "k")

### NUMBER OF PEAKS
n_peaks = len(idx_peaks)

### TIMES OF APs
t_peaks = AP_df['t_peaks'].to_numpy()

### ISIs
if len(idx_peaks) > 1: 
    ISIs = np.diff(idx_peaks) / (SR/1e3)
    instant_freqs = np.diff(idx_peaks) / (SR)
else:
    ISIs = np.nan


###ADAPTATION

# Adaptation factor ISI
ISI_adapt = ISIs[-1] / ISIs[0]

# Adaptation factor spike amplitude
v_amplitude = AP_df['v_amplitude'].to_numpy()

v_adapt = v_amplitude[-1] / v_amplitude[0]













