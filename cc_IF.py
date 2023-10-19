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
from useful_functions import calc_time_series, butter_filter, save_figures
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


test_step = potential_df[12] * 1e3

# test_AP = test_step[4900:5500]

test_AP_filtered = butter_filter(test_step, order=1, cutoff=1e3, sampling_rate=20e3)



v = test_AP_filtered[5300:5450]
t = calc_time_series(v, 20e3, 'ms')
dvdt = calc_dvdt(v, t)



v_derivatives, dv_axs = plt.subplots(2,1, layout = 'constrained', sharex=True)
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



# #plot_voltage_v_time(d2vdt[:-2], dv_axs[2], {'c':'r', 'label':'d2v'}, v_range=[-500, 500])

# lc = get_colorcode(t[:-2], d2vdt, dvdt, norm = norm, cmap='seismic')
# dv_axs[2].set_ylim([-500, 500])

# line = dv_axs[2].add_collection(lc)
# v_derivatives.colorbar(line, ax=dv_axs[2])

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

# ToDo
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
    
    
    
  
test_step = potential_df[12] * 1e3
test_AP = test_step
test_AP_filtered = butter_filter(test_AP, order=1, cutoff=1e3, sampling_rate=20e3)

#v = test_AP_filtered


v = test_AP_filtered[5300:5450]
t = calc_time_series(v, 20e3, 'ms')
dvdt = calc_dvdt(v, t)



v_derivatives, dv_axs = plt.subplots(2,1, layout = 'constrained', sharex=True)

norm = mtl.colors.CenteredNorm(0, dvdt.max())
lc = get_colorcode(t, v, dvdt, norm = norm, cmap='seismic')
line_v = dv_axs[0].add_collection(lc)
dv_axs[0].set_ylim([-100, 20])
v_derivatives.colorbar(line_v, ax=dv_axs[0])

lc = get_colorcode(t[:-1], dvdt, dvdt, norm = norm, cmap='seismic')
dv_axs[1].set_ylim([-150, 250])
line = dv_axs[1].add_collection(lc)

v_derivatives.colorbar(line, ax=dv_axs[1])

dv_axs[1].set_xlim([t.min(),t.max()])


idx_peaks, peak_dict = sc.signal.find_peaks(v, prominence = 20)
APs_parameters = get_AP_parameters(v, idx_peaks)


import seaborn as sbn

def plot_AP_parameter_v_index(AP_parameter, e_range=None, label=''):
    fig_AP, AP_axs = plt.subplots(1,1, layout = 'constrained')
    
    #sbn.violinplot(ax=AP_axs, data=AP_parameter)
    sbn.swarmplot(ax=AP_axs, data=AP_parameter)
    
    if e_range is not None:
        AP_axs.set_ylim(e_range)
        
    AP_axs.set_ylabel(label)
    


    
    
    
    
    
    
# %%

def get_step_parameters(v, i, SR=20e3, prominence = 40, dvdt_threshold=25):
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
            instant_freq
            ISI_adapt
            v_amplitude_adapt
    """
     
    ###PEAK IDX
    idx_peaks, peak_dict = sc.signal.find_peaks(v, prominence=prominence)
    
    AP_df = get_AP_parameters(v, idx_peaks, dvdt_threshold=dvdt_threshold)
    
    
    ### NUMBER OF PEAKS
    n_peaks = len(idx_peaks)
    
    ### ISIs
    if len(idx_peaks) > 1: 
        ISIs = np.diff(idx_peaks) / (SR/1e3)
        instant_freqs = 1 / (ISIs / 1e3)
        if len(idx_peaks) > 3:
            instant_freq = np.mean(instant_freqs[:3])
        else:
            instant_freq = np.nan
    else:
        ISIs = np.nan
        instant_freqs = np.nan
        instant_freq = np.nan
    
    # Spike amplitude
    v_amplitude = AP_df['v_amplitude'].to_numpy()
    
    ###ADAPTATION
    
    # Adaptation factor ISI
    if len(idx_peaks) > 1: 
        ISI_adapt = ISIs[-1] / ISIs[0]
        
        # Adaptation factor spike amplitude
        v_adapt = v_amplitude[-1] / v_amplitude[0]
    else:
        ISI_adapt = np.nan
        v_adapt = np.nan
        
    ###RETURN DICtIONARY
    step_parameters = {'n_peaks' : n_peaks,
                       'AP_parameters' : AP_df,
                       'idx_peaks' : idx_peaks,
                       'ISIs' : ISIs,
                       'instant_freq' : instant_freq,
                       'instant_freqs' : instant_freqs,
                       'ISI_adapt' : ISI_adapt,
                       'v_adapt' : v_adapt}

    return step_parameters




# %%

test_step_index = 25

test_step = potential_df[test_step_index] * 1e3
test_step_current = current_df[test_step_index] * 1e12

test_step_filtered = butter_filter(test_step, order=3, cutoff=1e3, sampling_rate=20e3)

v = test_step_filtered
i = (np.mean(current_df[i].iloc[15]) - np.mean(current_df[i].iloc[15]))*1e12
t = calc_time_series(v, 20e3, 'ms')


step_Ps = get_step_parameters(v, i)

darkmode_bool = True

if darkmode_bool:
    plt.style.use('dark_background')
    prime_color = 'w'
    color1 = 'lightblue'
    color2 = 'magenta'
    color3 = 'red'
    cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","white","magenta"])
    plt.rcParams['axes.grid'] = False
    plot_dict = {'color':prime_color, 'linewidth' : 0.5}
elif darkmode_bool == False:
    plt.style.use('default')
    prime_color = 'k'
    color1 = 'blue'
    color2 = 'purple'
    color3 = 'red'
    cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","grey","red"])
    plt.rcParams['axes.grid'] = True
    plot_dict = {'color':prime_color, 'linewidth' : 0.5}

### figure

fig_one_step = plt.figure()
mm = 1/25.4

fig_one_step = plt.figure(layout='constrained',
                          figsize=(328.67*mm, 165.5*mm))

subfigs = fig_one_step.subfigures(1, 2)

axs_v = subfigs[0].subplots(3,1,gridspec_kw={'height_ratios': [1,1,5]},
                            sharex=True)


### current step

axs_v[0].plot(t,test_step_current, c=color3, linewidth=1)

axs_v[0].set_ylabel('Inj. current\n[pA]')
axs_v[0].set_ylim([-50,350])
axs_v[0].set_yticks(np.arange(0,251,250))
axs_v[0].set_yticks(np.arange(-50,351,50), minor = True)



### eventplot

axs_v[1].eventplot(step_Ps['AP_parameters']['t_peaks'], 
                lineoffsets=0, 
                linelengths=5, 
                color = prime_color, 
                linewidth = 1)

axs_v[1].set_ylabel('APs\n')
#axs_v[0].set_ylim([-50,350])
axs_v[1].set_yticks([])
axs_v[1].tick_params(axis = 'x', size = 0)



### voltage step

axs_v[2].plot(t,v, prime_color, linewidth=1)

axs_v[2].set_xlabel('Time in [ms]')
axs_v[2].set_xlim([0,1500])
axs_v[2].set_xticks(np.arange(0,1501,1500/3))

axs_v[2].set_ylabel('Membrane potential \n[mV]')
axs_v[2].set_ylim([-100,20])
axs_v[2].set_yticks(np.arange(-100,41,20))




subfigs[0].align_ylabels(axs_v[:])

### right side

axs_p = subfigs[1].subplots(3,1, sharex=True)

axs_p[0].tick_params(axis = 'x', size = 0)

# fig_one_step.set_constrained_layout_pads(hspace=0.0, h_pad=0.0)

peak_IDs = np.arange(0,step_Ps['n_peaks'])
ISI_IDs = np.arange(0.5,step_Ps['n_peaks']-1,1)

### Amplitude
axs_p[0].plot(peak_IDs, step_Ps['AP_parameters']['v_amplitude'], prime_color)

axs_p[0].set_ylim([0,100])
axs_p[0].set_ylabel('AP amplitude [mV]')


### rise time
axs_p[1].plot(peak_IDs, step_Ps['AP_parameters']['t_toPeak'], color1)

axs_p[1].set_ylim([0,2])
axs_p[1].set_ylabel('Time to peak\n[ms]')


### ISI

axs_p[2].plot(ISI_IDs, step_Ps['ISIs'], color2)

axs_p[2].set_ylim([0,100])
axs_p[2].set_ylabel('ISI [ms]')

axs_p[2].set_xlim([0,step_Ps['n_peaks']-1])
axs_p[2].set_xlabel('Spike number [#]')



save_figures(fig_one_step, 
              figure_name = 'E-004-example_step', 
              save_dir = 'C:/Users/nesseler/Desktop/DC_figs', 
              darkmode_bool = darkmode_bool)



# %%


test_step = potential_df[35] * 1e3
test_AP_filtered = butter_filter(test_step, order=1, cutoff=1e3, sampling_rate=20e3)


v = test_AP_filtered
idx_peaks, peak_dict = sc.signal.find_peaks(v, prominence = 20)
SR=20e3
dvdt_threshold = 25
t_pre = 2
t_post = 5

t = calc_time_series(v, SR, 'ms')


v_derivatives, dv_axs = plt.subplots(2,1, layout = 'constrained', sharex=True)

norm = mtl.colors.CenteredNorm(0, dvdt.max())
lc = get_colorcode(t, v, dvdt, norm = norm, cmap='seismic')
line_v = dv_axs[0].add_collection(lc)
dv_axs[0].set_ylim([-100, 30])
v_derivatives.colorbar(line_v, ax=dv_axs[0])

lc = get_colorcode(t[:-1], dvdt, dvdt, norm = norm, cmap='seismic')
dv_axs[1].set_ylim([-150, 250])
line = dv_axs[1].add_collection(lc)

v_derivatives.colorbar(line, ax=dv_axs[1])

dv_axs[1].set_xlim([t.min(),t.max()])


dv_axs[0].eventplot(idx_peaks / (SR/1e3), 
                lineoffsets=20, 
                linelengths=5, 
                color = 'k', 
                linewidth = 1)


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


    
  










import seaborn as sbn

def plot_AP_parameter_v_index(AP_parameter, e_range=None, label=''):
    fig_AP, AP_axs = plt.subplots(1,1, layout = 'constrained')
    
    #sbn.violinplot(ax=AP_axs, data=AP_parameter)
    sbn.swarmplot(ax=AP_axs, data=AP_parameter)
    
    if e_range is not None:
        AP_axs.set_ylim(e_range)
        
    AP_axs.set_ylabel(label)
    





# %%

#12, 25, 37, 42

test_step_index = 25
darkmode_bool = True

for test_step_index in [12, 25, 37, 42]:
    # for darkmode_bool in [True, False]:
    xlimits = [225, 375]
    # xlimits = [200, 1300]
    
    if darkmode_bool:
        plt.style.use('dark_background')
        prime_color = 'w'
        color1 = 'lightblue'
        color2 = 'magenta'
        color3 = 'red'
        cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","white","magenta"])
        plt.rcParams['axes.grid'] = False
    elif darkmode_bool == False:
        plt.style.use('default')
        prime_color = 'k'
        color1 = 'blue'
        color2 = 'purple'
        color3 = 'red'
        cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","grey","red"])
        plt.rcParams['axes.grid'] = True
    
    
    
    test_step = potential_df[test_step_index] * 1e3
    test_step_current = current_df[test_step_index] * 1e12
    
    test_AP_filtered = butter_filter(test_step, order=3, cutoff=1e3, sampling_rate=20e3)
    
    
    v = test_AP_filtered
    idx_peaks, peak_dict = sc.signal.find_peaks(v, prominence = 20)
    SR=20e3
    dvdt_threshold = 25
    t_pre = 2
    t_post = 5
    
    t = calc_time_series(v, SR, 'ms')
    
    dvdt = calc_dvdt(v, t)
    
    fig_AP_detect, axs_AP_detect = plt.subplots(5,1, layout = 'constrained',
                                                gridspec_kw={'height_ratios': [2,1,5,1,5]},
                                                sharex=True,
                                                figsize=(328.67*mm, 165.5*mm))
    
    ### current step
    
    axs_AP_detect[0].plot(t,test_step_current, c=color3, linewidth=1)
    
    axs_AP_detect[0].set_ylabel('Inj. current\n[pA]')
    axs_AP_detect[0].set_ylim([-50,350])
    axs_AP_detect[0].set_yticks(np.arange(0,251,250))
    axs_AP_detect[0].set_yticks(np.arange(-50,351,50), minor = True)
    
    
    
    
    ### eventplot
    
    axs_AP_detect[1].eventplot(idx_peaks / (SR/1e3), 
                    lineoffsets=0, 
                    linelengths=5, 
                    color = prime_color, 
                    linewidth = 1)
    
    # axs_AP_detect[1].set_ylabel('APs\n[#]')
    #axs_v[0].set_ylim([-50,350])
    axs_AP_detect[1].set_yticks([])
    axs_AP_detect[1].tick_params(axis = 'x', size = 0)
    
    
    
    
    ### voltage step
    
    cmap_max = 50
    
    norm = mtl.colors.CenteredNorm(0, cmap_max)
    lc = get_colorcode(t, v, dvdt, norm = norm, cmap=cmap)
    line_v = axs_AP_detect[2].add_collection(lc)
    axs_AP_detect[2].set_ylim([-100, 30])
    
    cbar1 = fig_AP_detect.colorbar(line_v, ax=axs_AP_detect[2],
                                   label='dvdt [mV/ms]',
                                   ticks=np.arange(-cmap_max,cmap_max+1,cmap_max))
    
    axs_AP_detect[2].set_ylabel('Membrane potential \n[mV]')
    axs_AP_detect[2].set_ylim([-100,20])
    axs_AP_detect[2].set_yticks(np.arange(-100,41,50))
    axs_AP_detect[2].set_yticks(np.arange(-100,31,10), minor = True)
    
    
    
    ### eventplot threshold crossing dvdt
    dvdt_above_th = np.where(dvdt >= dvdt_threshold, 1, 0)
    dvdt_change = np.diff(dvdt_above_th)
    idx_threshold = np.where(dvdt_change == 1)[0] + 1 #+1 for backwards derivative
    t_threshold = np.divide(idx_threshold, (SR / 1e3))
    
    axs_AP_detect[3].eventplot(t_threshold, 
                    lineoffsets=0, 
                    linelengths=5, 
                    color = prime_color, 
                    linewidth = 1)
    
    # axs_AP_detect[3].set_ylabel('Thresholds\n[#]')
    #axs_v[0].set_ylim([-50,350])
    axs_AP_detect[3].set_yticks([])
    axs_AP_detect[3].tick_params(axis = 'x', size = 0)
    
    
    
    
    ### first derivate
    
    lc = get_colorcode(t[:-1], dvdt, dvdt, norm = norm, cmap=cmap)
    axs_AP_detect[4].set_ylim([-150, 250])
    line = axs_AP_detect[4].add_collection(lc)
    
    cbar2 = fig_AP_detect.colorbar(line_v, ax=axs_AP_detect[4],
                                   label='dvdt [mV/ms]',
                                   ticks=np.arange(-cmap_max,cmap_max+1,cmap_max))
        
    axs_AP_detect[4].set_ylabel('Rate of\nmembrane potential\nchange [mV/ms]')
    
    axs_AP_detect[4].set_ylim([-50,50])
    # axs_AP_detect[4].set_yticks(np.arange(-150,251,150))
    # axs_AP_detect[4].set_yticks(np.arange(-150,251,50), minor = True)
    
        
    axs_AP_detect[4].hlines(dvdt_threshold, t.min(), t.max(), linestyle=':', color=prime_color)
    
    fig_AP_detect.align_ylabels(axs_AP_detect)
    
    
    
    axs_AP_detect[4].set_xlim(xlimits)
    # axs_AP_detect[4].set_xticks(np.arange(250, 1250+1,250))
    axs_AP_detect[4].set_xticks(np.arange(225, 375+1,25))
    axs_AP_detect[4].set_xlabel('Time [ms]')
    
    plt.show()
    
    save_figures(fig_AP_detect, 
                figure_name = 'E-004-example_step ' + str(test_step_index) + ' ' + str(xlimits[1]), 
                save_dir = 'C:/Users/nesseler/Desktop/DC_figs', 
                darkmode_bool = darkmode_bool)




# %%

darkmode_bool = True

v_range = [-100, 40]

if darkmode_bool:
    plt.style.use('dark_background')
    prime_color = 'w'
    color1 = 'lightblue'
    color2 = 'magenta'
    color3 = 'red'
    cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","white","magenta"])
    plt.rcParams['axes.grid'] = False
    plot_dict = {'color':prime_color, 'linewidth' : 0.5}
elif darkmode_bool == False:
    plt.style.use('default')
    prime_color = 'k'
    color1 = 'blue'
    color2 = 'purple'
    color3 = 'red'
    cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","grey","red"])
    plt.rcParams['axes.grid'] = True
    plot_dict = {'color':prime_color, 'linewidth' : 0.5}


#import patchview

import matplotlib.pylab as plt
testFile = r"C:\Users\nesseler\Desktop\local E-Phys\CtC-AMY-ctrl-008-2.dat"

## import HekaHelpers which is a wrapper of main file: HEKA_Reader_MAIN
from patchview.HekaIO.HekaHelpers import HekaBundleInfo

## read Heka .dat file into a object
# test file is the full path for your Heka .dat file
bundleTester = HekaBundleInfo(testFile)

## Get sample rate
traceIndex = [0,4,0,0] ## [Group, Series, Sweep, Trace]
#SR = bundleTester.getSeriesSamplingRate(traceIndex)

## Get stimuli information
#time, stim, stimInfo = bundleTester.getStim(traceIndex)

## Get data from a single sweep and single channel
data = bundleTester.getSeriesData(traceIndex)

SR = bundleTester.getSeriesSamplingRate(traceIndex)

v1 = data[:,0,0] * 1e3

mm = 1/25.4

fig1, axs1 = plt.subplots(3, 1, 
                         layout = 'constrained', 
                         sharex=True,
                         figsize=(328.67*mm, 165.5*mm))


t = calc_time_series(v1, SR, 's')

axs1[0].plot(t, v1, **plot_dict)

# axs1[1].set_xlabel(f'Time [s]')
axs1[0].set_xlim([t[0], t[-1]])
axs1[0].set_xticks(np.arange(0,30+1,5))

axs1[0].set_ylim(v_range)
axs1[0].set_ylabel('Cell #1\nMembrane\npotential [mV]')
axs1[0].set_yticks(np.arange(-100,30+1,50))
axs1[0].set_yticks(np.arange(-100,30+1,10), minor=True)







traceIndex = [1,4,0,0] ## [Group, Series, Sweep, Trace]
#SR = bundleTester.getSeriesSamplingRate(traceIndex)

## Get stimuli information
#time, stim, stimInfo = bundleTester.getStim(traceIndex)

## Get data from a single sweep and single channel
data = bundleTester.getSeriesData(traceIndex)

SR = bundleTester.getSeriesSamplingRate(traceIndex)

v2 = data[:,0,0] * 1e3

t = calc_time_series(v2, SR, 's')

axs1[1].plot(t, v2, **plot_dict)

# axs1[1].set_xlabel(f'Time [s]')
axs1[1].set_xlim([t[0], t[-1]])
axs1[1].set_xticks(np.arange(0,30+1,5))

axs1[1].set_ylim(v_range)
axs1[1].set_ylabel('Cell #2\nMembrane\npotential [mV]')
axs1[1].set_yticks(np.arange(-100,30+1,50))
axs1[1].set_yticks(np.arange(-100,30+1,10), minor=True)




testFile = r"C:\Users\nesseler\Desktop\local E-Phys\CtC-AMY-ctrl-007.dat"

## import HekaHelpers which is a wrapper of main file: HEKA_Reader_MAIN
from patchview.HekaIO.HekaHelpers import HekaBundleInfo

## read Heka .dat file into a object
# test file is the full path for your Heka .dat file
bundleTester = HekaBundleInfo(testFile)

## Get sample rate
traceIndex = [2,2,0,0] ## [Group, Series, Sweep, Trace]
#SR = bundleTester.getSeriesSamplingRate(traceIndex)

## Get stimuli information
#time, stim, stimInfo = bundleTester.getStim(traceIndex)

## Get data from a single sweep and single channel
data = bundleTester.getSeriesData(traceIndex)

SR = bundleTester.getSeriesSamplingRate(traceIndex)

v3 = data[:,0,0] * 1e3





t = calc_time_series(v3, SR, 's')

axs1[2].plot(t, v3, **plot_dict)

axs1[2].set_xlabel(f'Time [s]')
axs1[2].set_xlim([t[0], t[-1]])
axs1[2].set_xticks(np.arange(0,30+1,5))

axs1[2].set_ylim(v_range)
axs1[2].set_ylabel('Cell #3\nMembrane\npotential [mV]')
axs1[2].set_yticks(np.arange(-100,30+1,50))
axs1[2].set_yticks(np.arange(-100,30+1,10), minor=True)



save_figures(fig1, 
              figure_name = 'resting membrane potentials', 
              save_dir = 'C:/Users/nesseler/Desktop/DC_figs', 
              darkmode_bool = darkmode_bool)



# %%

import time

darkmode_bool = True

if darkmode_bool:
    plt.style.use('dark_background')
    prime_color = 'w'
    color1 = 'lightblue'
    color2 = 'magenta'
    color3 = 'red'
    cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","white","magenta"])
    plt.rcParams['axes.grid'] = False
    plot_dict = {'color':prime_color, 'linewidth' : 0.5}
elif darkmode_bool == False:
    plt.style.use('default')
    prime_color = 'k'
    color1 = 'blue'
    color2 = 'purple'
    color3 = 'red'
    cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","grey","red"])
    plt.rcParams['axes.grid'] = True
    plot_dict = {'color':prime_color, 'linewidth' : 0.5}


testFile = r"C:\Users\nesseler\Desktop\local E-Phys\CtC-AMY-ctrl-008-2.dat"
bundleTester = HekaBundleInfo(testFile)


cell_ID = 0
## Get sample rate
traceIndex = [cell_ID,5,0,0] ## [Group, Series, Sweep, Trace]
#SR = bundleTester.getSeriesSamplingRate(traceIndex)

## Get stimuli information
#time, stim, stimInfo = bundleTester.getStim(traceIndex)

## Get data from a single sweep and single channel
data = bundleTester.getSeriesData(traceIndex)

SR = bundleTester.getSeriesSamplingRate(traceIndex)

v = data[:,0,:] * 1e3
i = data[:,1,:] * 1e12
t = calc_time_series(v, SR, scale='ms')


fig_IF, axs_IF = plt.subplots(2, 1, 
                              layout = 'constrained',
                              figsize=(328.67*mm, 165.5*mm),
                              gridspec_kw={'height_ratios': [2,8]},
                              sharex='col')


n_steps = np.shape(v)[1]
step_idx = np.arange(0,n_steps)

# axs_IF[1][0].plot(t, v, c = prime_color)


vs = np.transpose(v)
i_s = np.transpose(i)

segs = [np.column_stack([t, vi]) for vi in vs]


# line_segments = mtl.collections.LineCollection(segs, array=step_idx, cmap="plasma", linewidth=0.5)

# axs_IF[1][0].add_collection(line_segments)
#cbar = fig_IF.colorbar(line_segments)




norm = mtl.colors.Normalize(vmin=0, vmax=n_steps)
cmap = mtl.cm.ScalarMappable(norm=norm, cmap="plasma")
cmap.set_array(step_idx)


fig_IF.colorbar(cmap, label = 'Number of steps', ax=axs_IF[1])


n_peaks = []

#example_steps = [0, 10, 25, 60]

example_steps = step_idx

for i, vi in enumerate(vs):
    
    if i in example_steps:
        axs_IF[1].plot(t, vi, c=cmap.to_rgba(i + 1))
        axs_IF[0].plot(t, i_s[i], c=cmap.to_rgba(i + 1))
    
    v_pulse = vi[(250*50):(1250*50)]
    i_pulse = i_s[i][(250*50):(1250*50)]
    
    idx_peaks, dict_peak = sc.signal.find_peaks(v_pulse, prominence = 20, width = 50)
    
    n_peaks.append(len(idx_peaks))
    
    #axs_IF[1].plot(n_peaks, step_idx[:i+1])
    

axs_IF[0].set_ylim([-100,300])

axs_IF[0].set_yticks(np.arange(-50,251,50))
axs_IF[0].set_ylabel('Inj. current\n[pA]')



axs_IF[1].set_ylim([-150,50])
axs_IF[1].set_ylabel('Membrane potential\n[mV]')

axs_IF[1].set_xlabel('Time [ms]')

axs_IF[1].set_xlim([0,1500])
axs_IF[1].set_xticks(np.arange(0,1501,250))



save_figures(fig_IF, 
            figure_name = 'step IF example' + str(cell_ID) + ' ' + str(len(example_steps)), 
            save_dir = 'C:/Users/nesseler/Desktop/DC_figs', 
            darkmode_bool = darkmode_bool)







# potential_df = pd.DataFrame()
# current_df = pd.DataFrame()
# time_in_ms_df = pd.DataFrame()

# for i in range(45):
#     current_df[i] = pd.DataFrame(np.transpose(test['Trace_2_6_' + str(i+1) + '_2'])[1])
#     potential_df[i] = pd.DataFrame(np.transpose(test['Trace_2_6_' + str(i+1) + '_1'])[1])



pulse_df = pd.DataFrame()

mean_ls = []

for i in range(44):
    #sampling rate in Hz
    SR = 50e3
    
    #pre, post, and pulse duration in s
    pre_post_dur = 0.250
    pulse_dur = 1
    
    pre_idx = np.arange(0, int(pre_post_dur * SR))
    pulse_idx = np.arange(pre_idx[-1]+1, int((pre_post_dur + pulse_dur) * SR))
    post_idx = np.arange(pulse_idx[-1]+1, int((pre_post_dur + pulse_dur + pre_post_dur) * SR))
    





peak_times = []
inj_cur = []
n_peaks = []
mean_ISI = []

SR = 50e3

for i in range(n_steps):
    mV = vs[i] * 1e3
    
    idx_peaks, dict_peak = sc.signal.find_peaks(vs[i][pulse_idx], prominence = 40, width = 50)
    
    peak_times.append(idx_peaks / (SR/1e3))

    inj_cur.append((np.mean(i_s[i][pulse_idx]) - np.mean(i_s[i][pre_idx])))

    n_peaks.append(len(idx_peaks))

    if len(idx_peaks) > 1: 
        mean_ISI.append(np.mean(np.diff(idx_peaks/(SR/1e3))))
    else:
        mean_ISI.append(np.nan)




IF_fig, IF_axs = plt.subplots(1,3, sharey='row',
                              gridspec_kw={'width_ratios': [0.7, 0.2, 0.1]},
                              layout = 'constrained',
                              figsize=(328.67*mm, 165.5*mm))


IF_axs[0].eventplot(peak_times, orientation = 'horizontal', lineoffsets=inj_cur, linelengths=5, color = "w")
IF_axs[0].set_xlim([0, 1000])
IF_axs[0].set_xlabel('Time\n[ms]')


IF_axs[1].plot(n_peaks, inj_cur, color='w')
IF_axs[1].set_xlabel('AP. freq\n[Hz]')
IF_axs[1].set_xlim([0, 100])


IF_axs[2].plot(mean_ISI, inj_cur, color='r')
IF_axs[2].set_xlabel('Mean ISI\n[ms]')
IF_axs[2].set_xlim([0, 100])


IF_axs[0].set_ylim([-50, 250])
IF_fig.supylabel('Injected current [pA]')



save_figures(IF_fig, 
            figure_name = 'step IF overview' + str(cell_ID), 
            save_dir = 'C:/Users/nesseler/Desktop/DC_figs', 
            darkmode_bool = darkmode_bool)


