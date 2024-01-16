# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 14:31:40 2023

@author: nesseler
"""

import numpy as np
from functions_useful import calc_time_series, calc_dvdt
import pandas as pd
import matplotlib as mtl
import matplotlib.pyplot as plt


def get_AP_parameters(v, idx_peaks, SR=20e3, dvdt_threshold=20, t_pre=2, t_post=5):
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
    
   
        
    
    if len(idx_peaks) > 0:
        t = calc_time_series(v, SR, 'ms')
        
        dvdt = calc_dvdt(v, t)
        
        # d2vdt = calc_dvdt(dvdt, t[:-1])
        
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
        
        # multiple positiv threshold crossing -> choosing idx closest to spike
        if len(idx_threshold) > len(idx_peaks):
            for i_peaks in idx_peaks:
                diff_list = [abs(th-i_peaks) for th in idx_threshold]
                diff_argmin = np.argmin(diff_list)
                idx_threshold = [idx_threshold[diff_argmin]]
        
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
        t1_HM = np.zeros_like(idx_peaks, dtype=float)
        t2_HM = np.zeros_like(idx_peaks, dtype=float)
        
        HM = v_amplitude / 2
        v_HM = v_threshold + HM
        
        
        for idx, i_peak in enumerate(idx_peaks):
            #limit timeframe to look for the FWHM
            pre_idx = int(i_peak - (t_pre * (SR/1e3)))
            
            # fail-save: pre_idx below 0 when pre timeframe exceeds the step
            if pre_idx < 0:
                pre_idx = 0
            
            post_idx = int(i_peak + (t_post * (SR/1e3)))
            v_AP = v[pre_idx:post_idx+1]
        
            #use thresholding to find data points above half maximum
            v_AP_above_HM = np.where(v_AP >= v_HM[idx], 1, 0)
            v_change = np.diff(v_AP_above_HM)
        
            idx_change = np.where(v_change != 0)[0]
        
            # v_change = v_AP[idx_change]
            t_change = np.add(np.divide(idx_change, (SR/1e3)), (i_peak/(SR/1e3)-t_pre)) 

            ### opt: PLOT ###            
            # plt.plot(v)
            # plt.eventplot(idx_peaks, color = 'r', lineoffsets=30, linelengths=5)
            # plt.ylim([-100, 50])
            # plt.hlines(v_HM, 0, 1000)
            # plt.show() 
                
            FWHM[idx] = np.diff(idx_change) / (SR/1e3)
            t1_HM[idx] = t_change[0]
            t2_HM[idx] = t_change[1]
        
        
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
        

        
        dataframe_dict = {'v_peaks' : v_peaks,
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
                          'FWHM' : FWHM,
                          'v_HM' : v_HM,
                          't1_HM' : t1_HM,
                          't2_HM' : t2_HM}
        
        # print(dataframe_dict)
    
    elif len(idx_peaks) == 0:
        
        keys_ls = ['v_peaks', 't_peaks', 'v_threshold', 't_threshold', 'idx_threshold',
                   'v_amplitude', 't_toPeak', 'v_AHP', 't_AHP', 'idx_AHP',
                   'v_AHP_amplitude', 't_to_AHP', 't_rise', 'FWHM', 'v_HM', 't1_HM',
                   't2_HM', 'idx_step']
        
        dataframe_dict = dict.fromkeys(keys_ls, [np.nan])
    
    APs_dataframe = pd.DataFrame(dataframe_dict)

    return APs_dataframe


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

