#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:03:14 2023

@author: moritznesseler
"""

import numpy as np
import matplotlib as mtl
import matplotlib.pyplot as plt
import pandas as pd


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




from HekaIO.HekaHelpers import HekaBundleInfo

from functions.functions_useful import calc_time_series


def get_sampling_rate(bundleTester, traceIndex):
    
    #from patchview.HekaIO.HekaHelpers import HekaBundleInfo
    
    SR = bundleTester.getSeriesSamplingRate(traceIndex)
    
    SR = int(round(SR))
        
    return SR


def get_IF_data(file_path, traceIndex, scale):
    '''
    Function gets path to HEKA file and traceIndex and returns pandas dataframes,
    for current, voltage, and time, where each column represents a single step.
    The function also returns the sampling rate.
    Parameters:
        file_path : full path for your Heka .dat file
        traceIndex : Index in HEKA tree structure [2,6,0,0] 
                      [Group, Series, Sweep, Trace]
        scale : {'ms', 's'} Time scale in which the time series is being 
                calculated. 'ms' Milliseconds is default.    
    Returns:
        i : current
        v : voltage
        t : time series in specified scale (s, ms)
        SR : sampling rate
        n_steps 
    '''
    
    # read Heka .dat file into a object
    bundleTester = HekaBundleInfo(file_path)
    
    # get sampling rate
    SR = get_sampling_rate(bundleTester, traceIndex)
    
    # Get data from a series
    data = bundleTester.getSeriesData(traceIndex)
    
    # get data shape
    data_shape = np.shape(data)
    n_points, n_channels, n_steps = zip(data_shape)
    
    # assign current and voltage to dataframes
    v = data[:,0,:] * 1e3
    i = data[:,1,:] * 1e12
    
    v = np.transpose(v)
    i = np.transpose(i)
    
    t = calc_time_series(v, SR, scale=scale)
    
    return i, v, t, SR, n_steps[0]







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
    
    from useful_functions import calc_time_series, calc_dvdt
    
    t = calc_time_series(v, SR, 'ms')
    
    dvdt = calc_dvdt(v, t)
    
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