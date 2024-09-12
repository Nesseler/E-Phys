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


# def get_colorcode(x, y, data_fc, norm=None, cmap='seismic', plot_dict={'c':'k'}):
#     points = np.array([x, y]).T.reshape(-1, 1, 2)
#     segments = np.concatenate([points[:-1], points[1:]], axis=1)

#     return_bool = False

#     if norm is None:
#         # Create a continuous norm to map from data points to colors
#         norm = plt.Normalize(data_fc.min(), data_fc.max())
#         return_bool = True
        
        
#     lc = mtl.collections.LineCollection(segments, cmap=cmap, norm=norm)
#     # Set the values used for colormapping
#     lc.set_array(data_fc)
    
#     return (norm, lc) if return_bool else lc




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



