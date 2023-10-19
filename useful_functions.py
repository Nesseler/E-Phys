# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 19:22:27 2023

@author: nesseler
"""

import scipy as sc
import numpy as np


def butter_filter(data, order = 1, cutoff = 4e3, sampling_rate = 20e3):
    b, a = sc.signal.butter(order, cutoff, fs=sampling_rate)

    data_filtered = sc.signal.lfilter(b, a, data)  

    return data_filtered


def bessel_filter(data, order = 1, cutoff = 4e3, sampling_rate = 20e3):
    b, a = sc.signal.bessel(order, cutoff, fs=sampling_rate)

    data_filtered = sc.signal.lfilter(b, a, data)  

    return data_filtered


def calc_time_series(data, sampling_rate = 20e3, scale = 'ms'):
    '''
    Calculate a time series array from given data with sampling rate.
    Parameters:
        data : One-dimensional array with voltage in mV
        sampling_rate : Sampling rate in Hz. Default is 20 kHz.
        scale : {'ms', 's'} Time scale in which the time series is being 
                calculated. 'ms' Milliseconds is default.
    Returns:
        t : Time series array.
    
    '''
    #caculate time for data in ms
    if scale == 'ms':
        sampling_rate = sampling_rate / 1e3
    elif scale == 's':
        sampling_rate = sampling_rate
    else:
        raise ValueError('Choose between millisecond (ms) or second time scale (s)')
    
    t_total = len(data) / sampling_rate
    t = np.arange(t_total, step=1/sampling_rate)
    
    return t

#saving the figure
def save_figures(figure, figure_name, save_dir, darkmode_bool):
    
    import os.path
    
    if darkmode_bool == True:
        figure_name += " dark"
    else:
        figure_name += " light"
    
    figure.savefig(os.path.join(save_dir, os.path.normpath(figure_name + ".png")), format = 'png')
    figure.savefig(os.path.join(save_dir, os.path.normpath(figure_name + ".svg")), format = 'svg')