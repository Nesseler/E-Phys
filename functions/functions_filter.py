# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:43:23 2025

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


def cheby_filter(data, order = 1, rp = 4, cutoff = 4e3, sampling_rate = 20e3):
    b, a = sc.signal.cheby1(order, rp, cutoff, fs=sampling_rate)

    data_filtered = sc.signal.lfilter(b, a, data)  

    return data_filtered


def merge_filter_split_steps(data, SR):
    '''
    This function gets data in steps, concatenates these steps, to then apply
    a filter and split the concatenated steps back into their original form. 
    NOTE: The filter artifact in the beginning of the first step will be 
    replaced by np.nan values in the first 100 datapoints.
    '''
    
    if type(data) != np.ndarray:
        raise ValueError('Provided data does not match expected data structure!')
    
    # get number of steps
    n_steps = data.shape[0]
    
    # merge
    # concatenate steps
    data = data.flatten('C')
    
    # filter
    # filter all data with 1kHz cutoff
    data = butter_filter(data, order=3, cutoff=1e3, sampling_rate=SR)
    
    # replace first values with nans to eliminate filter artifact
    data[:100] = np.nan
    
    # split
    # split back into steps
    data = np.array_split(data, n_steps)
    
    return data