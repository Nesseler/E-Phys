# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 19:22:27 2023

@author: nesseler
"""

import scipy as sc
import numpy as np


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


def calc_dvdt(v,t):
    #calculate the first derivate dv/dt
    dv = np.diff(v)
    dt = np.diff(t)
    dvdt = dv / dt
    
    return dvdt

def calc_dvdt_padded(v,t):
    #calculate the first derivate dv/dt
    dv = np.diff(v)
    dt = np.diff(t)
    dvdt = dv / dt
    
    # pad with first value as nan
    dvdt = np.pad(dvdt,
                  pad_width = (1, 0),
                  mode = 'constant',
                  constant_values = (np.nan,))
    
    return dvdt

    

def calc_normed_hist(hist):
    '''
    Function calculate the normed histogram to the max number of points.
    Parameters:
        hist : Array of values that described the number of occurances per bin.
    Returns:
        normed list of values
    '''
    hist_max = np.max(hist)
    return [n / hist_max for n in hist]
    

def single_gaussian(x, amp1, cen1, sigma1):
    '''
    https://en.wikipedia.org/wiki/Normal_distribution
    '''
    return amp1 * (1 / (sigma1 * np.sqrt(2 * np.pi))) * (np.exp( (-1./2.) * ((x - cen1) / (sigma1))**2 ))



def double_gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2):
    '''
    https://en.wikipedia.org/wiki/Normal_distribution
    '''
    return amp1 * (1 / (sigma1 * np.sqrt(2 * np.pi))) * (np.exp( (-1./2.) * ((x - cen1) / (sigma1))**2 )) + amp2 * (1 / (sigma2 * np.sqrt(2 * np.pi))) * (np.exp( (-1./2.) * ((x - cen2) / (sigma2))**2 ))

    
    
def round_to_base(number, base):
    return base * round(number/base)  


def round_up_to_base(number, base):
    return base * np.ceil(number/base)  


def round_down_to_base(number, base):
    return base * np.floor(number/base)

    
# define function of exponential fit
def exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c


def calc_rsquared_from_exp_fit(x_data, y_data, popt):
    # calculate the residuals form fit
    residuals = y_data - exp_func(x_data, *popt)
    
    # calculate residual sum of squares
    ss_res = np.sum(residuals**2)
    
    # calc total sum of squares
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    
    # calc r squared
    r_squared = 1 - (ss_res / ss_tot)
    
    return r_squared


# define function of liner fit
def linear_func(x, a, b):
    return a * x + b


def calc_PDF_fromTimepoints(data: np.array, bins: np.array):

    # calc inter-event interval (IEI)
    IEI = np.diff(data, prepend = 0)    

    # calc mean rate
    mean_rate = 1 / np.mean(IEI)
    
    # calculate the probability density function with bins given
    pdf = mean_rate * np.exp(- mean_rate * bins)

    return pdf


def calc_CDF_fromTimepoints(data: np.array, bins: np.array):
    
    pdf = calc_PDF_fromTimepoints(data, bins)

    cdf = [np.max(pdf) - p for p in pdf]

    return cdf


def minmax_norm(data: np.array, minzero: bool = False):
    
    
    # calc min max normalized version of array 
    if minzero:
        normed = (data - 0) / (np.max(data) - 0)
    else:
        normed = (data - np.min(data)) / (np.max(data)- np.min(data))
    
    return normed


def create_uniform_events(n_events: int = 100, t: float = 180) -> np.array:
    """
    
    Parameters
    ----------
    n_events : int, optional
        DESCRIPTION. The default is 100.
    t : float, optional
        DESCRIPTION. The default is 180.

    Returns
    -------
    uni_events : TYPE
        DESCRIPTION.

    """
    # uniformly distribute the same number of events as measured
    # by randomly draw from uniform distribution between t0 and t_end
    uni_events = sc.stats.uniform.rvs(loc = 0, scale = t, size = n_events)
    
    # sort values to perserve time series aspect
    uni_events = np.sort(uni_events)
        
    return uni_events    
    