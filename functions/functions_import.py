# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 18:57:17 2024

@author: nesseler
"""
import pandas as pd
import os
import numpy as np

from parameters.directories_win import table_dir, table_file, raw_data_dir

def get_traceIndex_n_file(PGF = 'ccth1AP', cell_ID = 'E-092', sheet_name = 'PGFs'):
    # excel sheet with PGF indices as lookup table
    lookup_table = pd.read_excel(table_file,
                                 sheet_name = sheet_name,
                                 index_col = 'cell_ID')
    
    # get indices of current cell with the dataframe containing all indices    
    group_idx = int(lookup_table.at[cell_ID, 'group'])-1
    series_idx = int(lookup_table.at[cell_ID, f'{PGF}'])-1

    # construct traceIndex with indices
    traceIndex = [group_idx, series_idx, 0, 0]

    # call on data file with indices from dataframe above
    current_file = lookup_table.at[cell_ID, 'file']

    data_file_path = os.path.join(raw_data_dir, current_file + '.dat')

    data_file_path_str = fr"{data_file_path}"
    
    return traceIndex, data_file_path_str


def get_onlyfiles_list(dir_path):

    from os import listdir
    from os.path import isfile, join  

    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]



# %% data import

from HekaIO.HekaHelpers import HekaBundleInfo

from functions.functions_useful import calc_time_series


def get_cc_data(file_path, traceIndex, scale):
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
    SR = bundleTester.getSeriesSamplingRate(traceIndex)
    SR = int(round(SR))
    
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



def get_vc_data(file_path, traceIndex, scale):
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
    SR = bundleTester.getSeriesSamplingRate(traceIndex)
    SR = int(round(SR))
    
    # Get data from a series
    data = bundleTester.getSeriesData(traceIndex)
    
    # get data shape
    data_shape = np.shape(data)
    n_points, n_channels, n_steps = zip(data_shape)
    
    # assign current and voltage to dataframes
    v = data[:,1,:] * 1e3
    i = data[:,0,:] * 1e12
    
    v = np.transpose(v)
    i = np.transpose(i)
    
    t = calc_time_series(v, SR, scale=scale)
    
    return i, v, t, SR, n_steps[0]