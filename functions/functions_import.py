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
    '''
    '''
    
    # excel sheet with PGF indices as lookup table
    lookup = pd.read_excel(table_file,
                           sheet_name = sheet_name,
                           index_col = 'cell_ID')
    
    # limit lookup table to specific PGF
    lookup_idx = lookup[PGF].dropna().index
    lookup = lookup.loc[lookup_idx, :]
    
    # get indices of current cell with the dataframe containing all indices    
    group_idx = int(lookup.at[cell_ID, 'group'])-1
    
    # get series index
    series_idx = lookup.at[cell_ID, f'{PGF}']
    
    # check for multiple protocols
    if type(series_idx) is int:
    
        # define single series index
        series_idx = series_idx - 1
        
    elif type(series_idx) is np.float64:
        
            # define single series index
            series_idx = int(series_idx) - 1
    
    elif type(series_idx) is str:
        
        # find series indices
        series_idx = [int(i) - 1 for i in series_idx.split(',') if i.isdigit()]
    
    else:
        print(type(series_idx))
        raise ValueError(f'{cell_ID} Series index not found. Neither int nor str!')

    # construct traceIndex with indices
    traceIndex = [group_idx, series_idx, 0, 0]

    # call on data file with indices from dataframe above
    current_file = lookup.at[cell_ID, 'file']

    data_file_path = os.path.join(raw_data_dir, current_file + '.dat')

    data_file_path_str = fr"{data_file_path}"
    
    return traceIndex, data_file_path_str


def get_onlyfiles_list(dir_path):

    from os import listdir
    from os.path import isfile, join  

    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]




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



def get_vc_leak_data(file_path, traceIndex, scale):
    '''
    Function gets path to HEKA file and traceIndex and returns pandas dataframes,
    for current, voltage, leak current and time, where each column represents a 
    single step. The function also returns the sampling rate.
    Parameters:
        file_path : full path for your Heka .dat file
        traceIndex : Index in HEKA tree structure [2,6,0,0] 
                      [Group, Series, Sweep, Trace]
        scale : {'ms', 's'} Time scale in which the time series is being 
                calculated. 'ms' Milliseconds is default.    
    Returns:
        i : current
        v : voltage
        ileak : leak current
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
    ileak = data[:,2,:] * 1e12 
    
    # transpose arrays
    v = np.transpose(v)
    i = np.transpose(i)
    ileak = np.transpose(ileak)
    
    t = calc_time_series(v, SR, scale=scale)
    
    return i, v, ileak, t, SR, n_steps[0]




def get_PSCs_steps(PGF = 'vc_Erest_3min_ctrl', cell_ID = 'E-304', sheet_name = 'PGFs_Syn'):
    '''
    This functions retrieves the steps of PSC recordings that are to be
    included into the analysis.
    Parameters:
        PGF : str, name of PGF and column in lookup table
        cell_ID : str, cell identifier (like: E-000)
        sheet_name : str, sheet name in lookup table
    Returns:
        step_idc: list of int, step indices to be included
    '''
    
    # add str to PGF
    PGF = PGF + '_steps'
    
    # excel sheet with PGF indices as lookup table
    lookup_table = pd.read_excel(table_file,
                                 sheet_name = sheet_name,
                                 index_col = 'cell_ID')
    
    # get listed steps
    steps_str = lookup_table.at[cell_ID, PGF]
    
    # step for multiple entries per cell
    if type(steps_str) == pd.Series:
        steps_str = steps_str[steps_str.notna()].values[0] # .to_string()  
 
    # convert to list of ints
    steps = [int(i) - 1 for i in steps_str.split(',') if i.isdigit()]
    
    return steps
    


def get_MetaData(cell_IDs = None):
    '''
    This function loads the Metadata sheet and limits it to only the cell_IDs
    specified in the cell_IDs list.
    Parameters:
        cell_IDs: list of str, list of cell_ID strings like 'E-313'
    Returns:
        MetaData: pandas Dataframe, containing the metadata
    '''
    
    from parameters.directories_win import table_file
    
    MetaData = pd.read_excel(table_file,
                             sheet_name="MetaData",
                             index_col='cell_ID')

    if cell_IDs:
        MetaData = MetaData.loc[cell_IDs, :]
    
    return MetaData
