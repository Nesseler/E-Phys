# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:15:07 2024

@author: nesseler
"""


def clean_OnPath_column_to_path_ID_n_label(coordinates_dataframe):
    
    onPath_column = coordinates_dataframe['On Path']

    for i, txt in enumerate(onPath_column):
        path_ID = [int(s.replace("-", "").strip("()")) for s in txt.split() if s.replace("-", "").strip("()").isdigit()]
        
        if len(path_ID) > 1:
            raise ValueError('More than one integer in path_ID')
        elif len(path_ID) < 1:
            raise ValueError('Less than one integer in path_ID')
        elif len(path_ID) == 1:
            path_ID = path_ID[-1]
        
        
        # txt_ls = txt.split()
            
        if 'axon' in txt:
            path_label = 'axon'
        elif 'soma' in txt:
            path_label = 'soma'
        else:
            path_label = 'dendrite'
            
        if path_ID == 1:
            path_label = 'soma'
    
        coordinates_dataframe.at[i, 'path_ID'] = path_ID
        coordinates_dataframe.at[i, 'path_label'] = path_label
        
    coordinates_dataframe.drop(columns = ['On Path'], inplace = True)

    return coordinates_dataframe


import numpy as np

def calc_polar_histo_binangles(n_bins = 8):
    '''
    Function calculates border angles of bins in rad for the polar plot histogram.
    Parameters: 
        n_bins (int): Number of bins in polar histogram. Default is 8.
    Returns:
        bin_angles (list): List of border angles in rad.
    '''
    
    # step size is set bin number of resulting bins and 2*pi
    bin_stepsize = (2 * np.pi) / n_bins
    
    # start point: half of step size
    # because of rotated polar bins
    bin_start = bin_stepsize / 2
    
    # calc bin borders
    bin_angles = np.arange(bin_start, np.pi * 2, bin_stepsize)

    # roll bin borders
    bin_angles = np.roll(bin_angles, 1)    

    return bin_angles
    
    
    
    