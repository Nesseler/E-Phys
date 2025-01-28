# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 18:10:15 2024

@author: nesseler
"""

import numpy as np
import pandas as pd

def construct_current_array(i_hold, n_steps, parameters_dict, SR_ms):
    i_start = parameters_dict['i_start']
    i_delta = parameters_dict['i_delta']
    
    i_steps = np.arange(i_start, i_start + (i_delta * n_steps), i_delta)
    
    i = [None] * n_steps
    
    for idx, i_stim in enumerate(i_steps):
        i_pre = np.full(int((SR_ms * parameters_dict['t_pre'])), i_hold)
        i_stim = np.full(int((SR_ms * parameters_dict['t_stim'])), i_stim)
        i_post = np.full(int((SR_ms * parameters_dict['t_post'])), i_hold)
        
        i_step = np.concatenate((i_pre, i_stim, i_post))
        i[idx] = i_step
    
    # calculate relative i_input to i_hold
    i_start_rel = i_start - i_hold
    i_input_rel = np.arange(i_start_rel, i_start_rel + (i_delta * n_steps), i_delta)
    
    return i, i_input_rel



def construct_I_array(cell_ID, n_steps, PGF = 'cc_IF', sheet_name = 'PGFs_Syn', 
                      parameters = {'t_pre' : 250, 't_stim' : 1000, 't_post' : 250, 'SR' : 50000}):
    '''
    This function constructs an array that rebuilds the current applied during
    one protocol.
    Parameter:
        cell_ID : str, like 'E-304', cell identifier
        n_steps : int, number of steps in protocol
        PGF : str, default is 'cc_IF', name of protocol
        sheet_name : str, default is 'PGFs_Syn', sheet name in ePhys database
        parameters : dict, dictionary that includes t_start, t_stim, t_stop 
                     in ms and sampling rate
    Returns:
        i : numpy ndarray, containing the applied current over each step
        i_steps : numpy ndarray, containing one value of current for each step
    '''
    
    # cell_ID = 'E-304'
    # n_steps = 24
    # PGF = 'cc_IF'
    # sheet_name = 'PGFs_Syn'
    # parameters = {'t_pre' : 250, #ms
    #                         't_stim' : 1000, #ms
    #                         't_post' : 250, #ms
    #                         't' : np.arange(0, 1500, 1 / (50000/1e3)),
    #                         'SR' : 50000}

    # lazy load table file
    from parameters.directories_win import table_file
        
    # excel sheet with PGF indices as lookup table
    lookup_table = pd.read_excel(table_file,
                                  sheet_name = sheet_name,
                                  index_col = 'cell_ID')
    
    # get cell
    lookup_cell = lookup_table.loc[cell_ID]
    
    # check for multiple entries per cell
    if type(lookup_cell) == pd.DataFrame:
        lookup_cell = lookup_cell.iloc[0]
    
    # get i_hold and i_step
    i_hold  = lookup_cell[PGF + '-hold']
    i_delta = lookup_cell[PGF + '-stepdelta']
    
    # calc i_steps
    i_steps = np.arange(i_hold, i_hold + (i_delta * n_steps), i_delta)
    
    # get duration and sampling rate and calculate the number of points
    dur = parameters['t_pre'] + parameters['t_stim'] + parameters['t_post']
    SR = parameters['SR']
    n_points = int(dur * (SR/1e3))
    n_points_stim = int(parameters['t_stim'] * (SR/1e3))
    
    # initialize i with all i_hold steps
    i = np.full((n_steps, n_points), i_hold)
    
    # get list of indices for stim period
    idc_stim = np.arange(start = parameters['t_pre'] * (SR/1e3),
                          stop = (parameters['t_pre'] + parameters['t_stim']) * (SR/1e3),
                          dtype = int)
    
    # iterate through steps and redefine i
    for step in range(n_steps):
        
        # define current trace in step
        i_step = np.full((n_points_stim), i_steps[step])
        
        # write to array
        i[step, idc_stim] = i_step
    
    return i, i_steps
    