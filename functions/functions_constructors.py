# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 18:10:15 2024

@author: nesseler
"""

import numpy as np
       
def construct_current_array(i_hold, n_steps, parameters_dict, SR_ms):
    i_start = parameters_dict['i_start']
    i_delta = parameters_dict['i_delta']
    
    i_start_rel = i_start - i_hold
    
    i_steps = np.arange(i_start_rel, i_start_rel + (i_delta * n_steps), i_delta)
    
    i = [None] * n_steps
    
    for idx, i_stim in enumerate(i_steps):
        i_pre = np.full(int((SR_ms * parameters_dict['t_pre'])), i_hold)
        i_stim = np.full(int((SR_ms * parameters_dict['t_stim'])), i_stim)
        i_post = np.full(int((SR_ms * parameters_dict['t_post'])-1), i_hold)
        
        i_step = np.concatenate((i_pre, i_stim, i_post))
        i[idx] = i_step
    
    return i, i_steps