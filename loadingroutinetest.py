# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:20:10 2024

@author: nesseler
"""




# loading routine for current clamp data


'''
    parameters: cell_ID, PGF, PGF_parameters, t_scale, filter_method, current_method

'''








## replace all this:
    
    # # get the traceIndex and the file path string for data import functions
    # traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID)
    
    # # get IF data form file
    # i, v, t, SR, n_steps = get_IF_data(file_path, traceIndex, 'ms')
    
    # # sampling rate in ms
    # SR_ms = SR / 1e3
    
    # # concatenate individual steps
    # n_points = int(np.shape(i)[0] * np.shape(i)[1])
    
    # v_concat = v.flatten() 
    
    # t_ms = calc_time_series(v_concat, SR)
    # t_s = calc_time_series(v_concat, SR, scale = 's')
    
    # # filter voltage (to vf)
    # vf = butter_filter(v_concat, 
    #                    order = 3,
    #                    cutoff = 1e3,
    #                    sampling_rate = SR)
    
    # ### construct current dataframe
    # i_hold = I_hold_table.at[cell_ID, PGF]
    
    # # calculate current steps relative to I_hold
    # ## rounded to nearest 5
    # i_hold_rounded = round_to_base(i_hold, 5)
    
    # # get current arrays and list of input current relative to i_hold
    # i, i_input = construct_current_array(i_hold = i_hold_rounded,
    #                                      n_steps = n_steps,
    #                                      parameters_dict = cc_IF_parameters,
    #                                      SR_ms = SR_ms)

    # step_dur = cc_IF_parameters['t_pre'] + cc_IF_parameters['t_stim'] + cc_IF_parameters['t_post']
    # step_points = step_dur * SR_ms

    # # loop through steps to limit voltage trace
    # for step_idx in np.arange(0, n_steps, 1):
        
    #     # calc start and stop indices for step
    #     start_idx = int(step_points * step_idx)
    #     stop_idx = int(start_idx + step_points)
        
    #     # set voltage trace of step
    #     v[step_idx] = vf[start_idx:stop_idx]