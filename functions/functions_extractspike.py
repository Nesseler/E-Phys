# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:15:15 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import numpy as np

from parameters.parameters import dvdt_threshold


def get_threshold_crossing_closest_to_peak(data, idx_peak, threshold, change_direction):
    '''
    Functions returns the index closest to and idx_peak given that crosses
    the given threshold.
    
    '''
    # datay = spike_dvdt_pre_peak
    # threshold = 5
    # change_direction = -1 # -1: below to above, 1: above to below
    

    # crossing of positive threshold 
    data_below_th = np.where(data <= threshold, 1, 0)
    data_th_crossings = np.diff(data_below_th)
    
    # np.where provides tuple with np array
    data_th_crossings = np.where(data_th_crossings == change_direction)[0]
    
    # calc distance to peak index
    dist_to_peak = np.abs(np.subtract(data_th_crossings, idx_peak))
    
    # get index with minimal distance to peak 
    idx_th = data_th_crossings[np.argmin(dist_to_peak)]

    # return index
    return idx_th    



def extract_spike(t, v, dvdt, idx_peak):
    '''
    Uses set dvdt thresholds to define spike and returns the indices, that 
    contain the spike.
    
    '''
    dvdt_p_threshold = dvdt_threshold
    dvdt_n_threshold = -1


    if False:
        plt.hlines([dvdt_n_threshold, dvdt_p_threshold], -100, 60, colors = 'gray', linestyle = '--', alpha = 0.5)
        plt.plot(v, dvdt, c = 'gray', linestyle = '-')


    ### rising phase - threshold ###
    
    # limit dvdt to before peak
    spike_dvdt_pre_peak = dvdt[:idx_peak]

    # get index of the threshold
    idx_th = get_threshold_crossing_closest_to_peak(spike_dvdt_pre_peak, len(spike_dvdt_pre_peak)-1, dvdt_p_threshold, -1)


    ### repolarisation phase - AHP ###
    spike_dvdt_post_peak = dvdt[idx_peak:]

    # error handling: 
        # problem: negative threshold not crossed
        # iteratively move threshold
        
    # get index of AP end point
    # which is equivalent to the AHP
    threshold_adapt = True
    iterations = 0
    while threshold_adapt:
        try: 
            idx_AHP = get_threshold_crossing_closest_to_peak(spike_dvdt_post_peak, 0, dvdt_n_threshold, -1) + idx_peak
            threshold_adapt = False
        except:
            iterations += 1
            threshold_change = -1 * iterations
            dvdt_n_threshold = dvdt_n_threshold + threshold_change
            threshold_adapt = True
            print(threshold_change)
            
    # construct spike index list
    spike_idc = np.arange(idx_th, idx_AHP, 1, dtype = int)
    
    # limit t, v, and dvdt to just the spike
    spike_t = t[spike_idc]
    spike_v = v[spike_idc]
    spike_dvdt = dvdt[spike_idc] 


    if False:
        plt.scatter(v[idx_th], dvdt[idx_th], marker = 'x', c = 'm')
        plt.scatter(v[idx_AHP], dvdt[idx_AHP], marker = 'x', c = 'c')
        plt.ylim([-100, 250])
        plt.xlim([-100, 60])

        plt.show()
    
    if False:
        plt.plot(t, v, 'k')
        plt.plot(spike_t, spike_v, 'm')
        plt.title(iterations)
        plt.ylim([-100, 60])
        
        plt.show()
        
    return spike_idc, spike_t, spike_v, spike_dvdt
    



# import pandas as pd
# from os.path import join
# import scipy as sc
# from parameters.directories_win import quant_data_dir
# from parameters.parameters import min_peak_distance, min_peak_prominence

# mypath = join(quant_data_dir, '1stAP')

# from os import listdir
# from os.path import isfile, join
# onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

        
# for file in onlyfiles:
#     ### spike extraction ###
#     AP_df = pd.read_excel(join(quant_data_dir, '1stAP', file))
    
    
#     SR = 100000
    
#     spike_t = AP_df['t']
#     spike_v = AP_df['v']
#     spike_dvdt = AP_df['dvdt']
    
    
#     idx_peak, dict_peak = sc.signal.find_peaks(spike_v, 
#                                                   prominence = min_peak_prominence, 
#                                                   distance = min_peak_distance * (SR/1e3))
    
#     idx_peak = idx_peak[0]


#     spike_idc = extract_spike(spike_t, spike_v, spike_dvdt, idx_peak)










