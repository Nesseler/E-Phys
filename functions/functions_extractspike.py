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



def get_threshold_crossing_closest_to_peak_and_below_value(data_v, data_dvdt, idx_peak, threshold, change_direction):
    '''

    
    '''
    # datay = spike_dvdt_pre_peak
    # threshold = 5
    # change_direction = -1 # -1: below to above, 1: above to below
    

    # crossing of threshold 
    data_below_th = np.where(data_dvdt <= threshold, 1, 0)
    data_th_crossings = np.diff(data_below_th)
    
    # np.where provides tuple with np array
    data_th_crossings = np.where(data_th_crossings == change_direction)[0]
    
    # calc distance to peak index
    dist_to_peak = np.abs(np.subtract(data_th_crossings, idx_peak))
    
    # get index with minimal distance to peak 
    idx_min = np.argmin(dist_to_peak)
    idx_th = data_th_crossings[idx_min]
    
    above_repol_end_threshold = True
    while above_repol_end_threshold:
        if data_v[idx_th] > -20:
            idx_th = data_th_crossings[idx_min-1]
        else:
            above_repol_end_threshold = False

    # return index
    return idx_th  



def extract_spike(t, v, dvdt, idx_peak):
    '''
    Uses set dvdt thresholds to define spike and returns the indices, that 
    contain the spike.
    
    '''
    dvdt_p_threshold = dvdt_threshold
    dvdt_n_threshold = -3


    if False:
        plt.hlines([dvdt_n_threshold, dvdt_p_threshold], -100, 60, colors = 'gray', linestyle = '--', alpha = 0.5)
        plt.plot(v, dvdt, c = 'gray', linestyle = '-')


    ### rising phase - threshold ###
    
    # limit dvdt to before peak
    spike_dvdt_pre_peak = dvdt[:idx_peak]

    
    try:
        # get index of the threshold
        idx_th = get_threshold_crossing_closest_to_peak(spike_dvdt_pre_peak, len(spike_dvdt_pre_peak)-1, dvdt_p_threshold, -1)
        
    except ValueError:
        raise ValueError('AP threshold not crossed')


    ### repolarisation phase - AHP ###
    spike_dvdt_post_peak = dvdt[idx_peak:]
    spike_v_post_peak = v[idx_peak:]

    # error handling: 
        # problem: negative threshold not crossed
        # iteratively move threshold
        
    # get index of AP end point
    # which is equivalent to the AHP
    threshold_adapt = True
    iterations = 0
    while threshold_adapt:
        try: 
            idx_AHP = get_threshold_crossing_closest_to_peak_and_below_value(spike_v_post_peak, spike_dvdt_post_peak, 0, dvdt_n_threshold, -1) + idx_peak
            threshold_adapt = False
        except:
            iterations += 1
            threshold_change = -1 * iterations
            dvdt_n_threshold = dvdt_n_threshold + threshold_change
            threshold_adapt = True
            print(threshold_change)
            
    
    # problem 2: ohmic repolarisation of membrane during downstroke
    
    
            
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
        plt.plot(t, v, 'gray')
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

import pandas as pd


def get_AP_parameters(t_spiketrain, v_spiketrain, dvdt_spiketrain, idc_spikes, SR=50e3):
    """
    Function calculates all parameters associated with an action potential (AP).
    Parameters:
        t_spiketrain : timeseries array of the spiketrain
        v_spiketrain : voltage trace of the spiketrain
        dvdt_spiketrain : first derivative of the voltage trace
        idc_spikes : list of spike/peak indices in trace
        SR : Sampling rate in Hz. Default is 50 kHz.
    Returns:
        AP_parameters: Pandas Dataframe of all parameters for the provided peaks.
            idx_peak_in_spiketrain : index of spike in spiketrain
            v_peaks : peak voltage
            t_peaks : time point of peak
            v_threshold : voltage at threshold
            t_threshold : time point of threshold
            idx_threshold : index of threshold
            v_amplitude : spike amplitude
            t_toPeak : time to peak (time between threshold and peak)
            t_rise : time between 20 % and 80 % of spike amplitude
            v_AHP : voltage at spike afterhyperpolarisation
            t_AHP : time point of spike afterhyperpolarisation
            idx_AHP : index of spike afterhyperpolarisation
            v_AHP_amplitude : amplitude of spike afterhyperpolarisation (v_threshold - v_AHP)
            t_toAHP : time to spike afterhyperpolarisation (t_AHP - t_peak)
            FWHM : spike full width at half maximum
            v_HM : spike half maximal amplitude
            t1_HM : first time point of spike at half maximum
            t2_HM : second time point of spike at half maximum
        spike_v : voltage trace of isolated spike(s)
    """     

    keys_ls = ['idx_peak_in_spiketrain',
               'v_peaks', 't_peaks', 
               'v_threshold', 't_threshold', 'idc_threshold',
               'v_amplitude', 
               't_toPeak', 't_rise', 
               'FWHM', 'v_HM', 't1_HM', 't2_HM',
               'v_AHP', 't_AHP', 'idx_AHP', 'v_AHP_amplitude', 't_toAHP']
    
    # dataframe_dict = dict.fromkeys(keys_ls, [np.nan])
    
    APs_dataframe = pd.DataFrame(columns = keys_ls)
    
    
    if len(idc_spikes) > 0:
        
        # set spikes indices
        APs_dataframe['idx_peak_in_spiketrain'] = np.arange(len(idc_spikes))
        
        spike_vs = [None] * len(idc_spikes)
     
        # loop through spikes in list of peaks
        for i, idx_peak in enumerate(idc_spikes):
            
            ### peak
            
            # get voltage at peak
            v_peak = v_spiketrain[idx_peak]
            APs_dataframe.at[i, 'v_peaks'] = v_peak
            
            # get time spike
            t_peak = t_spiketrain[idx_peak]
            APs_dataframe.at[i, 't_peaks'] = t_peak
                      
        
            # extract entire spike shapes
            spike_idc, spike_t, spike_v, spike_dvdt = extract_spike(t_spiketrain, v_spiketrain, dvdt_spiketrain, idx_peak)
            
            # convert spike_v to pandas Series
            spike_v = pd.Series(spike_v)
        
            ### threshold
            
            # get index of threshold
            idx_threshold = spike_idc[0]
            APs_dataframe.at[i, 'idc_threshold'] = idx_threshold
        
            # get voltage at threshold
            v_threshold = v_spiketrain[idx_threshold]
            APs_dataframe.at[i, 'v_threshold'] = v_threshold
            
            # get time at threshold
            t_threshold = t_spiketrain[idx_threshold]
            APs_dataframe.at[i, 't_threshold'] = t_threshold
            
            ### amplitude
            
            # get AP amplitude
            v_amplitude = np.abs(v_threshold - v_peak)
            APs_dataframe.at[i, 'v_amplitude'] = v_amplitude
            
            ### time to peak
            
            # get time to peak
            t_toPeak = t_peak - t_threshold
            APs_dataframe.at[i, 't_toPeak'] = t_toPeak
            
            ### rise time
            
            # rise time is calculated between 20 % and 80 % of the voltage amplitude.
            v_20perc = v_threshold + (v_amplitude * 0.2)
            v_80perc = v_threshold + (v_amplitude * 0.8) 
            
            # get voltage trace between 20 % and 80 % amplitude
            v_rise = v_spiketrain[idx_threshold:idx_peak]
            v_rise = np.where((v_rise > v_20perc) & (v_rise < v_80perc))[0]
            
            t_rise = len(v_rise) / (SR / 1e3)
            APs_dataframe.at[i, 't_rise'] = t_rise
            
            ### FWHM
            
            # get half maximum 
            HM = v_amplitude / 2
            
            # get voltage at half maximum            
            v_HM = v_threshold + HM
            APs_dataframe.at[i, 'v_HM'] = v_HM
            
            # use thresholding to find data points above half maximum
            spike_v_above_HM = np.where(spike_v >= v_HM, 1, 0)
            spike_v_change = np.diff(spike_v_above_HM)        
            idx_change = np.where(spike_v_change != 0)[0]
            
            # get start and end time points for half max AP
            t1_HM = spike_t[idx_change[0]]
            t2_HM = spike_t[idx_change[1]]
            
            # calc FWHM
            FWHM = t2_HM - t1_HM
            
            # populate dataframe
            APs_dataframe.at[i, 't1_HM'] = t1_HM
            APs_dataframe.at[i, 't2_HM'] = t2_HM
            APs_dataframe.at[i, 'FWHM'] = FWHM

            ### AHP
            
            # get index of AHP
            idx_AHP = spike_idc[-1]
            APs_dataframe.at[i, 'idx_AHP'] = idx_AHP
            
            # get voltage of AHP
            v_AHP = v_spiketrain[idx_AHP]
            APs_dataframe.at[i, 'v_AHP'] = v_AHP
            
            # get time of AHP
            t_AHP = t_spiketrain[idx_AHP]
            APs_dataframe.at[i, 't_AHP'] = t_AHP
            
            # get time to AHP
            t_toAHP = t_AHP - t_peak
            APs_dataframe.at[i, 't_toAHP'] = t_toAHP
            
            # get AHP amplitude
            v_AHP_amplitude = v_AHP - v_threshold
            APs_dataframe.at[i, 'v_AHP_amplitude'] = v_AHP_amplitude
            
            if len(idc_spikes) > 1:
                spike_vs[i] = spike_v
       
    else:
        APs_dataframe = pd.DataFrame(np.nan, index = [0], columns = keys_ls)
        spike_v = pd.Series([], dtype = float)
        
    if len(idc_spikes) > 1:
        spike_v = spike_vs

    return APs_dataframe, spike_v










