# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 13:32:13 2024

@author: nesseler
"""

import pandas as pd
from os.path import join
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np

from parameters.directories_win import quant_data_dir
from parameters.parameters import min_peak_distance, min_peak_prominence

from functions.functions_extractspike import extract_spike
from functions.functions_useful import calc_dvdt_padded

mypath = join(quant_data_dir, '1stAP')

from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

        
file = onlyfiles[12]

repol_v_df = pd.DataFrame()
repol_dvdt_df = pd.DataFrame()


for i, file in enumerate(onlyfiles):
    ### spike extraction ###
    AP_df = pd.read_excel(join(quant_data_dir, '1stAP', file))
    
    
    SR = 100000
    
    spike_t = AP_df['t']
    spike_v = AP_df['v']
    spike_dvdt = AP_df['dvdt']
    spike_d2vdt2 = calc_dvdt_padded(spike_dvdt, spike_t)
    
    
    idx_peak, dict_peak = sc.signal.find_peaks(spike_v, 
                                                  prominence = min_peak_prominence, 
                                                  distance = min_peak_distance * (SR/1e3))
    
    idx_peak = idx_peak[0]
    
    
    
    
    ranges = {'v' : [-100, 60], 'dvdt' : [-150, 250], 'd2vdt2' : [-1000, 1000]}
    
    
    def min_max_normalize(data):
    # min max normalization
        # calculated as (value - min) / (max - min)
    
        data_min = np.nanmin(data)
        data_max = np.nanmax(data)
        data_mm_normed = [(value - data_min) / (data_max - data_min) for value in data]
    
        return data_mm_normed
    
    
    
    
    test = min_max_normalize(spike_d2vdt2)
    
    # %%
    
    if True:
        spike_idc, spike_t, spike_v, spike_dvdt = extract_spike(spike_t, spike_v, spike_dvdt, idx_peak)
    
        spike_d2vdt2 = spike_d2vdt2[spike_idc]
    
    # %%
    
    # fig_AP, axs_AP = plt.subplots(nrows = 4, 
    #                               ncols = 1,
    #                               sharex = 'col',
    #                               layout = 'constrained')
    
    # fig_AP.suptitle(file)
    
    # # trace
    # axs_AP[0].plot(spike_t, spike_v)
    
    # # first derivative
    # axs_AP[1].plot(spike_t, spike_dvdt)
    
    # # second derivative
    # axs_AP[2].plot(spike_t, spike_d2vdt2)
    
    
    
    # [axs_AP[i].set_ylim(r) for i, r in enumerate(ranges.values())]
    
    
    # axs_AP[0].set_xlim([spike_t.iat[0], spike_t.iat[-1]])
    
    
    
    
    # plt.show()
    
    
    data = spike_dvdt
    
    mask = np.where(spike_dvdt < 0, 1, 0)
    
    v_b_th = np.where(mask == 1, spike_v, np.nan)
    v_b_th = v_b_th[~np.isnan(v_b_th)]
    
    dvdt_b_th = np.where(mask == 1, spike_dvdt, np.nan)
    d2vdt2_b_th = np.where(mask == 1, spike_d2vdt2, np.nan)
    
    
    dvdt_b_th = dvdt_b_th[~np.isnan(dvdt_b_th)]
    
    dvdt_b_th_mmnormed = pd.Series(min_max_normalize(dvdt_b_th))
    
    
    
    

    
    # # get indices of local max
    # local_max_idc = sc.signal.argrelextrema(dvdt_b_th, np.greater)
    # local_max_dvdt = dvdt_b_th[local_max_idc]
    
    
    # # get inflection point, i.e. local minimum d2vdt2 (should be 0 but approximation here)
    # local_infl_idc = sc.signal.argrelextrema(d2vdt2_b_th, np.less)
    # local_infl_dvdt = dvdt_b_th[local_infl_idc]
    
    

    
    
    # print(sc.integrate.simpson(dvdt_b_th))
    
    # plt.plot(d2vdt2_b_th)
    # plt.scatter(local_infl_idc, d2vdt2_b_th[local_infl_idc], marker = 'x', color = 'k')
    
    # plt.show()
    
    
# %%
    # # get indices of local minima
    local_min_idc = sc.signal.argrelextrema(dvdt_b_th, np.less)
    
    
    
    def exclude_idc_with_values_below_value(idc, values, threshold):
        # set indeces below threshold as nan
        idc = np.where(values[idc] < threshold, idc, np.nan)
        # exclude nan
        idc = idc[~np.isnan(idc)]
        #convert back to ints
        idc = [int(idx) for idx in idc]
        # return idc
        return idc
    
    # limit local minima to below - 10 mV/ms
    local_min_idc = exclude_idc_with_values_below_value(local_min_idc, dvdt_b_th, -10)
    


    # get dvdt values for local min
    local_min_dvdt = dvdt_b_th[local_min_idc]


    # get indices of local max
    local_max_idc = sc.signal.argrelextrema(dvdt_b_th, np.greater)
    
    # limit local maxima to below - 10 mV/ms
    local_max_idc = exclude_idc_with_values_below_value(local_max_idc, dvdt_b_th, -10)
    
    # get dvdt values for local max
    local_max_dvdt = dvdt_b_th[local_max_idc]

    # plt.plot(dvdt_b_th)
    # plt.scatter(local_min_idc, local_min_dvdt, marker = 'x', color = 'r')
    # plt.scatter(local_max_idc, local_max_dvdt, marker = 'x', color = 'm')
    # # # plt.scatter(local_infl_idc, local_infl_dvdt, marker = 'x', color = 'k')
    
    # plt.title(file)
    # plt.ylim([-150, 0])
    # plt.xlim([0, 600])
    
    # plt.show()
    
    
    if len(local_min_idc) == 2:
        distance_min_max = local_max_dvdt - local_min_dvdt[1]
    
        # if distance_min_max > 10:
        #     repol_dvdt_df = pd.concat([repol_dvdt_df, pd.DataFrame({i : dvdt_b_th})], axis = 1)
    
    if len(local_min_idc) == 1:
        
        d2vdt2_b_th = np.array(min_max_normalize(d2vdt2_b_th))
        # d2vdt2_b_th = [value * 100 for value in d2vdt2_b_th]
        local_max_idc, _ = sc.signal.find_peaks(d2vdt2_b_th, height = 0.9)
        # local_max_idc = [int(idx) for idx in local_max_idc]
        local_max_d2vdt2 = d2vdt2_b_th[local_max_idc]
        
        plt.title(file)
        plt.plot(d2vdt2_b_th)
        plt.scatter(local_max_idc, local_max_d2vdt2, marker = 'x', color = 'm')
        plt.show()
    
        # if len(local_max_idc) > 1:
        #     repol_dvdt_df = pd.concat([repol_dvdt_df, pd.DataFrame({i : dvdt_b_th})], axis = 1)       



# %%
    
# for i in repol_dvdt_df.columns.to_list():
#     plt.plot(repol_dvdt_df[i],
#               c = 'gray')

# plt.show()
    
    
    
    