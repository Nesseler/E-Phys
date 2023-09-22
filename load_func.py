#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collection of functions that load specific data

Created on Wed Sep 20 15:02:47 2023

@author: moritznesseler
"""

import pandas as pd
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np

# %%

#def load_mat2df(file_path):

test = sc.io.loadmat('test_data/AMY-20230905-E004-cc_IF_06.mat')

# %%

all_header = list(test.keys())

traces_header = [key for key in all_header if 'Trace' in key]

voltage_headers = [key for key in traces_header if key[-1] == '1']


#create two dataframe with voltage, current, time


#for v_head in voltage_headers:
#    print(v_head)

test3 = np.transpose(test['Trace_2_6_20_1'])

test4 = np.concatenate((np.transpose(test['Trace_2_6_20_1'])[1], np.transpose(test['Trace_2_6_21_1'])[1]))

plt.plot(test4)


plt.show()


# %%

for i in range(44):
    test2 = np.transpose(test['Trace_2_6_' + str(i+1) + '_1'])

    plt.plot(test2[1])
    
    plt.pause(0.05)
    
    print(i)
    
# %%  

step20 = np.transpose(test['Trace_2_6_45_1'])[1]

#convert V to mV
step20 = step20 * 10**3

step20std = np.mean(step20)

print(len(step20)/(20 * 10**3))




#250ms pre & post pulse at 20 kHz Sampling Rate
#indices are 0-5000, 5001-25000, 25001-30000

#sampling rate in Hz
SR = 20000

#pre, post, and pulse duration in s
pre_post_dur = 0.250
pulse_dur = 1

pre_idx = np.arange(0, int(pre_post_dur * SR))
pulse_idx = np.arange(pre_idx[-1]+1, int((pre_post_dur + pulse_dur) * SR))
post_idx = np.arange(pulse_idx[-1]+1, int((pre_post_dur + pulse_dur + pre_post_dur) * SR))

ylimits = [-100, +20]


pulse_mean = np.mean(step20[pulse_idx])
pulse_median = np.median(step20[pulse_idx])

pulse_std = np.std(step20[pulse_idx])


threshold_figure, axs = plt.subplots(2,2, sharey = 'row', tight_layout = True)

axs[0,0].plot(step20,
              label = 'cc_IF')
axs[0,0].set_ylim(ylimits)
axs[0,0].set_xlim([0, post_idx[-1]])
axs[0,0].fill_between(x = pulse_idx, 
                      y1 = ylimits[0],
                      y2 = ylimits[1],
                      color = 'grey',
                      alpha = 0.5,
                      label = 'Analysed time frame')
axs[0,0].set_xlabel('Time [ms]')
axs[0,0].set_ylabel('Voltage [mV]')

n_hist, bins_hist = np.histogram(step20[pulse_idx], 
                                 bins = np.arange(-100, +20, 10))

axs[0,1].barh(bins_hist[:-1], 
              n_hist, 
              height = 10, 
              align = 'edge',
              label = 'All points histogram')

hist_max_points = int(pulse_dur * SR)

axs[0,1].set_xlim([0, hist_max_points])
axs[0,1].plot([0, hist_max_points],
              [pulse_mean]*2,
              color = 'red')
axs[0,1].plot([0, hist_max_points],
              [pulse_median]*2,
              color = 'magenta')

axs[0,1].fill_between(x = [0, hist_max_points], 
                      y1 = [pulse_median-pulse_std*2]*2,
                      y2 = [pulse_median+pulse_std*2]*2,
                      color = 'grey',
                      alpha = 0.5,
                      label = '2 std around median')

axs[0,1].fill_between(x = [0, hist_max_points], 
                      y1 = [pulse_median-pulse_std]*2,
                      y2 = [pulse_median+pulse_std]*2,
                      color = 'grey',
                      alpha = 0.7,
                      label = 'Std around median')


#threshold_figure.legend(loc = 'outside')

plt.show()

    
    
    
    