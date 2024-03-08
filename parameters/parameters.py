# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 17:11:46 2023

@author: nesseler

Parameters files
    contains hard coded variables for analysis
"""


# cc_APs
cc_APs_t_post_stim = 10 # in ms


# set parameters to find peaks
min_peak_prominence = 40 #(mV)
min_peak_distance = 1 #ms

dvdt_threshold = 5


AP_parameters = ['v_peaks',
                 't_peaks',
                  'v_threshold',
                  't_threshold',
                  'idx_threshold',
                  'v_amplitude',
                  't_toPeak',
                  'v_AHP',
                  't_AHP',
                  'idx_AHP',
                  'v_AHP_amplitude',
                  't_to_AHP',
                  't_rise',
                  'FWHM',
                  'v_HM',
                  't1_HM',
                  't2_HM']

# %% cnt_rest

min_spikes_tobe_active = 2

min_spike_in_burst = 4

bin_size_ISI_poisson = 25e-3



