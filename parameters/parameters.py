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
dvdt_n_threshold = -3


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



# %% cc_IF
# number of spikes to calculate the initial instant. firing rate
n_APs_initial_inst_freq = 3

# time period for initial inst freq
t_init_inst_freq = 100 # ms

t_expo_fit = 150. #ms

# set a guess for exponential fit
popt_guess = [50, 0.01, -100]

# useful step r_squared threshold
r_squared_thresh = 0.9

# set parameters to find peaks
min_peak_prominence_ccIF = 30 #(mV)
min_peak_distance_ccIF = 1 #ms
min_max_peak_width_ccIF = [0.5, 10]


# cc_IF adapation

# initial guesses for linear fit in ISIs
adaptation_popt_guess_linear_ISIs = [-1, 25]

adaptation_n_lastspikes = 4



