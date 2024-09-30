#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 13:11:17 2024

@author: moritznesseler
"""

# %% cell_IDs that will be excluded from analysis



# reasons to drop
# E-061 - early experiments excluded
# E-065 - early experiments excluded
# E-067 - only BAOT cell with high R_input
# E-185
# E-188 

cell_IDs_toDrop = ['E-061', 'E-065', 'E-067', 'E-185', 'E-188']# 'E-126', 'E-158']



# parameters to exclude from hierarchical clustering
parameters_toDrop = ['sag_delta',
                     'n_reboundspikes',
                     'reboundspike_v_threshold',
                     'reboundspike_v_amplitude',
                     'reboundspike_t_toPeak',
                     'reboundspike_t_rise',
                     'reboundspike_FWHM',
                     'freq_adaptation_steadystate',
                     # drop 2
                     'rheobase_abs',
                     'n_restspikes',
                     'n_rheobasespikes',
                     'rheobasespike_trise',
                     'rheobasespike_ttoAHP']