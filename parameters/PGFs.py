# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 17:14:08 2023

@author: nesseler

PGF descriptions
"""

import pandas as pd
import numpy as np

# %% cc_APs_xHz
cc_APs_parameters = {
                     '1Hz' :  {'t_pre' : 495., 't_stim' : 10., 't_post' : 495.},
                     '5Hz' :  {'t_pre' : 95. , 't_stim' : 10., 't_post' : 95. },
                     '10Hz' : {'t_pre' : 45. , 't_stim' : 10., 't_post' : 45. },
                     '30Hz' : {'t_pre' : 12.5, 't_stim' : 10., 't_post' : 12.5},
                     '50Hz' : {'t_pre' : 5.  , 't_stim' : 10., 't_post' : 5.  },
                     '75Hz' : {'t_pre' : 2.5 , 't_stim' : 10., 't_post' : 2.5 }
                     }

cc_APs_n_stims = 100

# create list with interstimulus interval per stimulation frequency
cc_APs_ISIs = dict(zip(cc_APs_parameters.keys(), [sum(cc_APs_parameters[f].values()) for f in cc_APs_parameters.keys()]))
cc_APs_total_dur = dict(zip(cc_APs_parameters.keys(),[ISI * cc_APs_n_stims for ISI in cc_APs_ISIs.values()]))

# create dataframe with all stimuation time of every frequency
cc_APs_t_stims_df = pd.DataFrame()

for freq in cc_APs_parameters.keys():
    
    t_pre = cc_APs_parameters[freq]['t_pre']
    
    t_stims = [t_pre + cc_APs_ISIs[freq] * i for i in np.arange(cc_APs_n_stims)]

    cc_APs_t_stims_df[freq] = t_stims
    
    

# %% cc_th1AP

cc_th1Ap_parameters = {'t_pre' : 240.,
                       't_stim': 10.,
                       't_post': 250.,
                       'i_delta' : 10,
                       'i_start' : 0}


# %% cc_cnt_rest

cc_cntrest_parameters = {'t' : 600,
                         'i_hold': 0}


# %% cc_IF

cc_IF_parameters = {'t_pre' : 250, #ms
                    't_stim' : 1000, #ms
                    't_post' : 250, #ms
                    'i_delta' : 5, #pA
                    'i_start' : -50, #pA
                    'max_n_steps' : 71}


# %% cc_sag

cc_sag_parameters = {'t_pre' : 250, #ms
                     't_stim' : 1000, #ms
                     't_post' : 250, #ms
                     'i_delta' : 10, #pA
                     'i_start' : -100, #pA
                     'max_n_steps' : 21,
                     'v_hold_pre' : -85} #mV

# %% vc_rest_EPSC_parameters

# vc_rest_EPSC_parameters = {'t_pre' : 100, #ms
#                            't_stim' : 30000, #ms
#                            't_post' : 10, #ms
#                            'max_n_steps' : 3,
#                            'v_hold_pre' : -85} #mV


# %% vc_TTX_washin / vc_TTX+Cd_washin / vc_TTX_washin_leak / vc_TTX+Cd_washin_leak

vc_TTX_washin_parameters = {'t_pre' : 2490, #ms
                            't_stim': 20, #ms
                            't_post': 2490, #ms
                            'v_stim': 0 #mV
                            } 

vc_TTX_washin_leak_parameters = {'t_pre' : 10, #ms
                                 't_stim': 20, #ms
                                 't_post': 10, #ms
                                 'v_stim': 0 #mV
                                 } 


# %% vc_Erest

vc_Erest_parameters = {'SR' : 100000, #Hz
                       't' : np.arange(0,(6*30), 1/100000), #s
                       'n_steps' : 6,
                       'dur_steps' : 30} #s
