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