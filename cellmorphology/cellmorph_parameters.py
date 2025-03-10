# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:18:05 2024

@author: nesseler
"""

import numpy as np
import pandas as pd

# %% sholl plots

# define sholl radius step size
sholl_step_size = 1 # um
sholl_profile_range = np.arange(1, 500 +sholl_step_size, sholl_step_size)


    
# %% cell coordinates

field_of_view = 590.76 #µm

max_depth = 300 #µm


# %% cell_IDs that will be excluded from analysis

cell_IDs_toDrop = ['E-108', # very small (probably not filled)
                   'E-126', 
                   'E-158',
                   'E-065', # BAOT/MeA
                   'E-070'] # BAOT/MeA

