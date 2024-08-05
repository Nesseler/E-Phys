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



# %% polar plots

orientation_labels = ['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv']

    
# %% cell coordinates

cell_coordinates_field_of_view = 590.76 #µm