# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:18:05 2024

@author: nesseler
"""

from numpy import arange

# define sholl radius step size
sholl_step_size = 1 # um
sholl_profile_range = arange(1, 500 +sholl_step_size, sholl_step_size)