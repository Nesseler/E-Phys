# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 17:36:31 2024

@author: nesseler
"""

import matplotlib.pyplot as plt 
import seaborn as sbn
import pandas as pd
from os.path import join

from parameters.directories_win import cell_morph_traces_sholl_dir, table_file

from getter.get_onlyfiles_list import get_onlyfiles_list


# get onlyfiles list
onlyfiles = get_onlyfiles_list(cell_morph_traces_sholl_dir)

# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# extract only filenames with all coordinates
onlyfiles_sholl_profile = [f for f in onlyfiles if '_profile' in f]

# read sholl profile from csv
sholl_csv = pd.read_csv(join(cell_morph_traces_sholl_dir, onlyfiles_sholl_profile[0]))



 