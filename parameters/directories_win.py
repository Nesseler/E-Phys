# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 16:13:15 2023

@author: nesseler

list of all directories on windows machine

"""

### ePhys ###

# ePhys_parent = 'Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA'
ePhys_parent = '/Users/moritznesseler/ePhys-BAOT_MeA'

# directory that contains the excel file with indices and meta data
table_dir = ePhys_parent + '/'
table_file = ePhys_parent + '/ePhys-database.xlsx'
# table_dir = 'C:/Users/nesseler/Desktop/'

# directory that contains raw data
raw_data_dir = ePhys_parent + '/RAW_data'

# directory for analysed values  
quant_data_dir = ePhys_parent + '/quantified_data'

# directory for analysed values that describe one cell
cell_descrip_dir = ePhys_parent + '/cell_descriptors'

# directory for figures to be saved
figure_dir = ePhys_parent + '/figures'

# directory for verification plots to be saved
vplot_dir = ePhys_parent + '/vplots'


### cell morphology ###

# directory for cell morphology measures (SNT traces/Measurements)
# cell_morph_parent = 'Z:/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA'
cell_morph_parent = '/Users/moritznesseler/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA'
cell_morph_descrip_dir = cell_morph_parent + '/SNT_traces/'
cell_morph_measures_dir = cell_morph_parent + '/SNT_traces/Measurements'
cell_morph_figures_dir = cell_morph_parent + 'cellmorph_figures'
cell_morph_traces_coordinates_dir = cell_morph_parent + '/SNT_traces/traces_coordinates'
cell_morph_plots_dir = cell_morph_parent + '/Plots'
cell_morph_traces_sholl_dir = cell_morph_parent + '/SNT_traces/traces_sholl_tables'