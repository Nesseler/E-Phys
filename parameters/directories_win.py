# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 16:13:15 2023

@author: nesseler

list of all directories on windows machine

"""

### ePhys ###

ePhys_parent = 'Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA'
# ePhys_parent = '/Users/moritznesseler/ePhys-BAOT_MeA'

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
cell_descrip_syn_dir = ePhys_parent + '/cell_descriptors_syn'

# directory for figures to be saved
figure_dir = ePhys_parent + '/figures'

# directory for temporary figures
tempfigs_dir = ePhys_parent + '/figures/temp_figs'

# directory for verification plots to be saved
vplot_dir = ePhys_parent + '/vplots'

# directory for clustering
clustering_dir = ePhys_parent + '/clustering'


