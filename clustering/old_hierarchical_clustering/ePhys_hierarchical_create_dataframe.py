#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 13:05:40 2024

@author: moritznesseler
"""


import pandas as pd
import scipy as sc
import numpy as np
from os.path import join

from parameters.directories_win import table_file, cell_descrip_dir
from parameters.PGFs import cc_IF_parameters


# load metadata
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')


# %% load dataframes to describe cells

# cc_rest
activity_df = pd.read_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col = 'cell_ID')

# cc_IF
active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')
passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col = 'cell_ID')
fstAP_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-fst_AP_parameters.xlsx'), index_col = 'cell_ID')
IF_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF.xlsx'), index_col = 'i_input')
IF_inst_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst.xlsx'), index_col = 'i_input')
IF_inst_initial_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-IF_inst_initial.xlsx'), index_col = 'i_input')

# cc_sag
sag_df = pd.read_excel(join(cell_descrip_dir, 'cc_sag-sagdelta.xlsx'), index_col = 'cell_ID')

# cc_th1AP
th1AP_iinput = pd.read_excel(join(cell_descrip_dir, 'ccth1AP-fst_AP_i.xlsx'), index_col = 'cell_ID')
th1AP_parameters = pd.read_excel(join(cell_descrip_dir, 'ccth1AP-fst_AP_parameters.xlsx'), index_col = 'cell_ID')

# cc_APs
result_freqs_df = pd.read_excel(join(cell_descrip_dir, 'ccAPs-resul_freq.xlsx'), index_col = 'frequencies')


# cc_IF_adapatation
adaptation_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-spike_adaptation.xlsx'), index_col = 'cell_ID')




# %% create unified dataframe

# create dataframe
celldescriptors_df = activity_df.loc[:, ['v_rest', 'n_spikes']]

# add cc_IF passive and active parameters
celldescriptors_df = pd.concat([celldescriptors_df, passiv_properties_df, active_properties_df], axis = 1)

# calc delta vrest to v_thres of ccIF rheobase spike
# celldescriptors_df['delta_vrest_to_vthres_rheobasespike'] = celldescriptors_df['v_thres_rheobase_spike'] - celldescriptors_df['v_rest']

# add cc_IF rheobase spike parameter
celldescriptors_df = pd.concat([celldescriptors_df, fstAP_df[['v_amplitude', 'FWHM', 't_toPeak', 't_rise', 'v_AHP_amplitude', 't_to_AHP']]], axis = 1)

# add cc_sag
### celldescriptors_df = pd.concat([celldescriptors_df, sag_df[['sag_delta', 'n_reboundspikes']]], axis = 1)
celldescriptors_df = pd.concat([celldescriptors_df, sag_df[['sag_delta', 'n_reboundspikes', 'reboundspike_v_threshold', 'reboundspike_v_amplitude', 'reboundspike_t_toPeak', 'reboundspike_t_rise', 'reboundspike_FWHM']]], axis = 1)

# add cc_IF_adaptation
celldescriptors_df = pd.concat([celldescriptors_df, adaptation_df], axis = 1)


# recalculation
# celldescriptors_df['rheobasespike_ttospike'] = fstAP_df['t_peaks'] - cc_IF_parameters['t_pre']
celldescriptors_df['rheobasespike_ttospike'] = fstAP_df['t_threshold'] - cc_IF_parameters['t_pre']
celldescriptors_df['delta_vrest_to_vthres'] = fstAP_df['v_threshold'] - activity_df['v_rest']
celldescriptors_df['t_to_AHP'] = fstAP_df['t_AHP'] - fstAP_df['t_peaks']

celldescriptors_df.rename(columns = {'n_spikes' : 'n_restspikes',
                                     'v_thres_rheobase_spike' : 'rheobasespike_vthreshold',
                                     'v_amplitude' : 'rheobasespike_vamplitude',
                                     'FWHM' : 'rheobasespike_FWHM',
                                     't_toPeak' : 'rheobasespike_ttopeak',
                                     't_rise' : 'rheobasespike_trise',
                                     'v_AHP_amplitude' : 'rheobasespike_AHPvamplitude',
                                     't_to_AHP' : 'rheobasespike_ttoAHP',
                                     'v_thres_rheobase_spike' : 'rheobasespike_vthreshold'},
                          inplace = True)



# %% drop cell_ID

from clustering.ePhys_hierarchical_parameters import cell_IDs_toDrop

# remove cells that do not contain all analysed values
celldescriptors_df.drop(index = cell_IDs_toDrop, inplace = True)


# %% save dataframe

from parameters.directories_win import hierarchical_dir

celldescriptors_df.to_excel(join(hierarchical_dir, 'ePhys_celldescriptors.xlsx'), index_label = 'cell_ID')




