# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 10:25:44 2025

@author: nesseler
"""


from functions.initialize_packages import *


# %% load dataframes to describe cells

from parameters.directories_win import cell_descrip_syn_dir

# cc_rest
activity            = pd.read_excel(join(cell_descrip_syn_dir, 'cc_rest-syn-activity.xlsx'), index_col = 'cell_ID')

# cc_IF
IF_dict             = pd.read_excel(join(cell_descrip_syn_dir, 'cc_IF-syn-IF_dict.xlsx'), index_col = 'cell_ID')
IF_rheobase         = pd.read_excel(join(cell_descrip_syn_dir, 'cc_IF-syn-IF_rheobase.xlsx'), index_col = 'cell_ID')
adaptation          = pd.read_excel(join(cell_descrip_syn_dir, 'cc_IF-syn-adaptation.xlsx'), index_col = 'cell_ID')

# cc_sag
passive_properties  = pd.read_excel(join(cell_descrip_syn_dir, 'cc_sag-syn-passive_properties.xlsx'), index_col = 'cell_ID')
sag_properties      = pd.read_excel(join(cell_descrip_syn_dir, 'cc_sag-syn-sag_properties.xlsx'), index_col = 'cell_ID')

# cc_th1AP
rheobase            = pd.read_excel(join(cell_descrip_syn_dir, 'cc_th1AP-syn-rheobase.xlsx'), index_col = 'cell_ID')


# %% get cell_IDs & MetaData

cell_IDs = activity.index.to_list()

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)


# %% create one dataframe

celldescriptors = pd.concat([activity,
                             IF_dict,
                             IF_rheobase.drop(columns = ['rheobase_abs', 'rheobase_rel']),
                             adaptation,
                             passive_properties,
                             sag_properties.loc[:, 'sagdelta'], # maybe n_reboundspikes
                             rheobase],
                            axis = 1)
                                
# %% drop parameters

parameters_toDrop = ['n_spikes', 
                     't_spikes', 
                     'activity', 
                     'i_rheobase',
                     'idx_rheobase',
                     'i_halfmax',
                     'idx_maxfreq', 
                     'idx_halfmax', 
                     'idx_maxinitinstfreq',
                     'freq_adaptation_steadystate',
                     'rheobase_abs']

celldescriptors.drop(columns = parameters_toDrop, inplace = True)


# %% calc additional parameters

celldescriptors['delta_vrest2vthres'] = celldescriptors['rheobasespike_vthreshold'] - celldescriptors['v_rest']


# %% sort columns

sorter = ['r_input',
          'tau_mem',
          'c_mem',
          'v_rest',
          'sagdelta',
          'rheobase_rel',
          'rheobasespike_vthreshold',
          'rheobasespike_vamplitude',
          'rheobasespike_ttoPeak',
          'rheobasespike_FWHM',
          'rheobasespike_AHPamplitude',
          'rheobasespike_ttoAHP',
          'delta_vrest2vthres',
          'i_maxfreq',
          'i_maxinitinstfreq',
          'spike_amplitude_adaptation_ratio',
          'spike_FWHM_adaptation_ratio',
          'freq_adaptation_ratio']

celldescriptors = celldescriptors[sorter]


# %% save dataframe

from parameters.directories_win import clustering_dir
celldescriptors.to_excel(join(clustering_dir, 'ePhys_celldescriptors.xlsx'), index_label = 'cell_ID')
