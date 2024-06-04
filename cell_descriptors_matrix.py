# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:45:58 2024

@author: nesseler
"""

import pandas as pd
import seaborn as sbn
from os.path import join
import matplotlib.pyplot as plt


from parameters.directories_win import cell_descrip_dir, table_file


# %% load dataframes to describe cells

# cc_rest
activity_df = pd.read_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), index_col = 'cell_ID')

# cc_IF
active_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-active_properties.xlsx'), index_col = 'cell_ID')
active_properties_df.drop(columns = 'rheobase_step_idx', inplace = True)
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

# MetaData
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# %% create unified dataframe

# create dataframe
celldescriptors_df = activity_df.loc[:, ['v_rest', 'n_spikes']]

# add cc_IF passive and active parameters
celldescriptors_df = pd.concat([celldescriptors_df, passiv_properties_df, active_properties_df], axis = 1)

# calc delta vrest to v_thres of ccIF rheobase spike
celldescriptors_df['delta_vrest_vThresRheobaseSpike'] = celldescriptors_df['v_thres_rheobase_spike'] - celldescriptors_df['v_rest']

# add cc_IF rheobase spike parameter
celldescriptors_df = pd.concat([celldescriptors_df, fstAP_df[['v_amplitude', 'FWHM', 't_toPeak', 't_rise']]], axis = 1)

# # add cc_sag
# celldescriptors_df = pd.concat([celldescriptors_df, sag_df[['sag_delta', 'n_reboundspikes']]], axis = 1)
# # celldescriptors_df = pd.concat([celldescriptors_df, sag_df[['sag_delta', 'n_reboundspikes', 'reboundspike_v_threshold', 'reboundspike_v_amplitude', 'reboundspike_t_toPeak', 'reboundspike_t_rise', 'reboundspike_FWHM']]], axis = 1)


# # add resulting freqs
# for cell_ID in result_freqs_df.columns.to_list():
#     for freq in result_freqs_df.index.to_list():
#         celldescriptors_df.at[cell_ID, f'{freq}_resultfreq'] = result_freqs_df.at[freq, cell_ID]

# add MetaData
MetaData = MetaData.loc[celldescriptors_df.index.to_list(), :]
# celldescriptors_df = pd.concat([celldescriptors_df, MetaData[['Region', 'Hemisphere', 'Sex', 'Stage', 'reconstructed']]], axis = 1)
celldescriptors_df = pd.concat([celldescriptors_df, MetaData['Region']], axis = 1)

celldescriptors_df.sort_values(['Region', 'r_input'], inplace = True, ascending=[False, False])

# %% normalise all columns

# MinMax-Normalization (xi-x_min) / (x_max - x_min)
norm_celldesc_df = pd.DataFrame()

for col in celldescriptors_df.columns.to_list():


    if 'resultfreq' in col:
        freq = int(col.split('Hz_resultfreq')[0])
        
        norm_celldesc_df[col] = (celldescriptors_df[col] - 0) / (freq - 0)
        
    if 'Region' in col:
        norm_celldesc_df[col] = celldescriptors_df[col].replace({'MeA': 1, 'BAOT/MeA' : 0.5, 'BAOT' : 0})
        
    else:
        norm_celldesc_df[col] = (celldescriptors_df[col] - celldescriptors_df[col].min()) / (celldescriptors_df[col].max() - celldescriptors_df[col].min())
    
        
# %% heatmap with all

fig_heat, axs_heat = plt.subplots(nrows = 1,
                                  ncols = 1,
                                  layout = 'constrained')

sbn.heatmap(norm_celldesc_df, square = False, ax = axs_heat)    


# %%



for idx, region in enumerate(['MeA', 'BAOT']):
    
    region_cellIDs = MetaData[MetaData['Region'] == region].index.to_list()
    
    # region_cellIDs = MetaData.query('Region == "' + region + '" and reconstructed == 1').index.to_list()
    
    print(f'n_cell in {region}: {len(region_cellIDs)}')

    region_norm_celldesc_df = norm_celldesc_df.loc[region_cellIDs, :]
    
    region_norm_celldesc_df.sort_values('r_input', inplace = True, ascending = False)

    fig_heat_region, axs_heat_region = plt.subplots(nrows = 1,
                                                    ncols = 1,
                                                    layout = 'constrained')

    axs_heat_region.set_title(region)
    
    sbn.heatmap(region_norm_celldesc_df, square = False, ax = axs_heat_region)






















