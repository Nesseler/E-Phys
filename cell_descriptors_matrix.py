# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:45:58 2024

@author: nesseler
"""

import pandas as pd
import seaborn as sbn
from os.path import join
import matplotlib.pyplot as plt


from parameters.directories_win import cell_descrip_dir, table_file, figure_dir, cell_morph_descrip_dir

from parameters.PGFs import cc_IF_parameters

from functions.functions_plotting import get_colors, get_figure_size, set_font_sizes, save_figures


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

# sholl metrics
sholl_metrics_df = pd.read_excel(join(cell_morph_descrip_dir, 'sholl_metrics.xlsx'), index_col = 'cell_ID')
sholl_metrics_df = sholl_metrics_df.loc[[cell_ID for cell_ID in activity_df.index.tolist() if cell_ID in sholl_metrics_df.index.to_list()], :]

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
celldescriptors_df['delta_vrest_to_vthres_rheobasespike'] = celldescriptors_df['v_thres_rheobase_spike'] - celldescriptors_df['v_rest']

# add cc_IF rheobase spike parameter
celldescriptors_df = pd.concat([celldescriptors_df, fstAP_df[['v_amplitude', 'FWHM', 't_toPeak', 't_rise', 'v_AHP_amplitude', 't_to_AHP']]], axis = 1)

# add cc_sag
### celldescriptors_df = pd.concat([celldescriptors_df, sag_df[['sag_delta', 'n_reboundspikes']]], axis = 1)
# celldescriptors_df = pd.concat([celldescriptors_df, sag_df[['sag_delta', 'n_reboundspikes', 'reboundspike_v_threshold', 'reboundspike_v_amplitude', 'reboundspike_t_toPeak', 'reboundspike_t_rise', 'reboundspike_FWHM']]], axis = 1)

# recalculation
# celldescriptors_df['rheobasespike_ttospike'] = fstAP_df['t_peaks'] - cc_IF_parameters['t_pre']
celldescriptors_df['rheobasespike_ttospike'] = fstAP_df['t_threshold'] - cc_IF_parameters['t_pre']
celldescriptors_df['delta_vrest_to_vthres'] = fstAP_df['v_threshold'] - activity_df['v_rest']
celldescriptors_df['t_to_AHP'] = fstAP_df['t_threshold'] - fstAP_df['t_AHP']

celldescriptors_df.rename(columns = {'n_spikes' : 'n_restspikes',
                                     'v_thres_rheobase_spike' : 'rheobasespike_vthreshold',
                                     'v_amplitude' : 'rheobasespike_vamplitude',
                                     'FWHM' : 'rheobasespike_FWHM',
                                     't_toPeak' : 'rheobasespike_ttopeak',
                                     't_rise' : 'rheobasespike_trise',
                                     'v_AHP_amplitude' : 'rheobasespike_AHPvamplitude',
                                     't_to_AHP' : 'rheobasespike_ttoAHP',
                                     'v_thres_rheobase_spike' : 'rheobasespike_vthreshold',
                                     'v_thres_rheobase_spike' : 'rheobasespike_vthreshold'},
                          inplace = True)

# add resulting freqs
### for cell_ID in result_freqs_df.columns.to_list():
###     for freq in result_freqs_df.index.to_list():
###         celldescriptors_df.at[cell_ID, f'{freq}_resultfreq'] = result_freqs_df.at[freq, cell_ID]

# cell morphology
###celldescriptors_df = pd.concat([celldescriptors_df, sholl_metrics_df[['max_intersections', 'critical_radius', 'enclosing_radius']]], axis = 1)

# add MetaData
MetaData = MetaData.loc[celldescriptors_df.index.to_list(), :]
# celldescriptors_df = pd.concat([celldescriptors_df, MetaData[['Region', 'Hemisphere', 'Sex', 'Stage', 'reconstructed']]], axis = 1)
celldescriptors_df = pd.concat([celldescriptors_df, MetaData['Region']], axis = 1)

# celldescriptors_df.sort_values(['Region', 'r_input'], inplace = True, ascending=[False, False])

# celldescriptors_df.drop(columns = ['Region'], inplace = True)

# %% normalise all columns

# MinMax-Normalization (xi-x_min) / (x_max - x_min)
norm_celldesc_df = pd.DataFrame()

for col in celldescriptors_df.columns.to_list():


    if 'resultfreq' in col:
                
        freq = int(col.split('Hz_resultfreq')[0])

        norm_celldesc_df[col] = (celldescriptors_df[col] - 0) / (freq - 0)
        
    elif 'Region' in col:
        norm_celldesc_df[col] = celldescriptors_df[col].replace({'MeA': 1, 'BAOT/MeA' : 0.5, 'BAOT' : 0})
        
    else:
        norm_celldesc_df[col] = (celldescriptors_df[col] - celldescriptors_df[col].min()) / (celldescriptors_df[col].max() - celldescriptors_df[col].min())
    
        
# %% heatmap with all

darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)


fig_heat, axs_heat = plt.subplots(nrows = 1,
                                  ncols = 1,
                                  layout = 'constrained',
                                  figsize = get_figure_size())

set_font_sizes()

sbn.heatmap(norm_celldesc_df,
            vmin = 0,
            vmax = 1,
            square = False, 
            ax = axs_heat, 
            cmap="flare_r", 
            yticklabels=False,
            linewidth = 0)    

save_figures(fig_heat, 'heatmap_allregions', figure_dir, darkmode_bool)

###

# norm_celldesc_df.drop(columns = ['Region'], inplace = True)

for idx, region in enumerate(['MeA', 'BAOT']):
    
    region_cellIDs = MetaData[MetaData['Region'] == region].index.to_list()
    
    # region_cellIDs = MetaData.query('Region == "' + region + '" and reconstructed == 1').index.to_list()
    
    print(f'n_cell in {region}: {len(region_cellIDs)}')

    region_norm_celldesc_df = norm_celldesc_df.loc[region_cellIDs, :]
    
    # region_norm_celldesc_df.sort_values('r_input', inplace = True, ascending = False)
    
    colors_dict, region_colors = get_colors(darkmode_bool)

    fig_heat_region, axs_heat_region = plt.subplots(nrows = 1,
                                                    ncols = 1,
                                                    layout = 'constrained',
                                                    figsize = get_figure_size())

    set_font_sizes()

    axs_heat_region.set_title(region)
    
    sbn.heatmap(region_norm_celldesc_df,
                vmin = 0,
                vmax = 1,
                square = False, 
                ax = axs_heat_region,
                cmap="flare_r",
                yticklabels=False,
                linewidth = 0)

    save_figures(fig_heat_region, f'heatmap_{region}', figure_dir, darkmode_bool)

# %%

norm_celldesc_df = norm_celldesc_df[norm_celldesc_df['Region'] != 0.5]



# %% first try of clustering

from sklearn.preprocessing import StandardScaler

data_scaler = StandardScaler()

# celldescriptors_df.drop(columns = 'Region', inplace = True)

# norm_celldesc_df.drop(columns = 'Region', inplace = True)

# scaled_celldescript = data_scaler.fit_transform(celldescriptors_df)


from scipy.cluster.hierarchy import linkage, dendrogram

complete_clustering = linkage(norm_celldesc_df, method="ward", metric="euclidean")

index_cellIDs = region_norm_celldesc_df.index.to_list()

dendro = plt.figure(figsize = get_figure_size())

dendrogram = dendrogram(complete_clustering, labels = norm_celldesc_df.index, leaf_font_size= 8, color_threshold=2.8)

# leave_idc = 

leave_cell_IDs = dendrogram['ivl']

# dendrogram.set_xticks(labels = leave_cell_IDs)

plt.show()


# # 
save_figures(dendro, 'dendrogram', figure_dir, darkmode_bool)



# %% 

norm_celldesc_df = norm_celldesc_df.reindex(leave_cell_IDs)

# norm_celldesc_df = pd.concat([norm_celldesc_df, MetaData['Region']], axis = 1)

# norm_celldesc_df['Region'] = celldescriptors_df['Region'].replace({'MeA': 1, 'BAOT/MeA' : 0.5, 'BAOT' : 0})

fig_heat, axs_heat = plt.subplots(nrows = 1,
                                  ncols = 1,
                                  layout = 'constrained',
                                  figsize = get_figure_size(width = 300))

set_font_sizes()

sbn.heatmap(norm_celldesc_df,
            vmin = 0,
            vmax = 1,
            square = False, 
            ax = axs_heat, 
            cmap="flare_r", 
            yticklabels=False,
            linewidth = 0)    


plt.show()


save_figures(fig_heat, 'heatmap_allregions-post_heri_clustering', figure_dir, darkmode_bool)