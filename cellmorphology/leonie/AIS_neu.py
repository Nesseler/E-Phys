# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 15:07:19 2025

@author: buesching
"""

###axon daten einlesen###
import pandas as pd
import seaborn as sbn
# import matplotlib.pyplot as plt
from os.path import join, dirname
import numpy as np

def get_onlyfiles_list(dir_path):

    from os import listdir
    from os.path import isfile, join  

    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]

def get_figure_size(width = 328.67, height = 165.5):
    mm = 1/25.4
    figsize=(width*mm, height*mm)
    return figsize
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_analysis_dir, cellmorph_metrics_dir, cellmorph_figures_dir

all_color_1='#94475EFF'
all_color_2='#E5A11FFF'
all_color_3='#364C54FF'

axon_table_file = '//fileserver2/AG Spehr BigData/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA' + '/' + 'Morphology1.xlsx'

 
# ID_table_dir='//fileserver2/AG Spehr BigData/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA/cellmorph_data/Grayscale_mask'
# onlyfiles = get_onlyfiles_list(ID_table_dir)
# cell_IDs = ['E' + f_str[1:5] for f_str in onlyfiles]

table_file = '//fileserver2/AG Spehr BigData/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/' + 'ePhys-database.xlsx'
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# get all cell_IDs to be included
cell_IDs = MetaData[MetaData['NWG_poster_leonie'] == 1].index.to_list()

# filter MetaData to include only specified cell_IDs
MetaData = MetaData.loc[cell_IDs, :]

# %%

BAOT_morph_parameters= MetaData[MetaData['Region']=='BAOT']
MeA_morph_parameters= MetaData[MetaData['Region']=='MeA']




#%%
axon_Data = pd.read_excel(axon_table_file,
                         sheet_name="Axon-Identification",index_col='cell_ID')

region_df=MetaData[['Region']]

axon_Data = axon_Data.loc[cell_IDs, :]

axon_Data = pd.concat([axon_Data,region_df],axis=1)


#filtern von source spalte nach dendritic und somatic
source_filter= ['somatic', 'dendritic']
axon_data_filtered= axon_Data[axon_Data['source'].isin(source_filter)]

#neue spalte hinzufügen= length+dinstance to soma

axon_data_filtered['length + distance']= axon_data_filtered['length(µm)']+axon_data_filtered['distance to soma']

#alle drei parameter in eine spalte packen damit nurn noch zwischen dendritic und somatic unterschieden wird 
axon_data_filtered_melted = pd.melt(axon_data_filtered, id_vars=['source'], value_vars=['length(µm)','distance to soma', 'length + distance'],var_name='variable', value_name='value')

#%%

#filtern von Region spalte nach MeA und BAOT
region_filter= ['MeA', 'BAOT']
axon_data_filtered_region= axon_data_filtered[axon_data_filtered['Region'].isin(region_filter)]

#nur dendritic data
axon_data_filtered_region_dendritic = axon_data_filtered_region[axon_data_filtered_region['source']=='dendritic']
axon_data_filtered_region_dendritic_melt = pd.melt(axon_data_filtered_region_dendritic, id_vars=['Region'], value_vars=['length(µm)','distance to soma','length + distance'],var_name='variable', value_name='value')

axon_data_filtered_region_somatic= axon_data_filtered_region[axon_data_filtered_region['source']=='somatic']

# %%

# import half violin function
from functions.functions_plotting import plot_half_violin

# initialize plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *


fig, axs = plt.subplots(nrows = 1,
                        ncols = 2,
                        figsize=get_figure_size(width=260.334, height=90),
                        layout='constrained',
                        dpi=600,
                        sharey = True)

ax= [0,1]
for v_idx, variable in enumerate(['length(µm)', 'distance to soma', 'length + distance']):
    
    for r_idx, region, direction in zip([0, 1], ['BAOT','MeA',], [-1, 1]):
        

        
        data = axon_data_filtered_region_dendritic[axon_data_filtered_region_dendritic['Region'] == region].loc[:, variable]
        
        x = r_idx + (v_idx *2)
        
        # calc mean, etc.
        metrics_mean = data.mean()
        metrics_median = data.median()
        metrics_std = data.std()       
        
        # set swarm x
        swarm_x = v_idx*2 + r_idx

        # plot swarmplot
        sbn.swarmplot(x = [swarm_x] * len(data),
                      y = data,
                      color = colors_dict['primecolor'],
                      ax = axs[1],
                      size = 4,
                      zorder = 1)
    
        plot_half_violin(data = data,
                         ax = axs[1],
                         v_position = x,
                         v_direction = direction,
                         v_offset = -direction * 0.4,
                         v_lw = 1.5,
                         v_width = 0.8,
                         v_color = region_colors[region],
                         v_zorder = 2,
                         v_abs_cutoff = [0, np.nan])
        
        # calc violin position
        x = swarm_x - (direction*0.3)
        
        # errorbar
        axs[1].errorbar(x = x,
                        y = metrics_mean,
                        yerr = metrics_std,
                        fmt='_', 
                        markersize = 6,
                        markerfacecolor = 'none',
                        markeredgewidth = 1.5,
                        capsize = 2,
                        color = region_colors[region],
                        linewidth = 1.5,
                        label = '_nolegend_',
                        zorder = 3)
        
        # plot median
        axs[1].scatter(x = x,
                       y = metrics_median,
                       marker='D', 
                       s = 5,
                       color = region_colors[region],
                       linewidth = 1.5,
                       label = '_nolegend_',
                       zorder = 4)
        
        # x
        xdict = {'ax_min' : 0.5,
                 'ax_max' : 4.5,
                 'pad' : 1.3,
                 'step' : 2,
                 'stepminor' : 2,
                 'label' : ''}
        
        apply_axis_settings(axs[1], axis = 'x', **xdict)
        
        # y
        ydict = {'ax_min' : 0,
                 'ax_max' : 80,
                 'pad' : 0.8,
                 'step' : 20,
                 'stepminor' : 5,
                 'label' : ''}
        

        
        apply_axis_settings(axs[1], axis = 'y', **ydict)

    axs[1].set_xticks(ticks = [0.5, 2.5, 4.5], 
                      labels = ['Length \n(µm)', 'Distance \nto soma \n(µm)', 'Length + \ndistance \n(µm)'],
                      rotation = 90)


#%%
for v_idx, variable in enumerate(['length(µm)', 'distance to soma', 'length + distance']):
    
    for r_idx, region, direction in zip([0, 1], ['BAOT','MeA',], [-1, 1]):
        

        
        data = axon_data_filtered_region_somatic[axon_data_filtered_region_somatic['Region'] == region].loc[:, variable]
        
        x = r_idx + (v_idx *2)
        
        # calc mean, etc.
        metrics_mean = data.mean()
        metrics_median = data.median()
        metrics_std = data.std()       
        
        

        # set swarm x
        swarm_x = v_idx*2 + r_idx

        # plot swarmplot
        sbn.swarmplot(x = [swarm_x] * len(data),
                      y = data,
                      color = colors_dict['primecolor'],
                      ax = axs[0],
                      size = 4,
                      zorder = 1)
    
        plot_half_violin(data = data,
                         ax = axs[0],
                         v_position = x,
                         v_direction = direction,
                         v_offset = -direction * 0.4,
                         v_lw = 1.5,
                         v_width = 0.8,
                         v_color = region_colors[region],
                         v_zorder = 2,
                         v_abs_cutoff = [0, np.nan])
        
        # calc violin position
        x = swarm_x - (direction*0.3)
        
        # errorbar
        axs[0].errorbar(x = x,
                        y = metrics_mean,
                        yerr = metrics_std,
                        fmt='_', 
                        markersize = 6,
                        markerfacecolor = 'none',
                        markeredgewidth = 1.5,
                        capsize = 2,
                        color = region_colors[region],
                        linewidth = 1.5,
                        label = '_nolegend_',
                        zorder = 3)
        
        # plot median
        axs[0].scatter(x = x,
                       y = metrics_median,
                       marker='D', 
                       s = 5,
                       color = region_colors[region],
                       linewidth = 1.5,
                       label = '_nolegend_',
                       zorder = 4)
        # x
        xdict = {'ax_min' : 0.5,
                 'ax_max' : 4.5,
                 'pad' : 1.3,
                 'step' : 2,
                 'stepminor' : 2,
                 'label' : ''}
        
        apply_axis_settings(axs[0], axis = 'x', **xdict)
        
        # y
        ydict = {'ax_min' : 0,
                 'ax_max' : 80,
                 'pad' : None,
                 'step' : 20,
                 'stepminor' : 5,
                 'label' : '',
                 'pad_factor' : 0.05}
        

        
        apply_axis_settings(axs[0], axis = 'y', **ydict)



        

        
        # calc violin position
        x = swarm_x - (direction*0.3)
        

    
    axs[0].set_xticks(ticks = [0.5, 2.5, 4.5], 
                      labels = ['Length \n(µm)', 'Distance \nto soma \n(µm)', 'Length + \ndistance \n(µm)'],
                      rotation = 90)




    
# remove top and right spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

#axs[0].set_xticks(ticks = np.arange(0, 6, 1))
    # remove top and right spines
   # [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]

plt.show()

# save figure
save_figures(fig, 
              figure_name = 'AIS', 
              save_dir = join(cellmorph_figures_dir, 'AIS'),
              darkmode_bool=darkmode_bool, 
              figure_format='both')



