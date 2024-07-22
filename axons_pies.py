# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 15:14:47 2024

@author: nesseler
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:41:00 2024

@author: buesching
"""

###axon daten einlesen###
import pandas as pd
import seaborn as sbn
import matplotlib.pyplot as plt
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

axon_table_file = '//Fileserver/AG Spehr BigData/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA' + '/' + 'Morphology1.xlsx'

 
ID_table_dir='//Fileserver/AG Spehr BigData/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA/Grayscale_mask'
onlyfiles = get_onlyfiles_list(ID_table_dir)
cell_IDs = ['E' + f_str[1:5] for f_str in onlyfiles]

table_file = '//Fileserver/AG Spehr BigData/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/' + 'ePhys-database.xlsx'
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')
MetaData = MetaData.loc[cell_IDs, :]

BAOT_morph_parameters= MetaData[MetaData['Region']=='BAOT']
MeA_morph_parameters= MetaData[MetaData['Region']=='MeA']
#%%
axon_Data = pd.read_excel(axon_table_file,
                         sheet_name="Axon-Identification",index_col='cell_ID')

region_df=MetaData[['Region']]
axon_Data =axon_Data.loc[cell_IDs, :]
axon_Data = pd.concat([axon_Data,region_df],axis=1)
#%% pie chart
pie_chart_dendritic=axon_Data[axon_Data['source']=='dendritic']
dendritic=np.array(len(pie_chart_dendritic))
pie_chart_somatic=axon_Data[axon_Data['source']=='somatic']
somatic=np.array(len(pie_chart_somatic))
pie_chart_non=axon_Data[axon_Data['source']=='non']
non=np.array(len(pie_chart_non))

liste=[somatic, dendritic, non]
#liste.append(somatic)
#liste.append(dendritic)
#liste.append(non)

fig, ax = plt.subplots()
ax.pie(liste, labels=('somatic','dendritic','non'), autopct='%1.1f%%')

fig_dir = '//Fileserver/AG Spehr BigData/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA/plots'

fig.savefig(fig_dir + ' piechart_all' + '.svg', format = 'svg')




#%%pie chart with bars 


#%%BAOT pie chart
BAOT_color_1='#8268ee'
BAOT_color_2='#ad759a'
BAOT_color_3='#002198'
axon_data_BAOT= axon_Data[axon_Data['Region']=='BAOT']


pie_chart_BAOT_dendritic=axon_data_BAOT[axon_data_BAOT['source']=='dendritic']
dendritic_BAOT=np.array(len(pie_chart_BAOT_dendritic))
pie_chart_somatic_BAOT=axon_data_BAOT[axon_data_BAOT['source']=='somatic']
somatic_BAOT=np.array(len(pie_chart_somatic_BAOT))
pie_chart_non_BAOT=axon_data_BAOT[axon_data_BAOT['source']=='non']
non_BAOT=np.array(len(pie_chart_non_BAOT))

liste_BAOT=[somatic_BAOT, dendritic_BAOT, non_BAOT]


fig, ax = plt.subplots()
ax.pie(liste_BAOT, labels=('somatic','dendritic','non'), autopct='%1.1f%%', colors=[BAOT_color_1,BAOT_color_3,BAOT_color_2])


#%%pie chart MeA

axon_data_MeA= axon_Data[axon_Data['Region']=='MeA']
MeA_color_1='#ff8d00'
MeA_color_2='#ba5003'
MeA_color_3='#7d1905'


pie_chart_MeA_dendritic=axon_data_MeA[axon_data_MeA['source']=='dendritic']
dendritic_MeA=np.array(len(pie_chart_MeA_dendritic))
pie_chart_somatic_MeA=axon_data_MeA[axon_data_MeA['source']=='somatic']
somatic_MeA=np.array(len(pie_chart_somatic_MeA))
pie_chart_non_MeA=axon_data_MeA[axon_data_MeA['source']=='non']
non_MeA=np.array(len(pie_chart_non_MeA))

liste_MeA=[somatic_MeA, dendritic_MeA, non_MeA]


fig, ax = plt.subplots()
ax.pie(liste_MeA, labels=('somatic','dendritic','non'), autopct='%1.1f%%', colors=[MeA_color_1,MeA_color_3,MeA_color_2])


# %% combined pie charts

from functions.functions_plotting import get_figure_size, get_colors, set_font_sizes, save_figures

# get colors
darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)


fig_pies, axs_pies = plt.subplots(nrows = 1,
                                  ncols = 3,
                                  figsize = get_figure_size(width = 159.835, height = 77.1),
                                  layout = 'tight')

# all cells
axs_pies[0].pie(liste, 
                labels=('somatic','dendritic','non'), 
                autopct='%1.1f%%',
                textprops={'fontsize': 9})

# MeA
axs_pies[1].pie(liste_MeA, 
                labels=('somatic','dendritic','non'), 
                autopct='%1.1f%%', 
                colors=[MeA_color_1,MeA_color_3,MeA_color_2],
                textprops={'fontsize': 9})

# BAOT
axs_pies[2].pie(liste_BAOT, 
                labels=('somatic','dendritic','non'), 
                autopct='%1.1f%%', 
                colors=[BAOT_color_1,BAOT_color_3,BAOT_color_2],
                textprops={'fontsize': 9})


set_font_sizes(12)

temp_fig_dir = 'C:/Users/nesseler/Desktop/TAC-presentation_data/ePhys'

save_figures(fig_pies, 'AIS_pies', temp_fig_dir, darkmode_bool,
              figure_format= 'both')



plt.show()





#%%
# bee swarm plot with axon length


sbn.set_style('white')
sbn.set_context('paper',font_scale=1)
plt.figure(figsize=get_figure_size(width=100 , height=100))

sbn.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})  


#sbn.violinplot(data=axon_Data, x=None, y=axon_Data['length(µm)'],width=0.1, fill=False)
custom_palette = sbn.color_palette("Dark2",2)
custom_palette_2 = sbn.color_palette("PuRd_r",1)
#filtern von source spalte nach dendritic und somatic
source_filter= ['somatic', 'dendritic']
axon_data_filtered= axon_Data[axon_Data['source'].isin(source_filter)]

#neue spalte hinzufügen= length+dinstance to soma

axon_data_filtered['length + distance']= axon_data_filtered['length(µm)']+axon_data_filtered['distance to soma']

#alle drei parameter in eine spalte packen damit nurn noch zwischen dendritic und somatic unterschieden wird 
axon_data_filtered_melted = pd.melt(axon_data_filtered, id_vars=['source'], value_vars=['length(µm)','distance to soma', 'length + distance'],var_name='variable', value_name='value')


# Violinplots erstellen

sbn.violinplot(x='variable', y='value', hue='source', data=axon_data_filtered_melted, split=True,gap=0.05, inner='quartile', linewidth=0.1, palette=custom_palette)

#swarmplot erstellen
sbn.swarmplot(data=axon_data_filtered_melted,x='variable',  y='value',color='k', alpha=0.4)

#Mittelwert und Standardabweichung berechnen 
means = axon_data_filtered_melted.groupby(['variable', 'source'])['value'].mean().reset_index()
stds = axon_data_filtered_melted.groupby(['variable', 'source'])['value'].std().reset_index()


# Füge Mittelwerte und Standardabweichungen zu den Plots hinzu
for i, variable in enumerate(axon_data_filtered_melted['variable'].unique()):
    for j, source in enumerate(axon_data_filtered_melted['source'].unique()):
        mean_val = means[(means['variable'] == variable) & (means['source'] == source)]['value'].values[0]
        std_val = stds[(stds['variable'] == variable) & (stds['source'] == source)]['value'].values[0]
        plt.errorbar(i + (j - 0.5) * 0.2, mean_val, yerr=std_val, fmt='o', color='black')


plt.show()

#%%
#violin plots for BAOT


#figure style
sbn.set_style('white')
sbn.set_context('paper',font_scale=1)
plt.figure(figsize=get_figure_size(width=100 , height=100))

sbn.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})  


BAOT_palette = [BAOT_color_3,BAOT_color_1,BAOT_color_2]

#filtern von source spalte nach dendritic und somatic
source_filter= ['somatic', 'dendritic']
axon_data_BAOT_filtered= axon_data_BAOT[axon_data_BAOT['source'].isin(source_filter)]

#neue spalte hinzufügen= length+dinstance to soma

axon_data_BAOT_filtered['length + distance']= axon_data_BAOT_filtered['length(µm)']+axon_data_BAOT['distance to soma']

#alle drei parameter in eine spalte packen damit nurn noch zwischen dendritic und somatic unterschieden wird 
axon_data_BAOT_filtered_melted = pd.melt(axon_data_BAOT_filtered, id_vars=['source'], value_vars=['length(µm)','distance to soma', 'length + distance'],var_name='variable', value_name='value')


# Violinplots erstellen

violins = sbn.violinplot(data = axon_data_BAOT_filtered_melted, 
                         x='variable', y='value', hue='source', 
                         split=True, gap=0.05, inner=None,linewidth= 0.1, palette=BAOT_palette)

for violin in violins.collections:
    violin.set_alpha(0.5)
'''
for violin_idx, violin in enumerate(violins.collections):
    
    violin.set_facecolor('None')
    
    if violin_idx % 2:
        violin.set_edgecolor(BAOT_color_3)
    else:
        violin.set_edgecolor(BAOT_color_1)
'''
#swarmplot erstellen
sbn.swarmplot(data=axon_data_BAOT_filtered_melted,x='variable',  y='value', hue = 'source', dodge = True, alpha=1, 
              palette=BAOT_palette,
              size = 4,
              marker = "$\circ$")

#Mittelwert und Standardabweichung berechnen 
means = axon_data_BAOT_filtered_melted.groupby(['variable', 'source'])['value'].mean().reset_index()
stds = axon_data_BAOT_filtered_melted.groupby(['variable', 'source'])['value'].std().reset_index()


# Füge Mittelwerte und Standardabweichungen zu den Plots hinzu
for i, variable in enumerate(axon_data_BAOT_filtered_melted['variable'].unique()):
    for j, source in enumerate(axon_data_BAOT_filtered_melted['source'].unique()):
        mean_val = means[(means['variable'] == variable) & (means['source'] == source)]['value'].values[0]
        std_val = stds[(stds['variable'] == variable) & (stds['source'] == source)]['value'].values[0]
        plt.errorbar(i + (j - 0.5) * 0.2, mean_val, yerr=std_val, fmt='o', color='black')


#%%violin plots for MeA

#figure style
sbn.set_style('white')
sbn.set_context('paper',font_scale=1)
plt.figure(figsize=get_figure_size(width=100 , height=100))

sbn.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})  



#filtern von source spalte nach dendritic und somatic
source_filter= ['somatic', 'dendritic']
axon_data_MeA_filtered= axon_data_MeA[axon_data_MeA['source'].isin(source_filter)]

#neue spalte hinzufügen= length+dinstance to soma

axon_data_MeA_filtered['length + distance']= axon_data_MeA_filtered['length(µm)']+axon_data_MeA['distance to soma']

#alle drei parameter in eine spalte packen damit nurn noch zwischen dendritic und somatic unterschieden wird 
axon_data_MeA_filtered_melted = pd.melt(axon_data_MeA_filtered, id_vars=['source'], value_vars=['length(µm)','distance to soma', 'length + distance'],var_name='variable', value_name='value')


# Violinplots erstellen

sbn.violinplot(x='variable', y='value', hue='source', data=axon_data_MeA_filtered_melted, split=True,gap=0.05, inner='quartile', linewidth=0.1, palette=[MeA_color_3,MeA_color_1,MeA_color_2])

#swarmplot erstellen
#sbn.swarmplot(data=axon_data_MeA_filtered_melted,x='variable',  y='value',color='k', alpha=0.4)

#Mittelwert und Standardabweichung berechnen 
means = axon_data_MeA_filtered_melted.groupby(['variable', 'source'])['value'].mean().reset_index()
stds = axon_data_MeA_filtered_melted.groupby(['variable', 'source'])['value'].std().reset_index()


# Füge Mittelwerte und Standardabweichungen zu den Plots hinzu
for i, variable in enumerate(axon_data_MeA_filtered_melted['variable'].unique()):
    for j, source in enumerate(axon_data_MeA_filtered_melted['source'].unique()):
        mean_val = means[(means['variable'] == variable) & (means['source'] == source)]['value'].values[0]
        std_val = stds[(stds['variable'] == variable) & (stds['source'] == source)]['value'].values[0]
        plt.errorbar(i + (j - 0.5) * 0.2, mean_val, yerr=std_val, fmt='o', color='black')