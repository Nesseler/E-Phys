# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 18:32:53 2024

@author: nesseler (aber eigentlich büsching) [das will ich sehen]
"""

###venn_diagram neus

import pandas as pd
import seaborn as sbn
import numpy as np
import matplotlib.pyplot as plt
from os.path import join, dirname
from matplotlib_venn import venn2




# %% functions

def get_onlyfiles_list(dir_path):

    from os import listdir
    from os.path import isfile, join  

    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]


def get_cell_IDs_all_listofPGFs(table_file, PGF_names = ['cc_rest, cc_IF']):
    '''
    This function returns a list of cell IDs that were stimulated with all 
    frequencies (1 - 75 Hz) in the cc_APs protocols.
    Parameter:
        -
    Returns:
        cell_IDs : List of cell_IDs
    '''
    # get table
    table = pd.read_excel(table_file, sheet_name="PGFs", index_col='cell_ID')
    
    # loop to create string to include all frequencies in query
    query_str = ''
    
    
    for idx, PGF in enumerate(PGF_names):

        if idx > 0:
            query_str = query_str + ' and '
            
        query_str = query_str + f'{PGF}.notnull()'   
    
    # limit lookup table
    lookup_table = table.query(query_str)
    
    # cell IDs 
    cell_IDs = list(lookup_table.index)

    return cell_IDs



# %% script

from parameters.directories_win import table_file

from functions.functions_plotting import get_colors, set_font_sizes, get_figure_size
# table_file = '//Fileserver/AG Spehr BigData/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/' + 'ePhys-database.xlsx'

# load all MetaData
MetaData = pd.read_excel(table_file,
                         sheet_name="MetaData",
                         index_col='cell_ID')

# load sheet with PGF references
PGFs = pd.read_excel(table_file,
                     sheet_name="PGFs",
                     index_col='cell_ID')

# get all cell_IDs with all PGFs
cell_IDs_ePhys = get_cell_IDs_all_listofPGFs(table_file, ['cc_rest', 'cc_IF'])#, 'cc_th1AP', 'cc_APs_1Hz', 'cc_APs_5Hz', 'cc_APs_10Hz', 'cc_APs_30Hz', 'cc_APs_50Hz', 'cc_APs_75Hz'])

# turn into dataframe with cell_IDs as index
ePhys_bool = pd.DataFrame({'all_ePhys' : np.ones(len(cell_IDs_ePhys))}, index=cell_IDs_ePhys)


# get all reconstructed cells from MetaData as booleans
recon_bool = MetaData[MetaData['reconstructed'] == 1]['reconstructed']

# concatenate both for later handling
venn_bool_df = pd.concat([recon_bool, ePhys_bool], axis = 1)
venn_bool_df.sort_index(inplace = True)

venn_bool_df_cell_IDs = venn_bool_df.index.to_list()

venn_bool_df = pd.concat([venn_bool_df, MetaData.loc[venn_bool_df_cell_IDs, 'Region']], axis = 1)


def get_set_numbers(venn_bool_df, region = 'all'):
    
    if region != 'all':
        venn_bool_df = venn_bool_df[venn_bool_df['Region'] == region]
        
    # get number of reconstructed and ePhys cells
    n_morph_n_ephys = len(venn_bool_df.query("all_ePhys.notnull() and reconstructed.notnull()"))
    
    n_morph_only =  len(venn_bool_df.query("reconstructed.notnull()")) - n_morph_n_ephys
    n_ephys_only =  len(venn_bool_df.query("all_ePhys.notnull()")) - n_morph_n_ephys
        
                                             
    return n_morph_only, n_ephys_only, n_morph_n_ephys


# BAOT_color = '#7a66fc' #'#ff1b6b' #'#03C03C'
# MeA_color =  '#ff8d00' #'#45caff' #'#FF8F00'

# region_colors = {'BAOT' : BAOT_color,
#                  'MeA' : MeA_color,
#                  'BAOT/MeA' : 'gray'}

darkmode_bool = True

colors_dict, region_colors = get_colors(darkmode_bool)

# %% create venn diagramm for all

fig_venns, axs_venns = plt.subplots(nrows = 2,
                                    ncols = 2,
                                    layout = 'constrained')

axs_venns = axs_venns.flatten()

axs_venns[0].set_title('all')
venn_all = venn2(subsets=(get_set_numbers(venn_bool_df)), 
                 set_labels = ('reconstructed', 'all_ePhys'), 
                 ax = axs_venns[0])

axs_venns[1].set_title('BAOT')
venn_BAOT = venn2(subsets=(get_set_numbers(venn_bool_df, region = 'BAOT')), 
                  set_labels = ('reconstructed', 'all_ePhys'), 
                  ax = axs_venns[1])            

axs_venns[2].set_title('MeA')
venn_MeA = venn2(subsets=(get_set_numbers(venn_bool_df, region = 'MeA')), 
                 set_labels = ('reconstructed', 'all_ePhys'), 
                 ax = axs_venns[2])

axs_venns[3].set_title('BAOT/MeA')
venn_BAOTnMeA = venn2(subsets=(get_set_numbers(venn_bool_df, region = 'BAOT/MeA')), 
                      set_labels = ('reconstructed', 'all_ePhys'), 
                      ax = axs_venns[3])                                          
                                             
plt.show()                                        
                                             
                                             
# %% alternative venn diagramm figure

fig_venns, axs_venns = plt.subplots(nrows = 2,
                                    ncols = 3,
                                    layout = 'constrained')



axs_venns[0][0].set_title('all')
venn_all = venn2(subsets=(get_set_numbers(venn_bool_df)),set_labels = ('reconst.', 'ePhys'),
                 ax = axs_venns[0][0], set_colors=('darkgray', 'gray'), alpha=0.8)

axs_venns[1][0].set_title('BAOT')
venn_BAOT = venn2(subsets=(get_set_numbers(venn_bool_df, region = 'BAOT')),set_labels = ('reconst.', 'ePhys'),
                  ax = axs_venns[1][0], set_colors=('darkblue', 'lightblue'), alpha=0.8)            

axs_venns[1][1].set_title('MeA')
venn_MeA = venn2(subsets=(get_set_numbers(venn_bool_df, region = 'MeA')),set_labels = ('reconst.', 'ePhys'), 
                 ax = axs_venns[1][1],  set_colors=('darkorange', 'gold'), alpha=0.8)

axs_venns[1][2].set_title('BAOT/MeA')
venn_BAOTnMeA = venn2(subsets=(get_set_numbers(venn_bool_df, region = 'BAOT/MeA')), set_labels = ('reconst.', 'ePhys'), 
                      ax = axs_venns[1][2],  set_colors=('darkgray', 'gray'), alpha=0.8)       
                                   
plt.delaxes(ax = axs_venns[0][2])
plt.delaxes(ax = axs_venns[0][1])             
                                    
plt.show()                                                                                    
                                             
                                             
                                             
                                             
                                             
                                             


