# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 09:31:49 2025

@author: nesseler
"""

from functions.initialize_packages import *


# %% load dataframes to describe cells

from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_metrics_dir

height_width_depth = pd.read_excel(join(cellmorph_metrics_dir, 'height_width_depth' + '.xlsx'), index_col = 'cell_ID')
n_primary = pd.read_excel(join(cellmorph_metrics_dir, 'n_primary' + '.xlsx'), index_col = 'cell_ID')
n_terminal = pd.read_excel(join(cellmorph_metrics_dir, 'n_terminal' + '.xlsx'), index_col = 'cell_ID')
bifurcation_ratio = pd.read_excel(join(cellmorph_metrics_dir, 'bifurcation_ratio' + '.xlsx'), index_col = 'cell_ID')
total_cable_length = pd.read_excel(join(cellmorph_metrics_dir, 'total_cable_length' + '.xlsx'), index_col = 'cell_ID')
sholl_metrics = pd.read_excel(join(cellmorph_metrics_dir, 'sholl_metrics' + '.xlsx'), index_col = 'cell_ID')
axon_data = pd.read_excel(join(cellmorph_metrics_dir, 'axon_data' + '.xlsx'), index_col = 'cell_ID')
spines = pd.read_excel(join(cellmorph_metrics_dir, 'spines' + '.xlsx'), index_col = 'cell_ID')
circ_stats = pd.read_excel(join(cellmorph_metrics_dir, 'circ_stats' + '.xlsx'), index_col = 'cell_ID')


# %% get cell_IDs & MetaData

# get cell_IDs
cell_IDs = height_width_depth.index.to_list()

# limit 
cell_IDs = cell_IDs[:-20]

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# reindex axon data 
axon_data = axon_data.reindex(cell_IDs, fill_value = np.nan)
spines = spines.reindex(cell_IDs, fill_value = np.nan)

# %% create dataframes

cellmetrics = pd.concat([height_width_depth.loc[cell_IDs, [f'{ntype}-{metric}' for ntype in ['dendrites', 'axons'] for metric in ['height', 'width']]],
                         n_primary.loc[cell_IDs, [f'n_primary-{ntype}' for ntype in ['dendrites']]],
                         n_terminal.loc[cell_IDs, [f'n_terminal-{ntype}' for ntype in ['dendrites']]],
                         bifurcation_ratio.loc[cell_IDs, [f'bifurcation_ratio-{ntype}' for ntype in ['dendrites', 'axons']]],
                         total_cable_length.loc[cell_IDs, [f'total_cable_length-{ntype}' for ntype in ['dendrites', 'axons']]],
                         sholl_metrics.loc[cell_IDs, [f'{metric}-{ntype}' for ntype in ['dendrites', 'axons'] for metric in ['max_intersections', 'critical_radius', 'enclosing_radius']]],
                         axon_data.loc[cell_IDs, ['distance to soma']]],
                         axis = 1)

cellmetrics_aux = pd.concat([MetaData.loc[cell_IDs, 'Region'],
                             spines.loc[cell_IDs, 'Spines'],
                             circ_stats.loc[cell_IDs, [f'circmean_rad-{ntype}' for ntype in ['dendrites', 'axons']]],
                             axon_data.loc[cell_IDs, 'source']
                            ],
                            axis = 1)

# %% rename parameters

for ntype in ['dendrites', 'axons']:
    
    # rename columns
    cellmetrics.rename(columns = {'n_primary-dendrites'        : 'dendrites-n_primary',
                                  'n_terminal-dendrites'       : 'dendrites-n_terminal',
                                  f'bifurcation_ratio-{ntype}' : f'{ntype}-bifurcation_ratio',
                                  f'total_cable_length-{ntype}': f'{ntype}-total_cable_length',
                                  f'critical_radius-{ntype}'   : f'{ntype}-critical_radius', 
                                  f'enclosing_radius-{ntype}'  : f'{ntype}-enclosing_radius', 
                                  f'max_intersections-{ntype}' : f'{ntype}-max_intersections'},
                       inplace = True)
    
    # rename columns
    cellmetrics_aux.rename(columns = {f'circmean_rad-{ntype}'         : f'{ntype}-circular_mean'},
                           inplace = True)
    
# rename columns
cellmetrics.rename(columns = {'distance to soma' : 'AIS distance to soma'}, inplace = True)
cellmetrics_aux.rename(columns = {'Spines' : 'spinyness', 'source' : 'AIS root'}, inplace = True)


# %% sort columns

sorter = ['dendrites-width',
          'dendrites-height', 
          'dendrites-n_primary',
          'dendrites-n_terminal',
          'dendrites-bifurcation_ratio', 
          'dendrites-total_cable_length', 
          'dendrites-critical_radius', 
          'dendrites-enclosing_radius', 
          'dendrites-max_intersections',
          'axons-width', 
          'axons-height',
          'axons-bifurcation_ratio', 
          'axons-total_cable_length', 
          'axons-critical_radius', 
          'axons-enclosing_radius', 
          'axons-max_intersections',
          'AIS distance to soma']

cellmetrics = cellmetrics[sorter]

# %% drop axon-less cells

cellmetrics = cellmetrics.dropna(axis = 0, how = 'any')
cellmetrics_aux = cellmetrics_aux.dropna(axis = 0, how = 'any')


# %% save dataframe

# from parameters.directories_win import clustering_dir
cellmetrics.to_excel(join(cellmorph_metrics_dir, 'cellmorph_cellmetrics.xlsx'), index_label = 'cell_ID')
cellmetrics_aux.to_excel(join(cellmorph_metrics_dir, 'cellmorph_cellmetrics_auxillary.xlsx'), index_label = 'cell_ID')