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


# %%
def plot_half_violin(data, ax,
                      data_pad_factor = 1,
                      v_resolution = np.nan,
                      v_kde_cutoff = 0.01,
                      v_abs_cutoff = [np.nan, np.nan],
                      v_bandwidth = np.nan,
                      v_position = 0,
                      v_direction = -1,
                      v_offset = -0.05,
                      v_width = 0.25,
                      v_color = 'w',
                      v_lw = 1,
                      v_baseline = False,
                      v_fill = False,
                      v_fillcolor = 'w',
                      v_filllw = 1,
                      v_zorder = 0):
    
    """
    Function plots half violin at specified position.
    Parameters:
        data : Array of data that is used for plotting.
        ax : matplotlib axis to plot violin on.
        data_pad_factor : Float, default is 1. Factor that data_range is multiplied
                          with for estimation of kernel density.
        
        v_resolution : Float (data scale), default is np.nan. Resolution of violin in 
                        data native scale.
        v_kde_cutoff : Float (kde percentage), default is 0.01. Cutoff that is used
                        to limit violin.
        v_abs_cutoff : List of floats (data scale), default is [np.nan, np.nan]. 
                        Absolute values that can be used to limit violin. Shape of 
                        list is [Min, Max]. When only one value is specified, the 
                        other must be set to np.nan and here the v_kde_cutoff is used.
                        If both values are specified the kde_cutoff must be np.nan to 
                        take effect.
        v_bandwidth : Float (0-1), default is np.nan. Sets bandwidth of density.
        v_position : Float, default is 0. Position of violin on x axis.
        v_direction : -1 or 1, default is -1. Direction factor of violin on x
                      axis.
        v_offset : Float, default is -0.05. Offset of violin from v_position.
        v_width : Float, default is 0.25. Width of violin on x axis
        v_color : Str, default is 'w'. Color of violin.
        v_lw : Float, default is 1. Linewidth of violin.
        v_baseline : Boolean, default is False. Boolean to plot the baseline of
                      the violin.
        v_fill : Boolean, default is False. Boolean to plot the fill of the violin.
        v_fillcolor : Str, default is 'w'. Color of fill.
        v_filllw : Float, default is 1. Linewidth of fill.
        v_zorder : Integer, default is 0. Z-Order of violin.
    
    Returns:
        half_violin : List of Line and Fillobjects. Order is [violin, fill] 
    """
    
    from scipy.stats import gaussian_kde
    
    # # function inputs
    # data = [5,60,85,100,140,-10,15,20,35,25,15,165,85,75,80,165,15,15,100,25,45,5]
    # ax = plt.gca()
    
    # # defaults
    # data_pad_factor = 1
    # v_resolution = 0.1 # in native data scale
    # v_kde_cutoff = 0.01 # in % default is 0.01
    # v_abs_cutoff = [-50, np.nan] # in native data scale, list of [min, max] for violin cutoff
    # v_bandwidth = np.nan
    
    # v_position = 0
    # v_direction = -1
    # v_offset = -0.05
    # v_width = 1
    
    # v_color = 'w'
    # v_lw = 1
    # v_baseline = True
    # v_fill = False
    # v_fillcolor = v_color
    # v_filllw = v_lw
    # v_zorder = 0
    
    # get distribution metrics
    data_n = len(data)
    data_min = np.min(data)
    data_max = np.max(data)
    data_range = data_max - data_min
    
    # pad data range for violin that extends over range
    padded_min = data_min - (data_range * data_pad_factor)
    padded_max = data_max + (data_range * data_pad_factor)
    
    # get kernel density estimate
    density = gaussian_kde(data)
    
    # scale factor (bandwidth)
    # bandwidth is calculated by n**(-1./(d+4)), with n datapoints and d dimensions
    # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde
    if not np.isnan(v_bandwidth):
        density.set_bandwidth(v_bandwidth)
        
    # calc resolution of violin with 1000 points between min and max
    if np.isnan(v_resolution):
        v_resolution = data_range / 1000

    # get range of value for density estimation
    vs = np.arange(padded_min, padded_max + v_resolution, v_resolution)
    
    # get kde (kernel density esimate)
    kde = density(vs)
    
    # normalize kde to maximum
    kde_normed = [(k / np.max(kde)) for k in kde]
    
    # # # violin cutoff # # #
    
    # condition for set kde cutoff
    if (not np.isnan(v_kde_cutoff) and np.isnan(v_abs_cutoff).all()):
        # remove values below 1 % on edges
        kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff or v < data_max and v > data_min)]
        # print('0')
        
    # condition for set min and max absolute cutoff
    elif (np.isnan(v_kde_cutoff) and not np.isnan(v_abs_cutoff).all()):
        # remove values outside absolute cutoff
        kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (v > v_abs_cutoff[0] and v < v_abs_cutoff[1])]
        # print('1')
      
    # condition for set kde_cutoff and either abs min or max cutoff
    elif (not np.isnan(v_kde_cutoff) and np.isnan(v_abs_cutoff).any()):
        # condition for set min
        if not np.isnan(v_abs_cutoff[0]):
            # remove values below 1 % on edges
            kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff and v > v_abs_cutoff[0])]
            # print('2.1')
            
        # condition for set max 
        elif not np.isnan(v_abs_cutoff[1]):
            # remove values below 1 % on edges
            kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff and v < v_abs_cutoff[1])]
            # print('2.2')
        
    # redefine violin
    kde_clipped = [kv[0] for kv in kv_s]
    vs_clipped = [kv[1] for kv in kv_s]
    
    # add direction, scale, offset ,and x position
    kde_withOffset = [v_direction * k * v_width + v_offset + v_position for k in kde_clipped]
    
    # set coordinates for baseline
    vs_baseline = vs_clipped
    ks_baseline = [v_offset + v_position] * len(vs_clipped)
    
    if not v_baseline:
        ks = kde_withOffset  # x
        vs = vs_clipped      # y
    
    
    elif v_baseline:
        # concatenate two x and y arrays
        halflen_baseline = int(len(vs_baseline) / 2)
        
        vs = vs_baseline[0:halflen_baseline][::-1] + vs_clipped + vs_baseline[halflen_baseline+1:-1][::-1]
        ks = ks_baseline[0:halflen_baseline][::-1] + kde_withOffset + ks_baseline[halflen_baseline+1:-1][::-1]
        
    
    # plot half violin
    half_violin_list = ax.plot(ks, vs,
                               color = v_color,
                               linewidth = v_lw,
                               label = '_nolegend_',
                               zorder = 2)
    
    # violin return objects
    half_violin = half_violin_list
    
    
    
    if v_fill:
        v_face = ax.fill_betweenx(y = vs_clipped,
                                  x1 = kde_withOffset,
                                  x2 = ks_baseline,
                                  color = v_fillcolor,
                                  linewidth = v_filllw,
                                  zorder = v_zorder)
        
        # append to return object
        half_violin.append(v_face)
    

    return half_violin

# %%

# initialize plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *


fig, axs = plt.subplots(nrows = 1,
                        ncols = 2)

for v_idx, variable in enumerate(['length(µm)', 'distance to soma', 'length + distance']):
    
    for (r_idx, region), direction in zip(enumerate(['MeA', 'BAOT']), [-1, 1]):
        
        print(v_idx, variable, region, r_idx, direction)
        
        data = axon_data_filtered[axon_data_filtered['Region'] == region].loc[:, variable]
        
        x = r_idx + (v_idx *2)
        print(x)
    
        plot_half_violin(data = data,
                         ax = axs[0],
                         v_position = x,
                         v_direction = direction,
                         v_offset = -direction * 0.4,
                         v_width = 0.75,
                         v_color = region_colors[region],
                         v_zorder = 2,
                         v_abs_cutoff = [0, np.nan])


axs[0].set_xticks(ticks = np.arange(0, 6, 1))

plt.show()




