# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:02:40 2023
Updated on Tue Apr 15 2025

@author: moritznesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir
from parameters.parameters import min_peak_prominence, min_peak_distance

# custom functions
from functions.functions_filter import butter_filter
from functions.functions_useful import calc_time_series, calc_dvdt_padded
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol


# define protocol
PGF = 'cc_rest'
sheet_name = 'PGFs'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = sheet_name)

# get number of cells
n_cells = len(cell_IDs)

# define output
activity_df = pd.DataFrame(index = cell_IDs, columns = ['v_rest', 'n_spikes', 't_spikes'])


# %% initialize plotting and verificaiton plots

# init plotting
from functions.initialize_plotting import *

# verification plots
vplots = True
if vplots:
    # load plotting functions
    from analysis.analysis_celldescrip_Syn.plot_analyze_cc_rest_syn import create_cc_rest_vplot
    

# %% data loading

# initialize dataframes to populate in the loop
v_df = pd.DataFrame(columns = cell_IDs)
vf_df = pd.DataFrame(columns = cell_IDs)
SR_df = pd.DataFrame(index = cell_IDs, columns = ['SR'])


print('loading ...')

for cell_idx, cell_ID in enumerate(tqdm(cell_IDs)):
            
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name = sheet_name)
    
    # get data with file path & trace index
    i, v, t, SR, n_step = get_cc_data(file_path, traceIndex, scale='s')
    
    # get first and only step of protocol
    v = v[0]
    
    # edge case when file exceeds the 30 sec recording (E-069)
    if len(v) > (30 * SR):
        warnings.warn(str(cell_ID) + ' exceeds 30 sec and will be cut down.')
        v = v[0:(30*SR)]
        
    # filter all data with 1kHz cutoff
    vf = butter_filter(v, order=3, cutoff=1e3, sampling_rate=SR)
    
    # replace first values with nans to eliminate filter artifact
    vf[:100] = np.nan

    # populate the dataframes & lists
    v_df[cell_ID] = v
    vf_df[cell_ID] = vf
    SR_df.at[cell_ID, 'SR'] = SR
    
    
# %% v_rest

# check if all protocols have the same sampling rate
if len(SR_df['SR'].unique()) != 1:
    warnings.warn('Not all protocols have the same sampling rate!')

else:
    # calculate time
    t = calc_time_series(v, sampling_rate=SR, scale = 's')
    t_ms = calc_time_series(v, sampling_rate=SR, scale = 'ms')


print('calc...')

for cell_ID in tqdm(cell_IDs):

    # define dataframe
    spiketimes_df = pd.DataFrame(index = cell_IDs, columns=['t_spikes'])
    
    # get filtered voltage and sampling rate
    t_ms = t_ms
    vf = vf_df[cell_ID]
    SR = SR_df.at[cell_ID, 'SR']
    
    # find peaks
    idc_spikes, dict_peak = sc.signal.find_peaks(vf, 
                                                 prominence = min_peak_prominence, 
                                                 distance = min_peak_distance * (SR/1e3))
        
    # calculate spike times in seconds
    t_spikes = np.divide(idc_spikes, SR)
    n_spikes = len(t_spikes)
    
    # write to dataframe
    activity_df.at[cell_ID, 't_spikes'] = list(t_spikes)
    activity_df.at[cell_ID, 'n_spikes'] = n_spikes
    
    
    ### v_rest
    
    from functions.functions_extractspike import extract_spike
    
    # cut out spikes for resting membrane potential calc
    if n_spikes > 0:
        
        # calc first derivative
        dvdt = calc_dvdt_padded(vf, t_ms)
        
        # define negative dvdt threshold, i.e.: detection of fast AHP
        dvdt_n_threshold = -1
        
        # copy voltage trace to adjustable variable 
        vf_wo_spikes = np.copy(vf)
        
        # iterate through spikes
        for spike_idx in idc_spikes:
        
            spike_idc, _, _, _ = extract_spike(t = t_ms, 
                                               v = vf_df[cell_ID].to_numpy(), 
                                               dvdt = dvdt, 
                                               idx_peak = spike_idx)
        
            # replace spike values with nans
            vf_wo_spikes[spike_idc] = np.nan

        # calc v_rest as mean over trace
        v_rest = np.nanmean(vf_wo_spikes)
        
    else:
        # calc v_rest as mean over trace
        v_rest = np.mean(vf)
        
        # set vf_wo_spikes
        vf_wo_spikes = np.full_like(vf, np.nan)
        
    # write to dataframe
    activity_df.at[cell_ID, 'v_rest'] = v_rest
    
            
    # create verification plot
    if vplots:
        create_cc_rest_vplot(cell_ID, 
                             t, 
                             vf, 
                             vf_wo_spikes, 
                             t_spikes, 
                             v_rest)


# %% create dict with all, active and non-active cells

# add activity column for categorical plots´with silent as default value
activity_df['activity'] = 'silent'

# change activity value of spiking cells with n_spike > 0 to 'spiking'
activity_df.loc[activity_df['n_spikes'] > 0, 'activity'] = 'spiking'


# %% saving

# save activity dataframe
activity_df.to_excel(join(cell_descrip_dir, 'cc_rest-activity.xlsx'), 
                      index_label='cell_ID')


# %% update analyzed cells

from functions.update_database import update_analyzed_sheet
    
update_analyzed_sheet(cell_IDs, PGF = PGF)



    

