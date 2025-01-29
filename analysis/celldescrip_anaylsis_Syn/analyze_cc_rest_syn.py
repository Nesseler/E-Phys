# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:58:26 2025

@author: nesseler

Output files : activity_df.xlsx -> cell_descriptors_Syn

"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_syn_dir
from parameters.parameters import min_peak_prominence, min_peak_distance

# custom functions
from functions.functions_useful import butter_filter, calc_time_series, calc_dvdt_padded
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.get_cell_IDs import get_cell_IDs_one_protocol


# define protocol
PGF = 'cc_rest'
sheet_name = 'PGFs_Syn'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = sheet_name)

# get number of cells
n_cells = len(cell_IDs)

# init plotting
from functions.initialize_plotting import *

# define output
activity_df = pd.DataFrame(index = cell_IDs, columns = ['v_rest', 'n_spikes', 't_spikes'])


# %%

# load anaylsis worksheet
from parameters.directories_win import table_file
analyzed = pd.read_excel(table_file,
                         sheet_name = 'analyzed',
                         index_col = 'cell_ID')

# get list of cell_IDs already analyzed
analyzed_cell_IDs = analyzed.loc[analyzed[PGF].notna()][PGF].index.to_list()

# redefine cell_IDs list
cell_IDs = [cell_ID for cell_ID in cell_IDs if cell_ID not in analyzed_cell_IDs]

# raise error
if len(cell_IDs) == 0:
    raise ValueError('Nothing new to analyze!')


# %%
# verification plots
vplots = True

if vplots:
    
    # load plotting functions
    from analysis.celldescrip_anaylsis_Syn.plot_analyze_cc_rest_syn import create_cc_rest_vplot

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
    
# %%

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
                                               idx_peak = spike_idx,
                                               dvdt_n_threshold=dvdt_n_threshold)
        
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


export_filename = 'cc_rest-Syn-activity.xlsx'

try:
    loaded_activity_df = pd.read_excel(join(cell_descrip_syn_dir, export_filename),
                                       index_col = 'cell_ID')
    
    # find out how to combine both dataframes
    loaded_activity_df.loc[cell_IDs] = activity_df.loc[cell_IDs].values
    
    # save activity dataframe
    loaded_activity_df.to_excel(join(cell_descrip_syn_dir, export_filename), 
                                index_label='cell_ID')
    
except FileNotFoundError:
    # save activity dataframe
    activity_df.to_excel(join(cell_descrip_syn_dir, export_filename), 
                         index_label='cell_ID')



# %% update analyzed cells

from functions.update_database import update_analyzed_sheet
    

update_analyzed_sheet(cell_IDs, PGF = PGF)



print('Finished!')
    