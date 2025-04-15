# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 18:41:12 2024
Updated on Tue Apr 15 2025

@author: nesseler
"""

# import standard packages
from functions.initialize_packages import *

# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.functions_filter import merge_filter_split_steps
from functions.functions_useful import calc_dvdt_padded
from functions.functions_extractspike import get_AP_parameters, extract_spike

# spike detection
from parameters.parameters import min_peak_prominence, min_peak_distance, min_max_peak_width

# PGF specific
from parameters.PGFs import cc_th1Ap_parameters
t = cc_th1Ap_parameters['t']
SR = cc_th1Ap_parameters['SR']
    

# define protocol
PGF = 'cc_th1AP'
sheet_name = 'PGFs'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = sheet_name)

        
# %% define output

th1AP_rheobase = pd.DataFrame(columns = ['rheobase_abs', 'rheobase_rel'],
                              index = cell_IDs)
th1AP_rheobase.index.name = 'cell_ID'

from parameters.parameters import AP_parameters
rheospike_params = pd.DataFrame(columns = AP_parameters,
                                index = cell_IDs)
rheospike_params.index.name = 'cell_ID'


# %% initialize plotting and verificaiton plots

# init plotting
from functions.initialize_plotting import *

# verification plots
vplots = True
if vplots:
    # load plotting functions
    from analysis.analysis_celldescrip_Syn.plot_analyze_cc_th1AP_syn import plot_full_th1AP, plot_rheospike


# %% load

# cell_IDs = ['E-213']

for cell_ID in tqdm(cell_IDs):

    # load IF protocol
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name = sheet_name)
    
    # get data with file path & trace index
    i, v, _, _, n_steps = get_cc_data(file_path, traceIndex, scale='s')
    
    # filter steps
    v = merge_filter_split_steps(v, SR)
    
    # limit time and voltage to stimulation period
    v_full = np.copy(v)
    t_full = np.copy(t)

    # calc dvdt for all steps
    dvdt_full = np.empty_like(v_full)
    for step in range(n_steps):
        dvdt_full[step] = calc_dvdt_padded(v = v_full[step], t = t_full)
    
    # construct current array
    from functions.functions_constructors import construct_I_array
    i_calc, i_input = construct_I_array(cell_ID, n_steps,
                                        PGF = PGF,
                                        sheet_name = sheet_name,
                                        parameters = cc_th1Ap_parameters)
    
    if vplots:
        plot_full_th1AP(cell_ID, 
                        t = t_full, 
                        v = v_full, 
                        i = i_calc, 
                        i_input = i_input)
        
# %% spike detection

    idc_spikes = []
    
    # limit range for spike detection 
    idc_detection = np.arange(cc_th1Ap_parameters['t_pre'] * (SR / 1e3), cc_th1Ap_parameters['t_pre'] + cc_th1Ap_parameters['t_stim']*2, dtype=int)

    # iterate through steps
    for step in range(n_steps):
        
        # define voltage trace for step
        v_step = v_full[step]
        
        # limit range for spike detection
        v_step_spike = v_step[idc_detection]
        
        # define current 
        i_step = i_input[step]
        
        # find peaks
        idc_spike, dict_spike = sc.signal.find_peaks(v_step, 
                                                     prominence = min_peak_prominence, 
                                                     distance = min_peak_distance * (SR/1e3),
                                                     width = np.multiply(min_max_peak_width, (SR/1e3)))        
        
        # write to list
        if len(idc_spike) > 0:
            idc_spikes.append(idc_spike[0])
        else:
            idc_spikes.append(0)

# %%

    # get rheobase step
    idx_rheo = next(idx for idx, idc_spike in enumerate(idc_spikes) if idc_spike > 0)
    
    # get voltage and dvdt
    v_rheo = v_full[idx_rheo]
    dvdt_rheo = dvdt_full[idx_rheo]
    
    # get input currents
    i_hold = i_calc[0][0]
    i_rheo_abs = np.int64(i_input[idx_rheo])
    i_rheo_rel = np.int64(i_rheo_abs - i_hold)
    
    # write to dataframe
    th1AP_rheobase.loc[cell_ID, ['rheobase_abs', 'rheobase_rel']] = [i_rheo_abs, i_rheo_rel]

    # get AP parameters of first spike
    spike_params, _ = get_AP_parameters(t_spiketrain = t, 
                                        v_spiketrain = v_full[idx_rheo], 
                                        dvdt_spiketrain = dvdt_full[idx_rheo], 
                                        idc_spikes = [idc_spikes[idx_rheo]], 
                                        SR = SR)
    
    # get rheospike
    _, spike_t, spike_v, spike_dvdt = extract_spike(t = t_full, 
                                                    v = v_full[idx_rheo], 
                                                    dvdt = dvdt_full[idx_rheo],
                                                    idx_peak = idc_spikes[idx_rheo])
    
    # write to dataframe
    rheospike_params.loc[cell_ID, :] = spike_params.loc[0, :]
    
    # first spike vplot        
    if vplots:
        plot_rheospike(cell_ID, 
                       t = t_full, 
                       v = v_rheo, 
                       dvdt = dvdt_rheo, 
                       spike_t = spike_t, 
                       spike_v = spike_v, 
                       spike_dvdt = spike_dvdt)
        

# %% saving

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir

# save activity dataframe
th1AP_rheobase.to_excel(join(cell_descrip_dir, 'cc_th1AP-rheobase.xlsx'), 
                      index_label='cell_ID')

rheospike_params.to_excel(join(cell_descrip_dir, 'cc_th1AP-rheobasespike_parameters.xlsx'), 
                      index_label='cell_ID')


# %% update analyzed cells

from functions.update_database import update_analyzed_sheet
    
update_analyzed_sheet(cell_IDs, PGF = PGF)



        