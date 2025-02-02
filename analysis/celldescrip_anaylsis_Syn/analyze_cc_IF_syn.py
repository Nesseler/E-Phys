# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:27:42 2025

@author: nesseler
"""
# import standard packages
from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_syn_dir

# spike detection
from parameters.parameters import min_peak_prominence_ccIF, min_peak_distance_ccIF, min_max_peak_width_ccIF

# IF characteristation
from parameters.parameters import t_init_inst_freq


# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.functions_filter import merge_filter_split_steps
from functions.functions_useful import calc_dvdt_padded
from functions.functions_extractspike import get_AP_parameters, extract_spike, get_spiketrain_n_ISI_parameter

# PGF specific
from parameters.PGFs import cc_IF_syn_parameters
t = cc_IF_syn_parameters['t']
SR = cc_IF_syn_parameters['SR']
    

# define protocol
PGF = 'cc_IF'
sheet_name = 'PGFs_Syn'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = sheet_name)

        
# %% define output

IF           = pd.DataFrame(columns=cell_IDs, index = np.arange(-50, 1000, 1, int))
IF_inst      = pd.DataFrame(columns=cell_IDs, index = np.arange(-50, 1000, 1, int))
IF_inst_init = pd.DataFrame(columns=cell_IDs, index = np.arange(-50, 1000, 1, int))

# set index label
for df in [IF, IF_inst, IF_inst_init]:
    df.index.name = 'i_input'

IF_dict      = pd.DataFrame(columns = ['i_rheobase'       , 'idx_rheobase',
                                       'i_maxfreq'        , 'idx_maxfreq',
                                       'i_halfmax'        , 'idx_halfmax',
                                       'i_maxinitinstfreq', 'idx_maxinitinstfreq'],
                            index = cell_IDs)

IF_rheobase = pd.DataFrame(columns = ['rheobase_abs', 'rheobase_rel', 'v_thres_rheobase_spike'], 
                           index = cell_IDs)

adaptation = pd.DataFrame(index = cell_IDs)

# set index label
for df in [IF_dict, IF_rheobase, adaptation]:
    df.index.name = 'cell_IDs'


# %% check for new cells to be analyzed

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


# %% initialize plotting and verificaiton plots

# init plotting
from functions.initialize_plotting import *

# verification plots
vplots = True
if vplots:
    # load plotting functions
    from analysis.celldescrip_anaylsis_Syn.plot_analyze_cc_IF_syn import plot_full_IF, plot_IF_step_spike_detection, plot_rheobase, plot_adaptation

    
# %% load

# cell_IDs = ['E-300']

for cell_ID in tqdm(cell_IDs):

    # load IF protocol
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name = sheet_name)
    
    # get data with file path & trace index
    i, v, _, _, n_steps = get_cc_data(file_path, traceIndex, scale='s')
    
    # filter steps
    v = merge_filter_split_steps(v, SR)
    
    # limit time and voltage to stimulation period
    idc_stim = cc_IF_syn_parameters['idc_stim']
    v_full = np.copy(v)
    t_full = np.copy(t)
    v_stim = v_full[:, idc_stim]
    t_stim = t_full[idc_stim]
    
    # construct current array
    from functions.functions_constructors import construct_I_array
    i_calc, i_input = construct_I_array(cell_ID, n_steps,
                                        PGF = PGF,
                                        sheet_name = sheet_name,
                                        parameters = cc_IF_syn_parameters)
    
    if vplots:
        plot_full_IF(cell_ID, 
                     t = t_full, 
                     v = v_full, 
                     i = i_calc, 
                     i_input = i_input)
        

# %% spike detection

    # iterate through steps
    
    for step in range(n_steps):
        
        # define voltage trace for step
        v_step = v_stim[step]
        
        # define current 
        i_step = i_input[step]
        
        # find peaks
        idc_spikes, dict_spikes = sc.signal.find_peaks(v_step, 
                                                       prominence = min_peak_prominence_ccIF, 
                                                       distance = min_peak_distance_ccIF * (SR/1e3),
                                                       width = np.multiply(min_max_peak_width_ccIF, (SR/1e3)))
            
        # calculate spike times in seconds
        t_spikes = np.divide(idc_spikes, (SR/1e3))
        n_spikes = len(t_spikes)
        
        # calc freq as number of APs over 1000 ms
        freq = n_spikes / (cc_IF_syn_parameters['t_stim'] / 1e3)
    
           
        # calculate ISI
        if n_spikes >= 2:
            
            # get ISI as difference of spike times
            ISIs = np.diff(t_spikes)
            
            # calculate average ISI
            mean_ISI = np.mean(ISIs)
    
            # calculate the instantaneous firing frequency as inverse of the firing frequency
            inst_freq = (1 / mean_ISI ) * 1e3
            
            # calculate the instantaneous initial firing frequency
            # as the mean firing frequency in the first (100) ms
            init_spikes = [t_spike for t_spike in t_spikes if t_spike <= t_init_inst_freq]
            
            # check for spikes in initial time frame
            if len(init_spikes) > 1:
                # calc ISIs
                init_ISIs = np.diff(init_spikes)
            
                # calc mean
                mean_init_ISIs = np.mean(init_ISIs)
            
                # calc freq
                init_inst_freq = (1 / mean_init_ISIs ) * 1e3
            
            else:
                # set default value
                init_inst_freq = np.nan 
    
        else:
            # set default value
            inst_freq = np.nan
            init_inst_freq = np.nan 
        
        # write to dataframe
        IF.at[i_step, cell_ID] = freq
        IF_inst.at[i_step, cell_ID] = inst_freq
        IF_inst_init.at[i_step, cell_ID] = init_inst_freq
        
        if True:
            plot_IF_step_spike_detection(cell_ID,
                                         step, 
                                         t = t_stim,
                                         v = v_step, 
                                         t_spikes = t_spikes,
                                         freq = freq,
                                         inst_freq = inst_freq,
                                         init_inst_freq = init_inst_freq)
    
    # %% rheobase
    
    # get rheobase step
    idx_rheo = next(idx for idx, n_spikes in enumerate(IF[cell_ID].dropna()) if n_spikes > 0)
    
    # define voltage trace
    v_rheo = v_full[idx_rheo]
    
    # calc first derivative
    dvdt_rheo = calc_dvdt_padded(v_rheo, t_full)
    
    # get rheobase current
    i_hold = i_input[0]  
    i_rheo_abs = i_input[idx_rheo]
    i_rheo_rel = i_input[idx_rheo] - i_hold
    
    # find rheobase spike
    idc_spikes, dict_spikes = sc.signal.find_peaks(v_rheo, 
                                                   prominence = min_peak_prominence_ccIF, 
                                                   distance = min_peak_distance_ccIF * (SR/1e3),
                                                   width = np.multiply(min_max_peak_width_ccIF, (SR/1e3)))
    
    # get rheospike
    _, rheospike_t, rheospike_v, rheospike_dvdt = extract_spike(t_full, v_rheo, dvdt_rheo, idc_spikes[0])
    
    # measure the rheobase spike
    rheospike_params, _ = get_AP_parameters(t_full, v_rheo, dvdt_rheo, [idc_spikes[0]])
    
    if True:
        plot_rheobase(cell_ID, 
                      idx_rheo,
                      t = t_full,
                      v = v_rheo,
                      dvdt = dvdt_rheo,
                      rheospike_t = rheospike_t, 
                      rheospike_v = rheospike_v, 
                      rheospike_dvdt = rheospike_dvdt, 
                      rheospike_params = rheospike_params)
    
    # write to dataframe
    IF_rheobase.loc[cell_ID, ['rheobase_abs', 'rheobase_rel', 'v_thres_rheobase_spike']] = [i_rheo_abs, i_rheo_rel, rheospike_params.loc[0, 'v_threshold']]
    
    
    # %% adaptation
    
    # max freq
    idx_maxfreq  = IF[cell_ID].dropna().argmax()
    v_maxfreq    = v_stim[idx_maxfreq]
    dvdt_maxfreq = calc_dvdt_padded(v_maxfreq, t_stim)
    i_maxfreq    = i_input[idx_maxfreq]


    # max initial instant frequency    
    # edge-case: if initial instant freq has only one value
    if IF_inst_init[cell_ID].dropna().shape[0] == 1:
        i_maxinitinstfreq    = IF_inst_init[cell_ID].dropna().index.values[0]
        idx_maxinitinstfreq  = np.where(IF[cell_ID].dropna().index == i_maxinitinstfreq)[0][0]
    else:
        i_maxinitinstfreq   = IF_inst_init[cell_ID].dropna().index[IF_inst_init[cell_ID].dropna().argmax()]
        idx_maxinitinstfreq = np.where((IF[cell_ID].dropna().index == i_maxinitinstfreq) == True)[0][0]
        
    v_maxinitinstfreq    = v_stim[idx_maxinitinstfreq]
    dvdt_maxinitinstfreq = calc_dvdt_padded(v_maxinitinstfreq, t_stim)
    
    
    # step half-way between rheobase and max inst initial freq
    idx_halfmax  = int(((idx_maxinitinstfreq - idx_rheo) / 2) + idx_rheo)
    v_halfmax    = v_stim[idx_halfmax]
    dvdt_halfmax = calc_dvdt_padded(v_halfmax, t_stim)
    i_halfmax    = i_input[idx_halfmax]
    
    # write to dataframe
    IF_dict.at[cell_ID, 'i_rheobase']          = i_rheo_abs - i_hold
    IF_dict.at[cell_ID, 'idx_rheobase']        = idx_rheo
    IF_dict.at[cell_ID, 'i_maxfreq']           = i_maxfreq - i_hold
    IF_dict.at[cell_ID, 'idx_maxfreq']         = idx_maxfreq
    IF_dict.at[cell_ID, 'i_halfmax']           = i_halfmax - i_hold
    IF_dict.at[cell_ID, 'idx_halfmax']         = idx_halfmax
    IF_dict.at[cell_ID, 'i_maxinitinstfreq']   = i_maxinitinstfreq - i_hold
    IF_dict.at[cell_ID, 'idx_maxinitinstfreq'] = idx_maxinitinstfreq

        
    from functions.functions_extractspike import get_spiketrain_n_ISI_parameter
    
    # get all spike parameters
    rheobase_spikes, rheobase_ISIs               = get_spiketrain_n_ISI_parameter(t_stim, v_stim[idx_rheo], dvdt_rheo[idc_stim], SR)
    maxfreq_spikes, maxfreq_ISIs                 = get_spiketrain_n_ISI_parameter(t_stim, v_maxfreq, dvdt_maxfreq, SR)
    halfmax_spikes, halfmax_ISIs                 = get_spiketrain_n_ISI_parameter(t_stim, v_halfmax, dvdt_halfmax, SR)
    maxinitinstfreq_spikes, maxinitinstfreq_ISIs = get_spiketrain_n_ISI_parameter(t_stim, v_maxinitinstfreq, dvdt_maxinitinstfreq, SR)
    
        
    from parameters.parameters import adaptation_n_lastspikes, adaptation_popt_guess_linear_ISIs
    from functions.functions_useful import linear_func
    
    # set step used for adapatation comparisons
    # as step with highest number of spikes
    adaptation_spikes = maxfreq_spikes
    adaptation_ISIs = maxfreq_ISIs
    
    # get number of spikes in max freq step
    n_spikes_maxfreq = adaptation_spikes.shape[0]
    
    #
    if n_spikes_maxfreq < adaptation_n_lastspikes - 1:
        n_lastspikes = n_spikes_maxfreq -1
    else:
        n_lastspikes = adaptation_n_lastspikes
    
    # # # spike amplitude adaptation # # #
    
    # get amplitude of first spike
    fst_spike_vamplitude = adaptation_spikes.at[0, 'v_amplitude']
    
    # get amplitude of last 4 spikes & calc mean
    lst_spike_vamplitude = adaptation_spikes.tail(n_lastspikes)['v_amplitude'].mean()
    
    # calc ratio of first to last spikes
    spike_amplitude_adaptation = lst_spike_vamplitude / fst_spike_vamplitude
    
    
    # # # spike width (FWHM) adatation # # #
    
    # compare first and mean of last few spikes
    fst_spike_FWHM = adaptation_spikes.at[0, 'FWHM']
    lst_spike_FWHM = adaptation_spikes.tail(n_lastspikes)['FWHM'].mean()
    spike_FWHM_adaptation = lst_spike_FWHM / fst_spike_FWHM
    
    
    # # # spike frequency adaptation # # #
    
    # get last spikes
    n_last_ISIs = adaptation_n_lastspikes - 1
    
    # compare first and mean of last few ISIs
    fst_ISI = adaptation_ISIs.at[1, 'inst_freq']
    lst_ISIs = adaptation_ISIs.tail(n_lastspikes)['inst_freq'].mean()
    freq_adaptation_ratio = lst_ISIs / fst_ISI
    
    
    # # # spike frequency adaptation steady state # # #
    
    ## linear fit to ISIs between 500 and 1250 ms
    # limit dataframe to ISIs to ISIs later in step
    ISIs_for_linear_fit = adaptation_ISIs.query('t_ISI > 500')
    
    # linear fit
    lst_tISIs = ISIs_for_linear_fit['t_ISI']
    lst_inst_freqs = ISIs_for_linear_fit['inst_freq']
    
    if len(lst_inst_freqs) > 3:
        # fit linear curve
        popt_adapfreq, pcov = sc.optimize.curve_fit(linear_func, lst_tISIs, lst_inst_freqs, 
                                                    p0 = adaptation_popt_guess_linear_ISIs, 
                                                    maxfev = 5000)
        
        # plot exponential fit
        t_linfit = np.arange(start = ISIs_for_linear_fit['t_ISI'].iat[0], 
                              stop = ISIs_for_linear_fit['t_ISI'].iat[-1], 
                              step = 1)
        
        # get incline of linear fit
        freq_adaptation_incline_linearfit = popt_adapfreq[0] # Hz / ms
    
        # convert to Hz / s
        freq_adaptation_incline_linearfit = freq_adaptation_incline_linearfit * 1e3
    
    else:
        print(f'{cell_ID} linear fit not achieved. Number of spikes after 250 ms to low.')
        freq_adaptation_incline_linearfit = np.nan
        t_linfit = np.nan
        popt_adapfreq = np.nan
        
    ## save values to dataframe
    adaptation.at[cell_ID, 'spike_amplitude_adaptation_ratio'] = spike_amplitude_adaptation
    adaptation.at[cell_ID, 'spike_FWHM_adaptation_ratio'] = spike_FWHM_adaptation
    adaptation.at[cell_ID, 'freq_adaptation_ratio'] = freq_adaptation_ratio
    adaptation.at[cell_ID, 'freq_adaptation_steadystate'] = freq_adaptation_incline_linearfit
    
    if True:
        plot_adaptation(cell_ID, t_full, v_full, 
                        idx_rheo, idx_maxfreq, idx_halfmax, idx_maxinitinstfreq, 
                        IF, IF_inst_init, 
                        i_rheo_abs, i_maxfreq, i_halfmax, i_maxinitinstfreq,
                        n_lastspikes,
                        adaptation_spikes, 
                        rheobase_spikes, maxfreq_spikes, halfmax_spikes, maxinitinstfreq_spikes, 
                        adaptation_ISIs, 
                        rheobase_ISIs, maxfreq_ISIs, halfmax_ISIs, maxinitinstfreq_ISIs, 
                        fst_ISI, lst_ISIs, lst_inst_freqs, freq_adaptation_ratio, popt_adapfreq, t_linfit, freq_adaptation_incline_linearfit,
                        fst_spike_vamplitude, lst_spike_vamplitude, spike_amplitude_adaptation, 
                        fst_spike_FWHM, lst_spike_FWHM, spike_FWHM_adaptation)


# %% saving

print('saving...')

# tobe saved
export_vars = {'IF' : IF, 
               'IF_inst' : IF_inst, 
               'IF_inst_init' : IF_inst_init,
               'IF_dict' : IF_dict, 
               'IF_rheobase' : IF_rheobase, 
               'adaptation' : adaptation}

export_prefix = 'cc_IF-syn-'
export_extension = '.xlsx'


# def export_vars_to_celldescriptors(export_vars_dict, export_prefix, export_extension, index_label):
    
#     from parameters.directories_win import cell_descrip_syn_dir

for export_name, export_var in export_vars.items():
    
    rows = export_var.index.to_list()
    cols = export_var.columns.to_list()
    index_label = export_var.index.name
    
    # try loading and writing or create new file
    try:
        loaded_export_var = pd.read_excel(join(cell_descrip_syn_dir, export_prefix + export_name + export_extension),
                                          index_col = index_label)
        
        # find out how to combine both dataframes
        loaded_export_var.loc[rows, cols] = export_var.loc[rows, cols].values
        
        # save activity dataframe
        loaded_export_var.to_excel(join(cell_descrip_syn_dir, export_prefix + export_name + export_extension), 
                                   index_label=index_label)
    
    except FileNotFoundError:
        # save activity dataframe
        export_var.to_excel(join(cell_descrip_syn_dir, export_prefix + export_name + export_extension), 
                            index_label=index_label)


# export_vars_to_celldescriptors(export_vars_cols, export_prefix, export_extension, 'i_input')

# %% update analyzed cells

from functions.update_database import update_analyzed_sheet
    

update_analyzed_sheet(cell_IDs, PGF = PGF)
