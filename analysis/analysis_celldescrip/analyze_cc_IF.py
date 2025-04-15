# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:16:29 2024
Updated on Tue Apr 15 2025

@author: nesseler
"""

# import standard packages
from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_dir

# spike detection
from parameters.parameters import min_peak_prominence, min_peak_distance, min_max_peak_width

# IF characteristation
from parameters.parameters import t_init_inst_freq

# sag characteristation
from parameters.parameters import popt_guess, r_squared_thresh, v_expfit_thresh


# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.functions_filter import merge_filter_split_steps
from functions.functions_useful import calc_dvdt_padded, exp_func, calc_rsquared_from_exp_fit
from functions.functions_extractspike import get_AP_parameters, extract_spike, get_spiketrain_n_ISI_parameter

# PGF specific
from parameters.PGFs import cc_IF_parameters as PGF_parameters
t = PGF_parameters['t']
SR = PGF_parameters['SR']
    

# define protocol
PGF = 'cc_IF'
sheet_name = 'PGFs'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = sheet_name)


# %% define output

# passive properties
passive_properties = pd.DataFrame(columns=['r_input', 'tau_mem', 'c_mem'], 
                                  index = cell_IDs)

# IF relationship
IF           = pd.DataFrame(columns=cell_IDs, index = np.arange(-100, 1000, 1, int))
IF_inst      = pd.DataFrame(columns=cell_IDs, index = np.arange(-100, 1000, 1, int))
IF_inst_init = pd.DataFrame(columns=cell_IDs, index = np.arange(-100, 1000, 1, int))

# set index label
for df in [IF, IF_inst, IF_inst_init]:
    df.index.name = 'i_input'

IF_dict      = pd.DataFrame(columns = ['i_rheobase'       , 'idx_rheobase',
                                       'i_maxfreq'        , 'idx_maxfreq'        , 'maxfreq',      
                                       'i_halfmax'        , 'idx_halfmax',
                                       'i_maxinitinstfreq', 'idx_maxinitinstfreq', 'maxinitinstfreq'],
                            index = cell_IDs)

IF_rheobase = pd.DataFrame(columns = ['rheobase_abs', 'rheobase_rel'], 
                           index = cell_IDs)

adaptation = pd.DataFrame(index = cell_IDs)

# set index label
for df in [IF_dict, IF_rheobase, adaptation, passive_properties]:
    df.index.name = 'cell_ID'


# %% initialize plotting and verification plots

# init plotting
from functions.initialize_plotting import *

# verification plots
vplots = True
if vplots:
    # load plotting functions
    from analysis.analysis_celldescrip_Syn.plot_analyze_cc_IF_syn import plot_full_IF, plot_IF_step_spike_detection, plot_rheobase, plot_adaptation
    from analysis.analysis_celldescrip_Syn.plot_analyze_cc_sag_syn import plot_passive_props_calc, plot_sag_n_reboundspikes


# %% analysis

# cell_IDs = ['E-117']

for cell_ID in tqdm(cell_IDs):
    
    # %% load

    # load IF protocol
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name = sheet_name)
    
    # get data with file path & trace index
    i, v, _, _, n_steps = get_cc_data(file_path, traceIndex, scale='s')
    
    # filter steps
    v = merge_filter_split_steps(v, SR)
    
    # limit time and voltage to stimulation period
    idc_stim = PGF_parameters['idc_stim']
    
    v_full = np.copy(v)
    t_full = np.copy(t)
    v_stim = v_full[:, idc_stim]
    t_stim = t_full[idc_stim]
    
    # construct current array
    from functions.functions_constructors import construct_I_array
    i_calc, i_input = construct_I_array(cell_ID, n_steps,
                                        PGF = PGF,
                                        sheet_name = sheet_name,
                                        parameters = PGF_parameters)
    
    # get i_hold
    i_hold = i_calc[0][0]
    
    if vplots:
        plot_full_IF(cell_ID, 
                     t = t_full, 
                     v = v_full, 
                     i = i_calc, 
                     i_input = i_input)


    # %% passive properties

    # find step closest to I_hold (0 rel input)
    i_closest_2_hold = min(i_input, key=lambda x:abs(x-i_hold))
    
    # find index
    stepidx_hold = list(i_input).index(i_closest_2_hold)

    passive_props_calcs_steps = pd.DataFrame(columns = ['r_input', 'tau_mem', 'c_mem', 'v_min', 'idx_vmin', 't_vmin', 'popt', 'r_squared', 'useable'],
                                             index = np.arange(stepidx_hold, 0, -1, dtype = int))
    
    for step in np.arange(stepidx_hold-1, 0, -1, dtype = int):
    
        # define idc for pre & post time period
        idc_pre  = np.arange(50 * (SR/1e3), 200 * (SR/1e3), dtype = int)
        idc_post = np.arange(300 * (SR/1e3), 1200 * (SR/1e3), dtype = int)
    
        # get v and i pre
        v_pre = np.mean(v_full[step][idc_pre])
        i_pre = np.mean(i[step][idc_pre])
        
        # get i post
        i_post = np.mean(i[step][idc_post])
        
        # calc i_delta
        delta_i = np.absolute(i_post - i_pre)
           
        # find minimum and its index in step
        v_min = np.min(v_stim[step])
        idx_vmin = np.argmin(v_stim[step])
        t_vmin = (idx_vmin / (SR / 1e3)) + PGF_parameters['t_pre']
        
        # get delta v guess
        delta_vguess = np.absolute(v_min - v_pre)
        
        # limit voltage trace until minimum
        v_expfit = v_stim[step][:idx_vmin]
        
        # create x dimension for fit
        x_expFit = np.arange(0, len(v_expfit))
    
        # # fit exponential curve
        popt, pcov = sc.optimize.curve_fit(exp_func, x_expFit, v_expfit, p0 = [delta_vguess, *popt_guess[1:]], maxfev = 5000)
    
        # get r_squared as goodness of fit
        r_squared = calc_rsquared_from_exp_fit(x_expFit, v_expfit, popt)
        
        # get delta v from fit of exponential curve
        delta_v = popt[0]
        
        # calc voltage after hyperpolarisation from fit
        v_post_expfit = exp_func(30000, *popt)


        ### INPUT RESISTANCE    
        # U = R * I -> R = ΔU / ΔI
    
        # calculate r input in MOhm
        r_input = (delta_v / delta_i) * 1e3
        
        
        ### MEMBRANE TIME CONSTANT
        # tau_mem
        # time it takes the potential to reach 1 - (1/e) (~63%) of the max voltage
    
        # calc 1 - 1/e
        tau_perc_value = 1-(1/np.exp(1))
        
        # calc max voltage delta
        delta_v_63 = delta_v * tau_perc_value
        
        # calc 63 % of max voltage
        v_tau = v_pre + delta_v_63
        
        # calc time (indices first) it takes to reach v_tau
        idx_63 = - (np.log((v_tau - popt[2]) / (popt[0]))) / (popt[1])
        tau_mem = - idx_63 / (SR / 1e3) 
        
        
        ### MEMBRANE CAPACITANCE
        # c_mem = tau_mem / r_input in pF
        c_mem = (tau_mem / r_input) * 1e3
        
        
        # set step as useable for passive properties or not
        useable_step = 1
        
        if r_squared < r_squared_thresh:
            useable_step = 0
        elif v_post_expfit < v_expfit_thresh:
            useable_step = 0
        
        # write measurements to dataframe
        passive_props_calcs_steps.loc[step, :] = [r_input, tau_mem, c_mem, v_min, idx_vmin, t_vmin, popt, r_squared, useable_step]
        
    # get first 3 useful steps
    passive_props_calcs = passive_props_calcs_steps.query('useable == 1').head(3)
    
    # raise warning if less than three steps are usable
    if passive_props_calcs.shape[0] < 3:
        warnings.warn(f'{cell_ID} used less than 3 hyperpolarising steps for passive properties!')

    # get passive properties as means of first 3 useable steps
    passive_properties.at[cell_ID, 'r_input'] = passive_props_calcs['r_input'].mean()
    passive_properties.at[cell_ID, 'tau_mem'] = passive_props_calcs['tau_mem'].mean()
    passive_properties.at[cell_ID, 'c_mem']   = passive_props_calcs['c_mem'].mean()
    
    if vplots:
        plot_passive_props_calc(cell_ID,
                                t_full,
                                v_full,
                                t_stim,
                                v_stim,
                                passive_props_calcs,
                                SR)


    # %% active properties
    
    # add 100 ms after stimulation period to include APs at the end of the step
    t_additional = 50 #ms
    idc_additional = np.arange(idc_stim[-1] +1, idc_stim[-1] +1 + (t_additional * (SR/1e3)), dtype = int)
    idc_stim = np.append(idc_stim, idc_additional)
    
    # redefine traces in stimulation period
    v_stim = v_full[:, idc_stim]
    t_stim = t_full[idc_stim]
    
    
    # %% spike detection

    # iterate through steps
    for step in range(n_steps):
        
        # define voltage trace for step
        v_step = v_stim[step]
        
        # define current 
        i_step_rel = np.int64(i_input[step] - i_hold)
        
        # find peaks
        idc_spikes, dict_spikes = sc.signal.find_peaks(v_step, 
                                                       prominence = min_peak_prominence, 
                                                       distance = min_peak_distance * (SR/1e3),
                                                       width = np.multiply(min_max_peak_width, (SR/1e3)))
            
        # calculate spike times in seconds
        t_spikes = np.divide(idc_spikes, (SR/1e3))
        n_spikes = len(t_spikes)
        
        # calc freq as number of APs over 1000 ms
        freq = n_spikes / (PGF_parameters['t_stim'] / 1e3)
    
           
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
            if len(init_spikes) > 2:
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
        IF.at[i_step_rel, cell_ID] = freq
        IF_inst.at[i_step_rel, cell_ID] = inst_freq
        IF_inst_init.at[i_step_rel, cell_ID] = init_inst_freq
        
        if False:
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
    i_rheo_abs = np.int64(i_input[idx_rheo])
    i_rheo_rel = np.int64(i_input[idx_rheo] - i_hold)
    
    # find rheobase spike
    idc_spikes, dict_spikes = sc.signal.find_peaks(v_rheo, 
                                                   prominence = min_peak_prominence, 
                                                   distance = min_peak_distance * (SR/1e3),
                                                   width = np.multiply(min_max_peak_width, (SR/1e3)))
    
    # get rheospike
    _, rheospike_t, rheospike_v, rheospike_dvdt = extract_spike(t_full, v_rheo, dvdt_rheo, idc_spikes[0])
    
    # measure the rheobase spike
    rheospike_params, _ = get_AP_parameters(t_full, v_rheo, dvdt_rheo, [idc_spikes[0]])
    
    if vplots:
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
    IF_rheobase.loc[cell_ID, ['rheobase_abs', 'rheobase_rel']] = [i_rheo_abs, i_rheo_rel]
    
    # write rheobaseespike to dataframe
    IF_rheobase.loc[cell_ID, ['rheobasespike_vthreshold']]   = rheospike_params.loc[0, 'v_threshold']
    IF_rheobase.loc[cell_ID, ['rheobasespike_vamplitude']]   = rheospike_params.loc[0, 'v_amplitude']
    IF_rheobase.loc[cell_ID, ['rheobasespike_ttoPeak']]      = rheospike_params.loc[0, 't_toPeak']
    IF_rheobase.loc[cell_ID, ['rheobasespike_FWHM']]         = rheospike_params.loc[0, 'FWHM']
    IF_rheobase.loc[cell_ID, ['rheobasespike_AHPamplitude']] = rheospike_params.loc[0, 'v_AHP_amplitude']
    IF_rheobase.loc[cell_ID, ['rheobasespike_ttoAHP']]       = rheospike_params.loc[0, 't_toAHP']


    # %% adaptation
    
    # max freq
    maxfreq      = IF[cell_ID].dropna().max()
    idx_maxfreq  = IF[cell_ID].dropna().argmax()
    v_maxfreq    = v_stim[idx_maxfreq]
    dvdt_maxfreq = calc_dvdt_padded(v_maxfreq, t_stim)
    i_maxfreq    = i_input[idx_maxfreq]


    # max initial instant frequency    
    # edge-case: if initial instant freq has only one value
    if IF_inst_init[cell_ID].dropna().shape[0] == 1:
        i_maxinitinstfreq    = IF_inst_init[cell_ID].dropna().index.values[0]
        maxinitinstfreq      = IF_inst_init[cell_ID].dropna().values[0]
        idx_maxinitinstfreq  = np.where(IF[cell_ID].dropna().index == i_maxinitinstfreq)[0][0]
    
    # condition for E-223
    elif IF_inst_init[cell_ID].dropna().shape[0] == 0:
        i_maxinitinstfreq    = IF_inst[cell_ID].dropna().index.values[0]
        maxinitinstfreq      = IF_inst[cell_ID].dropna().values[0]
        idx_maxinitinstfreq  = np.where(IF[cell_ID].dropna().index == i_maxinitinstfreq)[0][0]
        
    else:
        i_maxinitinstfreq   = IF_inst_init[cell_ID].dropna().index[IF_inst_init[cell_ID].dropna().argmax()]
        maxinitinstfreq     = IF_inst_init[cell_ID].dropna().max()
        idx_maxinitinstfreq = np.where((IF[cell_ID].dropna().index == i_maxinitinstfreq) == True)[0][0]
        
    v_maxinitinstfreq    = v_stim[idx_maxinitinstfreq]
    dvdt_maxinitinstfreq = calc_dvdt_padded(v_maxinitinstfreq, t_stim)
        
    
    # step half-way between rheobase and max inst initial freq
    idx_halfmax  = int(((idx_maxinitinstfreq - idx_rheo) / 2) + idx_rheo)
    v_halfmax    = v_stim[idx_halfmax]
    dvdt_halfmax = calc_dvdt_padded(v_halfmax, t_stim)
    i_halfmax    = i_input[idx_halfmax]
    
    # write to dataframe
    IF_dict.at[cell_ID, 'i_rheobase']          = np.int64(i_rheo_abs - i_hold)
    IF_dict.at[cell_ID, 'idx_rheobase']        = np.int64(idx_rheo)
    IF_dict.at[cell_ID, 'i_maxfreq']           = np.int64(i_maxfreq - i_hold)
    IF_dict.at[cell_ID, 'idx_maxfreq']         = np.int64(idx_maxfreq)
    IF_dict.at[cell_ID, 'i_halfmax']           = np.int64(i_halfmax - i_hold)
    IF_dict.at[cell_ID, 'idx_halfmax']         = np.int64(idx_halfmax)
    IF_dict.at[cell_ID, 'i_maxinitinstfreq']   = np.int64(i_maxinitinstfreq - i_hold)
    IF_dict.at[cell_ID, 'idx_maxinitinstfreq'] = np.int64(idx_maxinitinstfreq)
    
    # write max firing frequencies to dataframe
    IF_dict.at[cell_ID, 'maxfreq']             = np.int64(maxfreq)
    IF_dict.at[cell_ID, 'maxinitinstfreq']     = np.float64(maxinitinstfreq)

        
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
        # print(f'\n\t{cell_ID} linear fit not achieved. Number of spikes after 250 ms to low.')
        freq_adaptation_incline_linearfit = np.nan
        t_linfit = np.nan
        popt_adapfreq = np.nan
        
    ## save values to dataframe
    adaptation.at[cell_ID, 'spike_amplitude_adaptation_ratio'] = spike_amplitude_adaptation
    adaptation.at[cell_ID, 'spike_FWHM_adaptation_ratio'] = spike_FWHM_adaptation
    adaptation.at[cell_ID, 'freq_adaptation_ratio'] = freq_adaptation_ratio
    adaptation.at[cell_ID, 'freq_adaptation_steadystate'] = freq_adaptation_incline_linearfit
    
    if vplots:
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

export_vars = {'passive_properties' : passive_properties,
               'IF' : IF, 
               'IF_inst' : IF_inst, 
               'IF_inst_init' : IF_inst_init,
               'IF_dict' : IF_dict, 
               'IF_rheobase' : IF_rheobase, 
               'adaptation' : adaptation}

# iterate through dataframes dict
for export_name, export_var in export_vars.items():
    
    # get index label
    index_label = export_var.index.name
        
    # export as excel file
    export_var.to_excel(join(cell_descrip_dir, 'cc_IF-' + export_name + '.xlsx'), 
                        index_label=index_label)


# %% update analyzed cells

from functions.update_database import update_analyzed_sheet
    
update_analyzed_sheet(cell_IDs, PGF = PGF)




