# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:38:16 2025

@author: nesseler
"""

# import standard packages
from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_syn_dir

# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.functions_filter import merge_filter_split_steps
from functions.functions_useful import calc_dvdt_padded, exp_func, calc_rsquared_from_exp_fit
from functions.functions_extractspike import extract_spike, get_AP_parameters

# PGF specific
from parameters.PGFs import cc_sag_syn_parameters as PGF_parameters
t = PGF_parameters['t']
SR = PGF_parameters['SR']

# sag characteristation
from parameters.parameters import popt_guess, r_squared_thresh, v_expfit_thresh

# define protocol
PGF = 'cc_sag'
sheet_name = 'PGFs_Syn'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = sheet_name)

        
# %% define output

# passive properties
passive_properties = pd.DataFrame(columns=['r_input', 'tau_mem', 'c_mem'], 
                                  index = cell_IDs)
passive_properties.index.name = 'cell_IDs'

# sag properties
sag_properties = pd.DataFrame(columns=['sagdelta', 
                                       'n_reboundspikes', 
                                       'reboundspike_t_peak', 
                                       'reboundspike_t_threshold',
                                       'reboundspike_v_threshold',
                                       'reboundspike_v_amplitude',
                                       'reboundspike_t_toPeak',
                                       'reboundspike_t_rise',
                                       'reboundspike_FWHM',
                                       'reboundspike_AHP_amplitude'], 
                              index = cell_IDs)
sag_properties.index.name = 'cell_IDs'


# %% check for new cells to be analyzed


# %% initialize plotting and verificaiton plots

# init plotting
from functions.initialize_plotting import *

# verification plots
vplots = True
if vplots:
    # load plotting functions
    from analysis.celldescrip_anaylsis_Syn.plot_analyze_cc_sag_syn import plot_full_sag, plot_passive_props_calc, plot_sag_n_reboundspikes

  
    
# %% load

# to check: E-231, 236, 248, 290, 213, 258
# cell_IDs = ['E-214']

for cell_ID in tqdm(cell_IDs):

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
    
    if vplots:
        plot_full_sag(cell_ID, 
                      t = t_full, 
                      v = v_full, 
                      i = i_calc, 
                      i_input = i_input)
    
# %% passive properties

    passive_props_calcs_steps = pd.DataFrame(columns = ['r_input', 'tau_mem', 'c_mem', 'v_min', 'idx_vmin', 't_vmin', 'popt', 'r_squared', 'useable'],
                                             index = np.arange(1, n_steps, dtype = int))
    
    for step in np.arange(1, n_steps, dtype = int):
    
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
        
# %% sag potential 

    # sag deltas
    sagdeltas = pd.DataFrame(columns=['v_min', 't_vmin', 'sagdelta', 'v_steadystate'], 
                             index = np.arange(1, n_steps, dtype = int))
    sagdeltas.index.name = 'step_idx'

    # steady state
    from parameters.parameters import perc_step_steadystate
    
    # get indices for steady state
    idc_steadystate = idc_stim[-(int(len(idc_stim)*perc_step_steadystate)):]

    # calculate sag delta for each step
    for step in np.arange(1, n_steps, dtype = int):
        
        # get vmin
        v_min = passive_props_calcs_steps.at[step, 'v_min']
        tv_min = passive_props_calcs_steps.at[step, 't_vmin']
        
        # get steady-state voltage
        v_steady = v_full[step][idc_steadystate]
        
        # get mean of steady-state
        v_steadystate = np.mean(v_steady)
        
        # calc sag delta
        sagdelta = np.absolute(v_steadystate - v_min)
        
        # write to dataframe
        sagdeltas.loc[step, ['v_min', 't_vmin', 'sagdelta', 'v_steadystate']] = [v_min, t_vmin, sagdelta, v_steadystate]
    
    # get sagdelta for cell
    from parameters.parameters import min_potential_for_sag, cc_sag_holding_potential

    # get input resistance of cell
    r_input = passive_properties.at[cell_ID, 'r_input'] 
    
    # # calc current necessary to hyperpolarize cell to min_potential_for_sag (-130 mV)
    # # U = R * I -> ΔI = ΔU / R
    # delta_i_sag = (np.absolute(cc_sag_holding_potential - min_potential_for_sag) / r_input) *1e3 # pA
    
    # get first step where the local minimum of the step surpasses min_potential_for_sag
    sag_step = passive_props_calcs_steps[passive_props_calcs_steps['v_min'] < min_potential_for_sag].index.values[0]

    # get sagdelta
    sagdelta = sagdeltas.at[sag_step, 'sagdelta']
    
    # write to dataframe
    sag_properties.at[cell_ID, 'sagdelta'] = sagdelta

    
# %% reboundspike

    from parameters.parameters import min_peak_prominence_ccIF, min_peak_distance_ccIF, min_max_peak_width_ccIF
    from parameters.parameters import AP_parameters
    
    # sag_step = 5

    # get indices for post-stim period
    idc_post = np.arange((PGF_parameters['t_pre'] + PGF_parameters['t_stim']) * (SR/1e3),
                         (PGF_parameters['t_pre'] + PGF_parameters['t_stim'] + PGF_parameters['t_post']) * (SR/1e3),
                         dtype = int)

    # get voltage trace for post period
    vsag_post = v_full[sag_step][idc_post]
    t_post = t_full[idc_post]
    dvdt_post = calc_dvdt_padded(vsag_post, t_post)

    # detect reboundspike
    idc_spikes, dict_spikes = sc.signal.find_peaks(vsag_post, 
                                                   prominence = min_peak_prominence_ccIF, 
                                                   distance = min_peak_distance_ccIF * (SR/1e3),
                                                   width = np.multiply(min_max_peak_width_ccIF, (SR/1e3)))
    
    # calculate spike times in seconds
    t_spikes = np.divide(idc_spikes, (SR/1e3)) + PGF_parameters['t_pre'] + PGF_parameters['t_stim']

    
    # measure reboundspike if present
    if len(idc_spikes) > 0:
        
        # get reboundspike
        _, reboundspike_t, reboundspike_v, reboundspike_dvdt = extract_spike(t_post, vsag_post, dvdt_post, idc_spikes[0])
        
        # measure the rheobase spike
        reboundspike_params, _ = get_AP_parameters(t_post, vsag_post, dvdt_post, [idc_spikes[0]])
        
    else:
        # create the same empty dataframe
        reboundspike_params = pd.DataFrame(columns = AP_parameters, index = [0])
        reboundspike_t = np.nan 
        reboundspike_v = np.nan
        reboundspike_dvdt = np.nan
        
    # write to dataframe
    sag_properties.at[cell_ID, 'n_reboundspikes']            = len(idc_spikes)
    sag_properties.at[cell_ID, 'reboundspike_t_peak']        = reboundspike_params.at[0, 't_peaks']
    sag_properties.at[cell_ID, 'reboundspike_t_threshold']   = reboundspike_params.at[0, 't_threshold']
    sag_properties.at[cell_ID, 'reboundspike_v_threshold']   = reboundspike_params.at[0, 'v_threshold']
    sag_properties.at[cell_ID, 'reboundspike_v_amplitude']   = reboundspike_params.at[0, 'v_amplitude']
    sag_properties.at[cell_ID, 'reboundspike_t_toPeak']      = reboundspike_params.at[0, 't_toPeak']
    sag_properties.at[cell_ID, 'reboundspike_t_rise']        = reboundspike_params.at[0, 't_rise']
    sag_properties.at[cell_ID, 'reboundspike_FWHM']          = reboundspike_params.at[0, 'FWHM']
    sag_properties.at[cell_ID, 'reboundspike_AHP_amplitude'] = reboundspike_params.at[0, 'v_AHP_amplitude']
    
    # verification plot
    if vplots:
        plot_sag_n_reboundspikes(cell_ID, 
                                 t_full, 
                                 v_full, 
                                 sag_step, 
                                 sagdeltas, 
                                 t_spikes, 
                                 idc_post, 
                                 vsag_post, 
                                 dvdt_post,
                                 reboundspike_t, 
                                 reboundspike_v,
                                 reboundspike_dvdt)
                                        
    # save sagdeltas per cell
    from parameters.directories_win import quant_data_dir
    sagdeltas.to_excel(join(quant_data_dir, 'cc_sag-sagdeltas', f'{cell_ID}-sagdeltas.xlsx'),
                       index_label = 'step_idx')
    
# %% saving

print('\nsaving...')

# tobe saved
export_vars = {'passive_properties' : passive_properties, 
               'sag_properties' : sag_properties}

export_prefix = 'cc_sag-syn-'
export_extension = '.xlsx'


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
