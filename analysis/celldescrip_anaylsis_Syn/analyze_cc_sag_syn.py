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

# %% check for new cells to be analyzed


# %% initialize plotting and verificaiton plots

# init plotting
from functions.initialize_plotting import *

# verification plots
vplots = True
if vplots:
    # load plotting functions
    from analysis.celldescrip_anaylsis_Syn.plot_analyze_cc_sag_syn import plot_full_sag, plot_passive_props_calc

  
    
# %% load

# to check: E-231, 236, 248, 290
# cell_IDs = ['E-317']

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

    passive_props_calcs = pd.DataFrame(columns = ['r_input', 'tau_mem', 'c_mem', 'v_min', 'idx_vmin', 't_vmin', 'popt', 'r_squared', 'useable'],
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
        # U=R*I -> R = ΔU / ΔI
    
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
        passive_props_calcs.loc[step, :] = [r_input, tau_mem, c_mem, v_min, idx_vmin, t_vmin, popt, r_squared, useable_step]
        
    # get first 3 useful steps
    passive_props_calcs = passive_props_calcs.query('useable == 1').head(3)
    
    # raise warning if less than three steps are usable
    if passive_props_calcs.shape[0] < 3:
        warnings.warn(f'{cell_ID} used less than 3 hyperpolarising steps for passive properties!')

    # get passive properties as means of first 3 useable steps
    passive_properties['r_input'] = passive_props_calcs['r_input'].mean()
    passive_properties['tau_mem'] = passive_props_calcs['tau_mem'].mean()
    passive_properties['c_mem']   = passive_props_calcs['c_mem'].mean()
    
    if vplots:
        plot_passive_props_calc(cell_ID,
                                t_full,
                                v_full,
                                t_stim,
                                v_stim,
                                passive_props_calcs,
                                SR)
        
# %% 