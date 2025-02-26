# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:07:52 2024

@author: nesseler
"""

# import standard packages
from functions.initialize_packages import *

# define protocol to analyze
PGF = 'vc_TTX_washin'

# get cell IDs
from functions.get_cell_IDs import get_cell_IDs_one_protocol
cell_IDs = get_cell_IDs_one_protocol(PGF, sheet_name = 'PGFs_Syn')
cell_IDs_leak = get_cell_IDs_one_protocol(PGF + '_leak', sheet_name = 'PGFs_Syn')

# load parameters
from parameters.PGFs import vc_TTX_washin_parameters, vc_TTX_washin_leak_parameters

# initialize plotting
from functions.initialize_plotting import *

# verification plots
vplots = False


# %% set peak detection parameters

peak_prominence = 40 #pA
peak_width = 50 # points 


# %% load and process

from functions.functions_import import get_traceIndex_n_file, get_vc_data, get_vc_leak_data
from functions.functions_useful import calc_time_series
from parameters.directories_win import table_file, vplot_dir
import re

# init dataframe for measurements
peakCurrents = pd.DataFrame(columns = cell_IDs + cell_IDs_leak)

# cell_ID = 'E-268' #E-268

# load lookup table
lookup_table = pd.read_excel(table_file,
                             sheet_name = 'PGFs_Syn',
                             index_col = 'cell_ID')


for cell_ID in tqdm(cell_IDs + cell_IDs_leak):
           
    # check if cell has leak subtraction
    if cell_ID in cell_IDs_leak:
        
        PGF = 'vc_TTX_washin_leak'
        PGF_params = vc_TTX_washin_leak_parameters
        
        # get the traceIndex and the file path string for data import functions
        traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name='PGFs_Syn')
        
        # load data
        i, v, ileak, t, SR, n_steps = get_vc_leak_data(file_path, traceIndex, scale = 'ms')
        
    else: 
        
        PGF = 'vc_TTX_washin'
        PGF_params = vc_TTX_washin_parameters
        
        # get the traceIndex and the file path string for data import functions
        traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name='PGFs_Syn')
        
        # load data
        i, v, t, SR, n_steps = get_vc_data(file_path, traceIndex, scale = 'ms')
        


    # convert sampling rate
    SR_ms = SR/1000
    
    # load which steps to use
    ## get steps from lookup table
    steps_str = lookup_table.at[cell_ID, PGF+'_steps']

    # find integers in string
    steps_start_stop = [int(i) - 1 for i in steps_str.split('-') if i.isdigit()]
    
    # convert to indices
    steps_start_stop[0] = steps_start_stop[0]-1
    
    # create iteratable list
    steps_to_use = np.arange(steps_start_stop[0], steps_start_stop[1])
    
    # redefine number of steps
    n_steps = len(steps_to_use)
    
    # iterate through steps
    for step_idx in range(n_steps):
        
        # limit data to stimulation period
        step = i[step_idx][int(PGF_params['t_pre'] * SR_ms):int(PGF_params['t_post'] * SR_ms + PGF_params['t_stim'] * SR_ms)]
    
        # find peaks
        peaks, properties = sc.signal.find_peaks(-step, prominence = peak_prominence, width = peak_width)
    
        # get peak measurements
        peak_idx = peaks[0]
        peak_current = step[peak_idx]
        t_peak = peak_idx / SR_ms
        
        # get measurements of peak base
        leftbase_idx = properties['left_bases'][0]
        leftbase = step[leftbase_idx]
        t_leftbase = leftbase_idx / SR_ms
        
        # calculate peak current data
        delta = leftbase - peak_current
        
        # write to dataframe
        peakCurrents.at[step_idx, cell_ID] = delta
    
        # %
        # create verification plot
        if vplots:
        
            vfig, vax = plt.subplots(nrows = 1,
                                     ncols = 1,
                                     figsize = get_figure_size(width = 150, height = 120),
                                     dpi = 100,
                                     layout = 'constrained')
            
            vfig.suptitle(f'{cell_ID} step: {step_idx}')
            
            # get trace
            current_trace = i[step_idx][int((PGF_params['t_pre']-5) * SR_ms):int((PGF_params['t_post']+5) * SR_ms + PGF_params['t_stim'] * SR_ms)]
        
            # calc time series
            t_plot = np.arange(-5, 25, step = 1/SR_ms)
        
            vax.plot(t_plot,
                     current_trace,
                     c = colors_dict['primecolor'],
                     lw = 1)
            
            # plot peak
            vax.scatter(t_peak, peak_current,
                        marker = 'x',
                        c = 'r',
                        s = 20,
                        zorder = 2)
            
            # plot leftbase
            vax.scatter(t_leftbase, leftbase,
                        marker = 'x',
                        c = 'g',
                        s = 20,
                        zorder = 2)
            
            # y axis
            apply_axis_settings(ax = vax, 
                                axis = 'y', 
                                ax_min = -6000, 
                                ax_max = 2000, 
                                pad = 50, 
                                step = 2000, 
                                stepminor = 250, 
                                label = 'Current [pA]')

            # x axis
            apply_axis_settings(ax = vax, 
                                axis = 'x', 
                                ax_min = -5, 
                                ax_max = 25, 
                                pad = 1, 
                                step = 5, 
                                stepminor = 1, 
                                label = 'Time [ms]' )
            
            # despine
            [vax.spines[spine].set_visible(False) for spine in ['top', 'right']]
            
            # display plot
            plt.show()
            
            # set directory for figure
            vfig_dir = join(vplot_dir, 'TTX_washin')

            # save figure
            save_figures(vfig, f'figure-TTXwashin-{cell_ID}-step{step_idx}', 
                          save_dir = vfig_dir,
                          darkmode_bool = darkmode_bool,
                          figure_format = 'both')

            
            

    
# %%

from parameters.directories_win import quant_data_dir, figure_dir

peakCurrents.to_excel(join(quant_data_dir, 'TTX_washin', 'TTX_washin-peakCurrents.xlsx'), index_label = 'step_idc')



