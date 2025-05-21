# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 14:00:29 2025

@author: nesseler

!! must be run in virtual environment containing all necessary packages !!

0
install virutal environment package
> use 'pip install virtualenv'

1
creation of virtual environment
> create folder to host virtual environment
> use anaconda prompt to navigate to folder
> use 'python<version> -m venv <virtual-environment-name>' to create virtual env
    > 'python3.10 -m venv miniML_venv' OR 'conda create -n miniML_venv python=3.10'

2
de-/activate virtual environment
> To activate this environment, use
    - conda activate miniML_venv
> To deactivate an active environment, use
    - conda deactivate

3
install kernel for spyder in virtual environment
> conda install spyder-kernels=2.5 (or 3.0)

4 
install dependencies
- use anaconda prompt to navigate to miniML folder
- activate virtual environment
- run 'pip install -r requirements.txt' to install packages

5
install additional packages (using pip install in venv)
> tqdm
> openpyxl

"""

from scipy.ndimage import maximum_filter1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mtl
import numpy as np
import pandas as pd
import scipy as sc
import sys
from miniML import MiniTrace, EventDetection
from miniML_plot_functions import miniML_plots
from tqdm import tqdm

# init plotting
mtl.rcParams.update({'font.size': 9, 'font.family' : 'Arial'})

# experimental setting
hold = 'Erest'
treat = 'ctrl'

# parameter
scaling = 1e12
unit = 'pA'
SR = 100_000

# analysis parameter
filter_freq = 750
filter_order = 3
direction = 'negative'
factor = 19
model_th = 0.5

# output
plot = True
save = True

# paths
ePhys_parent = 'Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA'
rawData_path = 'Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/RAW_data/'
syn_path = 'Z:/n2021_MOS_AOS_Integration/ePhys-BAOT_MeA/synaptic_currents/'
data_path = syn_path + 'miniML_dtc-validation/'

# load xlsx sheet
lookup = pd.read_excel(ePhys_parent + '/ePhys-database.xlsx', 
                       sheet_name = 'Syn_testdataset', 
                       index_col = 'cell_ID')

# get cell_IDs
# cell_IDs = lookup.index.to_list()
# cell_IDs = cell_IDs[2:]
# cell_IDs.remove('E-303')
cell_IDs = ['E-298']

# define multiple factors
# factors = np.arange(2, 50+.1, 2, dtype = int)
# factors = np.append(factors, 19)
# factors = [int(f/6) for f in [36, 96, 114, 276]]
factors = [int(f/6) for f in [12, 24]]

ths = [0.5] #[0.5, 0.75, 0.9]

# newly analyzed
newly_analyzed = str()

# %%

for cell_ID in cell_IDs:
    
    # get info from ePhys xlsx sheet    
    rawData_filename = lookup.at[cell_ID, 'file'] + '.dat'
    rectype = lookup.at[cell_ID, 'Protocol']
    PGF_idx = lookup.at[cell_ID, 'PGF'] -1
    PGF_max = lookup.at[cell_ID, 'PGF_max'] 

    # load trace to MiniTrace
    trace = MiniTrace.from_heka_file(filename=(rawData_path + rawData_filename),
                                     rectype = rectype,
                                     group = lookup.at[cell_ID, 'group'],
                                     exclude_series = np.delete(np.arange(0, PGF_max), PGF_idx),
                                     scaling = scaling,
                                     unit = unit)
    
    # filter trace
    b, a = sc.signal.bessel(filter_order, filter_freq, fs = SR)
    data_filtered = sc.signal.lfilter(b, a, trace.data)

    # create new miniTrace object
    trace = MiniTrace(data = data_filtered,
                      sampling_interval = 1 / SR,
                      y_unit = unit,
                      filename = 'None')

    if plot:
        trace.plot_trace()
        
        
    for factor in factors:
        for model_th in ths:
            
            print(f'\n{cell_ID} - factor : {factor} - th : {model_th}\n')

            # # set settings for event detection
            eventdetection_settings = {'window_size' : 600 * factor,
                                       'model_threshold' : model_th,
                                       'batch_size' : 512,
                                       'event_detection_peakw' : 5,
                                       'stride' : 30,
                                       'rel_prom_cutoff' : 0.25,
                                       'convolve_win' : 20 * factor}
        
            # run prediction
            detection = EventDetection(data=trace,
                                       model_path='C:/Users/nesseler/miniML/models/GC_lstm_model.h5',
                                       window_size = eventdetection_settings['window_size'],
                                       model_threshold = eventdetection_settings['model_threshold'],
                                       batch_size = eventdetection_settings['batch_size'],
                                       event_direction=direction,
                                       compile_model=True,
                                       verbose=2)
        
            # detect events
            detection.detect_events(eval=True,
                                    stride = eventdetection_settings['stride'],
                                    peak_w = eventdetection_settings['event_detection_peakw'],
                                    rel_prom_cutoff = eventdetection_settings['rel_prom_cutoff'],
                                    convolve_win = eventdetection_settings['convolve_win'],
                                    resample_to_600 = True)
        
            # opt: miniML plot
            if plot:
                MiniPlots = miniML_plots(data=detection)
                MiniPlots.plot_prediction(include_data=True, 
                                          plot_filtered_prediction=True, 
                                          plot_filtered_trace=True,
                                          plot_event_params=True)
                MiniPlots.plot_event_overlay()
                MiniPlots.plot_event_histogram(plot='amplitude', 
                                               cumulative=False)
        
            if save:
                # calc window size and convert to string
                winsize_ms = (600 * factor) / (SR/1e3)
                winsize_str = str(int(winsize_ms))
                
                # get model threshold
                th_str = str(model_th).replace('.', 'p')
                
                # create filename
                filename = f'miniMLdetect_{cell_ID}_{hold}_{treat}_{winsize_str}_{th_str}' 
                
                if (factor in [int(f/6) for f in [36, 96, 114, 276]]):
                    incl_bool = True
                else:
                    incl_bool = False
                
                # save to pickle file
                detection.save_to_pickle(filename = data_path + filename + '.pickle', 
                                         include_prediction = incl_bool, 
                                         include_data = incl_bool)
            
                newly_analyzed = newly_analyzed + f' {cell_ID}_{hold}_{treat}_{winsize_str}_{th_str}' + '\n'


# %% write txt file

from datetime import datetime

if save:

    # get date and time
    datetime_str = datetime.now().strftime('%Y%m%d_%H%M')
    
    # set path and filename
    settingsoutput_path = syn_path + '/miniML_validation/'
    txt_filename = 'miniML_dtc-validation_settings-' + datetime_str
    
    with open(settingsoutput_path + txt_filename + '.txt', 'w+') as txt:
        
        # filter
        txt.write(f'Filter\n Filter type: bessel\n Filter freq: {filter_freq}\n Filter order: {filter_order}')
    
        # PGF
        txt.write(f'\n\nPGF\n eHold: {hold}\n Treatment: {treat}\n Scaling: {scaling}\n Unit: {unit}\n SR [Hz]: {SR}')
            
        # analysis dict
        analy_s = '\n\n'
        for k, v in eventdetection_settings.items():
            analy_s = analy_s + ' ' + k.capitalize() + ': ' + str(v) + '\n'
        txt.write(analy_s)
            
        # analyzed cells
        txt.write('\nAnalyzed:\n' + newly_analyzed)
     
        
# %% final

print('\nDone!')
