# -*- coding: utf-8 -*-
"""
Created on Tue May 13 17:45:35 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, synaptic_dir
from parameters.parameters import PSC_bins

# set cell_ID
cell_IDs = ['E-298', 'E-301', 'E-302', 'E-303', 'E-309', 'E-310', 'E-314']
cell_ID = cell_IDs[3]

winsize = 114
th = 0.5

# create dataframe
allevents  = pd.DataFrame(columns = ['cell_ID', 'amplitude', 'risetime', 'halfdecaytime', 'charge', 'score', 'event',
                                     'baseline', 'onset_idx', 'peak_idx', 'event_idx',
                                     'rise_min_idx', 'rise_min_value', 'rise_max_idx', 'rise_max_value'])


# %% map to dataframe function

def map_detection_to_df(detection, df, cell_ID):
    """
    This function maps the relevant values from the detection dict to a dataframe.

    Parameters
    ----------
    detection : dict
        DESCRIPTION.
    df : pandas.DataFrame
        DESCRIPTION.
    cell_ID : str
        DESCRIPTION.

    Returns
    -------
    pandas.DataFrame

    """
    
    n_events = detection['individual_values']['amplitudes'].shape[0]
    
    temp_df = pd.DataFrame.from_dict({'cell_ID' : [cell_ID] * n_events,
                                      'amplitude' : detection['individual_values']['amplitudes'],
                                      'risetime' : detection['individual_values']['risetimes'],
                                      'halfdecaytime' : detection['individual_values']['half_decaytimes'],
                                      'charge' : detection['individual_values']['charges'],
                                      'score' : detection['event_location_parameters']['event_scores'],
                                      'event' : detection['events'].tolist(),
                                      'baseline' : detection['event_location_parameters']['event_baselines'], 
                                      'onset_idx' : detection['event_location_parameters']['event_onset_locations'], 
                                      'peak_idx' : detection['event_location_parameters']['event_peak_locations'], 
                                      'event_idx' : detection['event_location_parameters']['event_locations'],
                                      'rise_min_idx' : detection['event_location_parameters']['min_positions_rise'], 
                                      'rise_min_value' : detection['event_location_parameters']['min_values_rise'], 
                                      'rise_max_idx' : detection['event_location_parameters']['max_positions_rise'], 
                                      'rise_max_value' : detection['event_location_parameters']['max_values_rise']})
    
    # append to full dataframe
    df = pd.concat([df, temp_df], ignore_index=True)

    return df


# %% load

for cell_ID in tqdm(cell_IDs):
    
    # set filename
    filename = f'miniMLdetect_{cell_ID}_Erest_ctrl_{winsize}_{str(th).replace(".", "p")}' 
    
    # open a file, where you stored the pickled data
    file = open((synaptic_dir + f'/miniML_dtc-validation/' + filename + '.pickle'), 'rb')
    
    # dump information to that file
    detection = pickle.load(file)
    
    # close and remove (from memory) the file
    file.close()
    del file 
    gc.collect()

    # write to dataframe
    allevents = map_detection_to_df(detection, df = allevents, cell_ID = cell_ID)


# %% calc histogram

score_bins = PSC_bins['scores']
scores_hist, _ = np.histogram(allevents['score'], bins = score_bins)


# %% figure

# init plotting
from functions.initialize_plotting import *

plt.stairs(values = scores_hist, 
           edges = score_bins,
           fill = False,
           lw = 0.75,
           color = colors_dict['primecolor'],
           label = 'score')





