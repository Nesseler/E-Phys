# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:38:22 2025

@author: nesseler
"""

from pandas import DataFrame, concat

def map_detection_to_df(detection, df, cell_ID):
    """
    This function maps the relevant values from the detection dict to a dataframe.

    Parameters
    ----------
    detection : dict
        miniML detection object saved as pickle and read as dict.
    df : pandas.DataFrame
        DataFrame is will be appended & concatenated with the event values.
    cell_ID : str
        Unique cell identifier, like 'E-303'.

    Returns
    -------
    pandas.DataFrame

    """
    
    n_events = detection['individual_values']['amplitudes'].shape[0]
    
    temp_df = DataFrame.from_dict({'cell_ID' : [cell_ID] * n_events,
                                   'amplitude' : detection['individual_values']['amplitudes'],
                                   'risetime' : detection['individual_values']['risetimes'],
                                   'halfdecaytime' : detection['individual_values']['half_decaytimes'],
                                   'charge' : detection['individual_values']['charges'],
                                   'score' : detection['event_location_parameters']['event_scores'],
                                   'event' : detection['events'].tolist(),
                                   'baseline' : detection['event_location_parameters']['event_baselines'],
                                   'bsl_start' : detection['event_location_parameters']['event_baseline_starts'], 
                                   'bsl_end' : detection['event_location_parameters']['event_baseline_ends'], 
                                   'onset_idx' : detection['event_location_parameters']['event_onset_locations'], 
                                   'peak_idx' : detection['event_location_parameters']['event_peak_locations'], 
                                   'event_idx' : detection['event_location_parameters']['event_locations'],
                                   'halfdecay_idx' : detection['event_location_parameters']['half_decay_positions'],
                                   'rise_min_idx' : detection['event_location_parameters']['min_positions_rise'], 
                                   'rise_min_value' : detection['event_location_parameters']['min_values_rise'], 
                                   'rise_max_idx' : detection['event_location_parameters']['max_positions_rise'], 
                                   'rise_max_value' : detection['event_location_parameters']['max_values_rise']})
    
    # append to full dataframe
    df = concat([df, temp_df], ignore_index=True)

    return df