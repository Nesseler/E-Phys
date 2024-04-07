# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:15:07 2024

@author: nesseler
"""


def clean_OnPath_column_to_path_ID_n_label(coordinates_dataframe):
    
    onPath_column = coordinates_dataframe['On Path']

    for i, txt in enumerate(onPath_column):
        path_ID = [int(s.replace("-", "").strip("()")) for s in txt.split() if s.replace("-", "").strip("()").isdigit()]
        
        if len(path_ID) > 1:
            raise ValueError('More than one integer in path_ID')
        elif len(path_ID) < 1:
            raise ValueError('Less than one integer in path_ID')
        elif len(path_ID) == 1:
            path_ID = path_ID[-1]
        
        
        # txt_ls = txt.split()
            
        if 'axon' in txt:
            path_label = 'axon'
        elif 'soma' in txt:
            path_label = 'soma'
        else:
            path_label = 'dendrite'
            
        if path_ID == 1:
            path_label = 'soma'
    
        coordinates_dataframe.at[i, 'path_ID'] = path_ID
        coordinates_dataframe.at[i, 'path_label'] = path_label
        
    coordinates_dataframe.drop(columns = ['On Path'], inplace = True)

    return coordinates_dataframe


    
    
    
    