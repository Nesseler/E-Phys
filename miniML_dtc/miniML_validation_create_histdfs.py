# -*- coding: utf-8 -*-
"""
Created on Wed May 21 10:58:36 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, synaptic_dir
from parameters.parameters import PSC_bins

# import functions
from functions.functions_pscs import map_detection_to_df

# set validation cell_IDs
cell_IDs = ['E-298', 'E-301', 'E-302', 'E-303', 'E-309', 'E-310', 'E-314']

th = 0.5
winsizes = [36, 96, 114, 276]


# %% set output dataframe dict

hist_dfs = {'score' : pd.DataFrame(index = PSC_bins['score'][:-1], columns = winsizes),
            'amplitude' : pd.DataFrame(index = PSC_bins['amplitude'][:-1], columns = winsizes),
            'halfdecaytime' : pd.DataFrame(index = PSC_bins['halfdecaytime'][:-1], columns = winsizes),
            'risetime' : pd.DataFrame(index = PSC_bins['risetime'][:-1], columns = winsizes)}


# %% load
 
# create dataframe
allevents = pd.DataFrame()

for winsize in tqdm(winsizes):
    for cell_ID in tqdm(cell_IDs, leave=False):
        
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
        
        # get previouse length
        pre_len = allevents.shape[0]
    
        # write to dataframe
        allevents = map_detection_to_df(detection, df = allevents, cell_ID = cell_ID, include_events=False)
    
        # get number of events
        n_events = detection['individual_values']['amplitudes'].shape[0]
        
        # get current length of allevents
        allevents_len = allevents.shape[0]
        
        # indexing
        allevents_idx = np.arange(pre_len, allevents_len, dtype = int)
    
        # extend with windows size
        allevents.loc[allevents_idx, 'winsize'] = [winsize] * int(allevents_idx.shape[0])
        
        # update previouse length
        pre_len = allevents_len

# drop nan valued rows
allevents.dropna(axis = 'index', how = 'any', inplace = True)
    
    
# %% calc histogram

for measure in hist_dfs.keys():
    for winsize in winsizes:
        
        # get occurances and write to df
        hist_dfs[measure].loc[:, winsize], _ = np.histogram(allevents[allevents['winsize'] == winsize][measure], bins = PSC_bins[measure])
            

# %% save dataframes

for measure, hist_df in hist_dfs.items():
    # save dataframes
    hist_df.to_excel(synaptic_dir + '/miniML_validation' + f'/{measure}_hist_{str(th).replace(".", "p")}.xlsx', 
                     index_label = f'{measure}_bins')
    
