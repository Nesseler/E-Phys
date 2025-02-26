# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 10:08:04 2025

@author: nesseler
"""

# %% cc_rest

try:
    print('running cc_rest analyis ...')
    from analysis.analysis_celldescrip_Syn.analyze_cc_rest_syn import cell_IDs
    
except ValueError as errormesseage:
    
    if str(errormesseage) == 'Nothing new to analyze!': 
        print('\tcc_rest: ' + str(errormesseage))
        
    else:
        raise RuntimeError('cc_rest not analyzed!')
        
else:
    print('\tcc_rest analyzed: ', cell_IDs)
    

# %% cc_IF

try:
    print('\nrunning cc_IF analyis ...')
    from analysis.analysis_celldescrip_Syn.analyze_cc_IF_syn import cell_IDs

except ValueError as errormesseage:

    if str(errormesseage) == 'Nothing new to analyze!':
        print('\tcc_IF: ' + str(errormesseage))
        
    else:
        raise RuntimeError('cc_IF not analyzed!')
        
else:
    print('\tcc_IF analyzed: ', cell_IDs)
    
    
# %% cc_sag

try:
    print('\nrunning cc_sag analyis ...')
    from analysis.analysis_celldescrip_Syn.analyze_cc_sag_syn import cell_IDs

except ValueError as errormesseage:
    
    if str(errormesseage) == 'Nothing new to analyze!':
        print('\tcc_sag: ' + str(errormesseage))
        
    else:
        raise RuntimeError('cc_sag not analyzed!')
        
else:
    print('\tcc_sag analyzed: ', cell_IDs)
    
    
# %% cc_th1AP

try:
    print('\nrunning cc_th1AP analyis ...')
    from analysis.analysis_celldescrip_Syn.analyze_cc_th1AP_syn import cell_IDs

except ValueError as errormesseage:
    
    if str(errormesseage) == 'Nothing new to analyze!':
        print('\tcc_th1AP: ' + str(errormesseage))
        
    else:
        raise RuntimeError('cc_th1AP not analyzed!')
        
else:
    print('\tcc_th1AP analyzed: ', cell_IDs)
    
    
# %% take out the trash

import gc
gc.collect() 

