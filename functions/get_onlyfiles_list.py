# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:56:24 2024

@author: nesseler
"""

def get_onlyfiles_list(dir_path):

    from os import listdir
    from os.path import isfile, join  

    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]