#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:59:49 2023

@author: moritznesseler
"""

import matplotlib as mtl
import matplotlib.pyplot as plt


def get_colors(darkmode_bool=False):

    if darkmode_bool:
        plt.style.use('dark_background')
        primecolor = 'w'
        color1 = 'cyan'
        color2 = 'magenta'
        color3 = [0, 1, 0]
        cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","white","magenta"])
        #plt.rcParams['axes.grid'] = False
        plt.grid(False)
        plot_dict = {'color':primecolor, 'linewidth' : 0.5}
    elif darkmode_bool == False:
        plt.style.use('default')
        primecolor = 'k'
        color1 = 'blue'
        color2 = 'purple'
        color3 = 'red'
        cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","grey","red"])
        plt.rcParams['axes.grid'] = True
        plot_dict = {'color':primecolor, 'linewidth' : 0.5}
        
    colors_dict = {'primecolor': primecolor,
                   'color1': color1,
                   'color2': color2,
                   'color3': color3,
                   'cmap': cmap,
                   'plot_dict': plot_dict
                   }
        
    return colors_dict

def get_figure_size():
    mm = 1/25.4
    figsize=(328.67*mm, 165.5*mm)
    return figsize


def remove_x_ticks_between(axes, n_layers):
    for i in range(1,n_layers):
        axes[i].tick_params(axis = 'y', size = 0)