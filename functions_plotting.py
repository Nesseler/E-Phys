#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:59:49 2023

@author: moritznesseler
"""

import matplotlib as mtl
import matplotlib.pyplot as plt
import numpy as np


def get_colors(darkmode_bool=False):

    if darkmode_bool:
        plt.style.use('dark_background')
        primecolor = 'w'
        color1 = 'cyan'
        color2 = 'magenta'
        color3 = [0, 1, 0]
        cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","white","magenta"])
        #plt.rcParams['axes.grid'] = False
        # plt.grid(False)
        plot_dict = {'color':primecolor, 'linewidth' : 0.5}
        seccolor = 'k'
    elif darkmode_bool == False:
        plt.style.use('default')
        primecolor = 'k'
        color1 = 'blue'
        color2 = 'purple'
        color3 = 'red'
        cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","grey","red"])
        plt.rcParams['axes.grid'] = True
        plot_dict = {'color':primecolor, 'linewidth' : 0.5}
        seccolor = 'w'
        
    colors_dict = {'primecolor': primecolor,
                   'color1': color1,
                   'color2': color2,
                   'color3': color3,
                   'cmap': cmap,
                   'plot_dict': plot_dict,
                   'seccolor' : seccolor}
        
    return colors_dict


def get_figure_size(width = 328.67, height = 165.5):
    mm = 1/25.4
    figsize=(width*mm, height*mm)
    return figsize


def get_colorcode(x, y, data_fc, norm=None, cmap='seismic', plot_dict={'c':'k'}, return_bool = False):
    '''
    Function creates a specified colorcode form data and auxilliary data that 
    will be used to create the color-code.
    Parameters:
        x : data on x-axis
        y : data on y-axis
        data_fc : data used to create the colorcode
        norm : matplotlib normalisation object. Default is None (then norm object will be newly created).
        cmap : color-map. Default is 'seismic'
        plot_dict : plotting dictionary 
        return_bool : Boolean value wether norm object and line collection should be returned, otherwise
                      only the line collection will be returned. Default is False.
    Returns:
        (norm, lc) : Tuple of norm object and line collection.
        lc : Line collection.
    '''
    
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    if norm is None:
        # Create a continuous norm to map from data points to colors
        norm = plt.Normalize(data_fc.min(), data_fc.max())
        return_bool = True
        
        
    lc = mtl.collections.LineCollection(segments, cmap=cmap, norm=norm)
    # Set the values used for colormapping
    lc.set_array(data_fc)
    
    return (norm, lc) if return_bool else lc


def remove_x_ticks_between(axes, n_layers):
    for i in range(1,n_layers):
        axes[i].tick_params(axis = 'y', size = 0)
        
#saving the figure
def save_figures(figure, figure_name, save_dir, darkmode_bool):
    
    import os.path
    
    if darkmode_bool == True:
        figure_name += " dark"
    else:
        figure_name += " light"
    
    figure.savefig(os.path.join(save_dir, os.path.normpath(figure_name + ".png")), format = 'png')
    figure.savefig(os.path.join(save_dir, os.path.normpath(figure_name + ".svg")), format = 'svg')
    
    
def set_font_sizes(small_font_size = 14, large_font_size = 16):
    '''
    Function sets font sizes of select text elements in figure to provided sizes.
    Parameters:
        small_font_size : Small font size for regular text. Default is 14.
        large_font_size : Large font size for titles and headings. Default is 16.
    '''
    
    plt.rc('font', size = small_font_size)
    plt.rc('axes', titlesize = small_font_size, 
                   labelsize = small_font_size,
                   linewidth = 0.5)
    plt.rc('xtick', labelsize = small_font_size)
    plt.rc('ytick', labelsize = small_font_size)
    plt.rc('lines', linewidth = 2)


def return_segments(x, ys):
    '''
    Function returns segments from single x array and multiple y array to use
    with matplotlib.collections.LineCollection.
    Parameters:
        x : Single array of common x-coordinates.
        ys : 2D array / list of y-coordinates to be plotted on common x-coordinates.
    Returns:
        segs : Segmented lines ready to use with 'LineCollection(segs)'
    '''
    n_ys = len(ys)
    n_x = len(x)
    segs = np.zeros((n_ys, n_x, 2))
    segs[:, :, 1] = ys
    segs[:, :, 0] = x
    return segs

    
    
    
    
    
    