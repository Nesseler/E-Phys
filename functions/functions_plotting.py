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
        plt.style.use('default')
        plt.style.use('dark_background')
        primecolor = 'w'
        color1 = 'cyan'
        color2 = 'magenta'
        color3 = 'red'
        cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","white","magenta"])
        plt.rcParams['axes.grid'] = False
        # plt.grid(False)
        plot_dict = {'color':primecolor, 'linewidth' : 0.5}
        seccolor = 'k'

        BAOT_color = '#7a66fc' #'#ff1b6b' #'#03C03C'
        MeA_color =  '#ff8d00' #'#45caff' #'#FF8F00'
        BAOT_MeA_color = 'gray'
        
    elif darkmode_bool == False:
        plt.style.use('default')
        primecolor = 'k'
        color1 = 'blue'
        color2 = 'purple'
        color3 = 'red'
        cmap = mtl.colors.LinearSegmentedColormap.from_list("", ["blue","grey","red"])
        plt.rcParams['axes.grid'] = False
        plot_dict = {'color':primecolor, 'linewidth' : 0.5}
        seccolor = 'w'
        
        BAOT_color = '#43388a'
        MeA_color = '#ff7d00'
        BAOT_MeA_color = 'gray'
        
        
    colors_dict = {'primecolor': primecolor,
                   'color1': color1,
                   'color2': color2,
                   'color3': color3,
                   'cmap': cmap,
                   'plot_dict': plot_dict,
                   'seccolor' : seccolor}
    
    regions_c = {'BAOT' : BAOT_color,
                 'MeA' : MeA_color,
                 'BAOT/MeA' : BAOT_MeA_color}
        
    return colors_dict, regions_c


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
def save_figures(figure, figure_name, save_dir, 
                 darkmode_bool = None, figure_format = 'png', saving_feedback = False):
    '''
    Function to save figures.
    Parameters:
        figure (obj) : matplotlib figure object
        figure_name  (str): Name of figure as string. Will be used as filename.
        save_dir (str): Directory to save figure in. Str
        darkmode_bool (bool): Boolean to add descriptor in filename that specifies light- or darkmode. 
        Default is None type.
        figure_format (str): Choice of how to save figure. Default is png. (png or svg)
    '''
    
    # lazy load join and normpath
    from os.path import join, normpath
    
    # add descriptor to figure name if light or darkmode is specified.
    if darkmode_bool == True:
        figure_name += " dark"
    elif darkmode_bool == False:
        figure_name += " light"
        
        
    # saving figure in different formats
    if figure_format == 'png':
        figure_name = figure_name + ".png"
        
        figure.savefig(join(save_dir, normpath(figure_name)), format = 'png')
        
        if saving_feedback:
            print(f'"{figure_name}" saved at {save_dir}.')
        
    elif figure_format == 'svg':
        figure_name = figure_name + ".svg"
        
        figure.savefig(join(save_dir, normpath(figure_name)), format = 'svg')
        
        if saving_feedback:
            print(f'"{figure_name}" saved at {save_dir}.')
            
    elif figure_format == 'both':
        figure.savefig(join(save_dir, normpath(figure_name + ".png")), format = 'png')
        figure.savefig(join(save_dir, normpath(figure_name + ".svg")), format = 'svg')
     
        if saving_feedback:
            print(f'"{figure_name}" saved at {save_dir}.')
            
    else:
        raise Warning(f'"{figure_name}" not saved. Figure format not specified correctly!')

     
    
    
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




    
    
    
    
    