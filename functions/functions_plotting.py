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
        
        BAOT_color = '#7a66fc' 
        MeA_color =  '#ff8d00'
        BAOT_MeA_color = 'gray'
        
        # BAOT_color = '#43388a'
        # MeA_color = '#ff7d00'
        # BAOT_MeA_color = 'gray'
        
        
    colors_dict = {'primecolor': primecolor,
                   'color1': color1,
                   'color2': color2,
                   'color3': color3,
                   'cmap': cmap,
                   'plot_dict': plot_dict,
                   'seccolor' : seccolor,
                   'BAOT_lighter' : '#cac2fe',
                   'MeA_lighter' : '#ffd199'}
    
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
        


def add_measures_to_dataframe(dataframe_to_save, axis_for_calcs):
    '''
    Function to add parameters to the end of a DataFrame.
    Parameters:
        dataframe_to_save (Pandas DataFrame): DataFrame that includes all data plotted in the figure.
        axis_for_calcs (int): Axis along calculations in dataframe are made. Default is 0.
    Returns:
        dataframe_to_save (Pandas DataFrame): Modified DataFrame with included measurements.
    '''
    dataframe_to_save.loc['mean'] = dataframe_to_save.mean(axis = axis_for_calcs, numeric_only = True)
    dataframe_to_save.loc['median'] = dataframe_to_save.median(axis = axis_for_calcs, numeric_only = True)
    dataframe_to_save.loc['std'] = dataframe_to_save.std(axis = axis_for_calcs, numeric_only = True)
    dataframe_to_save.loc['quantile_0p25'] = dataframe_to_save.quantile(q = 0.25, axis = axis_for_calcs, numeric_only = True)
    dataframe_to_save.loc['quantile_0p50'] = dataframe_to_save.quantile(q = 0.50, axis = axis_for_calcs, numeric_only = True)
    dataframe_to_save.loc['quantile_0p75'] = dataframe_to_save.quantile(q = 0.75, axis = axis_for_calcs, numeric_only = True)

    return dataframe_to_save



#saving the figure
def save_figures(figure, figure_name, save_dir, 
                 darkmode_bool = None, figure_format = 'png', saving_feedback = False,
                 dataframe_to_save = None, index_label = 'cell_ID', add_measures = True, axis_for_calcs = 0,
                 groups_bool = False, groups = ['BAOT/MeA', 'MeA', 'BAOT'], groups_name = 'Region'):
    '''
    Function to save figures.
    Parameters:
        figure (obj) : matplotlib figure object
        figure_name  (str): Name of figure as string. Will be used as filename.
        save_dir (str): Directory to save figure in. Str
        darkmode_bool (bool): Boolean to add descriptor in filename that specifies light- or darkmode. 
        Default is None type.
        figure_format (str): Choice of how to save figure. Default is png. (png or svg)
        dataframe_to_save (Pandas DataFrame): DataFrame that includes all data plotted in the figure. Default is None type.
        index_label (str): Label for index column in dataframe.
        add_measures (bool): Boolean to add or omit creation of the dataframe measurements. Default is True.
        axis_for_calcs (int): Axis along calculations in dataframe are made. Default is 0.
        groups_bool (bool): Boolean to save and calculate the DataFrame by different groups. Default is False.
        groups (list (of strings)): List of keys that can be provided to divided the dataframe into groups. Default is list of regions.
        groups_name (str): String of column label that should be used for division of dataframe into groups. Dafault is 'Region' key.
        
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
        
        

    if dataframe_to_save is not None:
        # check for strings in groups list that are not compatible with paths
        for i, group in enumerate(groups):
            if '/' in group:
                groups[i] = groups[i].replace('/', '_')
    
        # start saving the dataframe
        if not groups_bool:
            if add_measures:
                # add mean, median, stdev, quartiles
                export_dataframe = add_measures_to_dataframe(dataframe_to_save, axis_for_calcs)
            elif not add_measures:
                export_dataframe = dataframe_to_save
            
            # save dataframe
            export_dataframe.to_excel(join(save_dir, figure_name + '.xlsx'), index_label = index_label)  
            
            if saving_feedback:
                print(f'"{figure_name}" dataframe saved at {save_dir}.')
            
        elif groups_bool:  
            # loop through groups
            for group in groups:
                # limit dataframe
                group_df = dataframe_to_save[dataframe_to_save[groups_name] == group]
                    
                if add_measures:
                    # add mean, median, stdev, quartiles
                    export_dataframe = add_measures_to_dataframe(group_df, axis_for_calcs)
                elif not add_measures:
                    export_dataframe = dataframe_to_save
                
                # save dataframe per group
                export_dataframe.to_excel(join(save_dir, figure_name + f'-{group}.xlsx'), index_label = index_label)    

                if saving_feedback:
                    print(f'"{figure_name}-{group}" dataframe saved at {save_dir}.')
            
    
    

     
    
    
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



# define function to change projection type of subplot specific subplot
def change_projection(fig, axs, ax_tochange, projection = 'polar'):

    rows, cols, start, stop = ax_tochange.get_subplotspec().get_geometry()

    if type(axs) == dict:
        axs[str(start)].remove()
        axs[str(start)] = fig.add_subplot(rows, cols, start+1, projection=projection) 
    else:
        axs.flat[start].remove()
        axs.flat[start] = fig.add_subplot(rows, cols, start+1, projection=projection)
    

    
    
    
    
    