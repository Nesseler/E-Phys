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


def get_blkr_colors(darkmode_bool=False):
    
    if darkmode_bool:
        blkr_colors = {'GBZ'          : '#7b7fae',
                       'AP5_NBQX'     : '#DB675D'}
    else:
        blkr_colors = {'GBZ'          : '#7b7fae',
                       'AP5_NBQX'     : '#DB675D'}
        
    return blkr_colors



def get_figure_size(width = 318.67, height = 160.5):
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
            
    
    

     
    
    
# def set_font_sizes(small_font_size = 14, large_font_size = 16):
#     '''
#     Function sets font sizes of select text elements in figure to provided sizes.
#     Parameters:
#         small_font_size : Small font size for regular text. Default is 14.
#         large_font_size : Large font size for titles and headings. Default is 16.
#     '''
    
#     plt.rc('font', size = small_font_size)
#     plt.rc('axes', titlesize = small_font_size, 
#                    labelsize = small_font_size,
#                    linewidth = 0.5)
#     plt.rc('xtick', labelsize = small_font_size)
#     plt.rc('ytick', labelsize = small_font_size)
#     plt.rc('lines', linewidth = 2)



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
    

    
    
def plot_half_violin(data, ax,
                      data_pad_factor = 1,
                      v_resolution = np.nan,
                      v_kde_cutoff = 0.01,
                      v_abs_cutoff = [np.nan, np.nan],
                      v_bandwidth = np.nan,
                      v_position = 0,
                      v_direction = -1,
                      v_offset = -0.05,
                      v_width = 0.25,
                      v_color = 'w',
                      v_lw = 1,
                      v_baseline = False,
                      v_fill = False,
                      v_fillcolor = 'w',
                      v_filllw = 1,
                      v_zorder = 0):
    
    """
    Function plots half violin at specified position.
    Parameters:
        data : Array of data that is used for plotting.
        ax : matplotlib axis to plot violin on.
        data_pad_factor : Float, default is 1. Factor that data_range is multiplied
                          with for estimation of kernel density.
        
        v_resolution : Float (data scale), default is np.nan. Resolution of violin in 
                        data native scale.
        v_kde_cutoff : Float (kde percentage), default is 0.01. Cutoff that is used
                        to limit violin.
        v_abs_cutoff : List of floats (data scale), default is [np.nan, np.nan]. 
                        Absolute values that can be used to limit violin. Shape of 
                        list is [Min, Max]. When only one value is specified, the 
                        other must be set to np.nan and here the v_kde_cutoff is used.
                        If both values are specified the kde_cutoff must be np.nan to 
                        take effect.
        v_bandwidth : Float (0-1), default is np.nan. Sets bandwidth of density.
        v_position : Float, default is 0. Position of violin on x axis.
        v_direction : -1 or 1, default is -1. Direction factor of violin on x
                      axis.
        v_offset : Float, default is -0.05. Offset of violin from v_position.
        v_width : Float, default is 0.25. Width of violin on x axis
        v_color : Str, default is 'w'. Color of violin.
        v_lw : Float, default is 1. Linewidth of violin.
        v_baseline : Boolean, default is False. Boolean to plot the baseline of
                      the violin.
        v_fill : Boolean, default is False. Boolean to plot the fill of the violin.
        v_fillcolor : Str, default is 'w'. Color of fill.
        v_filllw : Float, default is 1. Linewidth of fill.
        v_zorder : Integer, default is 0. Z-Order of violin.
    
    Returns:
        half_violin : List of Line and Fillobjects. Order is [violin, fill] 
    """
    
    from scipy.stats import gaussian_kde
    
    # # function inputs
    # data = [5,60,85,100,140,-10,15,20,35,25,15,165,85,75,80,165,15,15,100,25,45,5]
    # ax = plt.gca()
    
    # # defaults
    # data_pad_factor = 1
    # v_resolution = 0.1 # in native data scale
    # v_kde_cutoff = 0.01 # in % default is 0.01
    # v_abs_cutoff = [-50, np.nan] # in native data scale, list of [min, max] for violin cutoff
    # v_bandwidth = np.nan
    
    # v_position = 0
    # v_direction = -1
    # v_offset = -0.05
    # v_width = 1
    
    # v_color = 'w'
    # v_lw = 1
    # v_baseline = True
    # v_fill = False
    # v_fillcolor = v_color
    # v_filllw = v_lw
    # v_zorder = 0
    
    # get distribution metrics
    data_n = len(data)
    data_min = np.min(data)
    data_max = np.max(data)
    data_range = data_max - data_min
    
    # pad data range for violin that extends over range
    padded_min = data_min - (data_range * data_pad_factor)
    padded_max = data_max + (data_range * data_pad_factor)
    
    # get kernel density estimate
    density = gaussian_kde(data)
    
    # scale factor (bandwidth)
    # bandwidth is calculated by n**(-1./(d+4)), with n datapoints and d dimensions
    # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde
    if not np.isnan(v_bandwidth):
        density.set_bandwidth(v_bandwidth)
        
    # calc resolution of violin with 1000 points between min and max
    if np.isnan(v_resolution):
        v_resolution = data_range / 1000

    # get range of value for density estimation
    vs = np.arange(padded_min, padded_max + v_resolution, v_resolution)
    
    # get kde (kernel density esimate)
    kde = density(vs)
    
    # normalize kde to maximum
    kde_normed = [(k / np.max(kde)) for k in kde]
    
    # # # violin cutoff # # #
    
    # condition for set kde cutoff
    if (not np.isnan(v_kde_cutoff) and np.isnan(v_abs_cutoff).all()):
        # remove values below 1 % on edges
        kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff or v < data_max and v > data_min)]
        # print('0')
        
    # condition for set min and max absolute cutoff
    elif (np.isnan(v_kde_cutoff) and not np.isnan(v_abs_cutoff).all()):
        # remove values outside absolute cutoff
        kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (v > v_abs_cutoff[0] and v < v_abs_cutoff[1])]
        # print('1')
      
    # condition for set kde_cutoff and either abs min or max cutoff
    elif (not np.isnan(v_kde_cutoff) and np.isnan(v_abs_cutoff).any()):
        # condition for set min
        if not np.isnan(v_abs_cutoff[0]):
            # remove values below 1 % on edges
            kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff and v > v_abs_cutoff[0])]
            # print('2.1')
            
        # condition for set max 
        elif not np.isnan(v_abs_cutoff[1]):
            # remove values below 1 % on edges
            kv_s = [[k, v] for k, v in zip(kde_normed, vs) if (k > v_kde_cutoff and v < v_abs_cutoff[1])]
            # print('2.2')
        
    # redefine violin
    kde_clipped = [kv[0] for kv in kv_s]
    vs_clipped = [kv[1] for kv in kv_s]
    
    # add direction, scale, offset ,and x position
    kde_withOffset = [v_direction * k * v_width + v_offset + v_position for k in kde_clipped]
    
    # set coordinates for baseline
    vs_baseline = vs_clipped
    ks_baseline = [v_offset + v_position] * len(vs_clipped)
    
    if not v_baseline:
        ks = kde_withOffset  # x
        vs = vs_clipped      # y
    
    
    elif v_baseline:
        # concatenate two x and y arrays
        halflen_baseline = int(len(vs_baseline) / 2)
        
        vs = vs_baseline[0:halflen_baseline][::-1] + vs_clipped + vs_baseline[halflen_baseline+1:-1][::-1]
        ks = ks_baseline[0:halflen_baseline][::-1] + kde_withOffset + ks_baseline[halflen_baseline+1:-1][::-1]
        
    
    # plot half violin
    half_violin_list = ax.plot(ks, vs,
                               color = v_color,
                               linewidth = v_lw,
                               label = '_nolegend_',
                               zorder = 2)
    
    # violin return objects
    half_violin = half_violin_list
    
    
    
    if v_fill:
        v_face = ax.fill_betweenx(y = vs_clipped,
                                  x1 = kde_withOffset,
                                  x2 = ks_baseline,
                                  color = v_fillcolor,
                                  linewidth = v_filllw,
                                  zorder = v_zorder)
        
        # append to return object
        half_violin.append(v_face)
    

    return half_violin


# # # functions to change axes # # #

def apply_axis_settings(ax, axis = 'y', 
                        ax_min = 0, 
                        ax_max = 100, 
                        pad = 1, 
                        step = 10, 
                        stepminor = 10, 
                        label = 'Label [unit]',
                        start_at_0 = False,
                        limits_n_0 = False,
                        ticklabels = None,
                        rotation = None,
                        pad_factor = 0.01):
    """
    Function uses specified settings to change yaxis layout of specified subplot.
    Parameter:
        ax : matplotlib axis
        axis : str, default is 'y'. Specifies axis that is edited.
        min : float, default is 0
        max : float, default is 100
        pad : float, default is 1, if set to None pad will be calculated as 1 %
              of range between min and max
        step : float, default is 10,
        stepminor : float, default is 10,
        label : string, default is 'label [unit]'
        ticklabels : list of str, default is None, can be set to specify the
                     ticklabels
        rotation : float, default is None, can be set to rotate the ticklabels
    """
    
    if start_at_0:
        ticks = np.arange(0, ax_max+ stepminor/3, step)
    elif limits_n_0:
        ticks = [ax_min, 0, ax_max]
    else:
        ticks = np.arange(ax_min, ax_max+ stepminor/3, step)
        
    ticksminor = np.arange(ax_min, ax_max+ stepminor/3, stepminor)
    
    # set padding as 1 % 
    if not pad:
        pad = abs(ax_max - ax_min) * pad_factor
        
    
    if axis == 'y':
        ax.set_ylim([ax_min - pad, ax_max + pad])    
        ax.set_yticks(ticks = ticks)
        ax.set_yticks(ticks = ticksminor, minor = True)
        ax.spines['left'].set_bounds([ax_min, ax_max])
        
        if label:
            ax.set_ylabel(label)
        
        if (ticklabels and not rotation):
            ax.set_yticklabels(labels = ticklabels)
                               
        elif (ticklabels and rotation):
            ax.set_yticklabels(labels = ticklabels, rotation = rotation)
        
    elif axis == 'x':
        ax.set_xlim([ax_min - pad, ax_max + pad])
        ax.set_xticks(ticks = ticks)
        ax.set_xticks(ticks = ticksminor, minor = True)
        ax.spines['bottom'].set_bounds([ax_min, ax_max])
        
        if label:
            ax.set_xlabel(label)
        
        if (ticklabels and not rotation):
            ax.set_xticklabels(labels = ticklabels)
                               
        elif (ticklabels and rotation):
            ax.set_xticklabels(labels = ticklabels, rotation = rotation)
    
    elif axis == 'z':
        ax.set_zlim([ax_min - pad, ax_max + pad])
        ax.set_zticks(ticks = ticks)
        ax.set_zticks(ticks = ticksminor, minor = True)
        ax.set_zticks(ticks = np.arange(ax_min, ax_max+ stepminor, stepminor), minor = True)
        # ax.spines['bottom'].set_bounds([ax_min, ax_max])
        
        if label:
            ax.set_zlabel(label)
        
        if (ticklabels and not rotation):
            ax.set_zticklabels(labels = ticklabels)
                               
        elif (ticklabels and rotation):
            ax.set_zticklabels(labels = ticklabels, rotation = rotation)




def remove_spines_n_ticks(axs, axis = 'y'):
    '''
    This function removes the spines, mayor and minor ticks of the given axis.
    Parameters:
        axs: list of axes objects
        axis: str, default is y, defines the axis
    '''
    
    if axis == 'y':
        for ax in axs:
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis = 'y', size = 0)
            ax.tick_params(axis = 'y', which = 'minor', size = 0)
    
    elif axis == 'x':
        for ax in axs:
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(axis = 'x', size = 0)
            ax.tick_params(axis = 'x', which = 'minor', size = 0)
        
 
# gradient bar functions
        
def hex2rgb(hex, normalize=False):
    h = hex.strip('#')
    rgb = np.asarray(list(int(h[i:i + 2], 16) for i in (0, 2, 4)))
    return rgb

def draw_rectangle_gradient(ax, x1, y1, width, height, color1='white', color2='blue', alpha1=0.0, alpha2=0.5, n=100):
    # https://stackoverflow.com/questions/24976471/matplotlib-rectangle-with-color-gradient-fill
    
    color1 = hex2rgb(color1) / 255.  # np array
    color2 = hex2rgb(color2) / 255.  # np array


    # Create an array of the linear gradient between the two colors
    gradient_colors = []
    for segment in np.linspace(0, width, n):
        interp_color = [(1 - segment / width) * color1[j] + (segment / width) * color2[j] for j in range(3)]
        interp_alpha = (1 - segment / width) * alpha1 + (segment / width) * alpha2
        gradient_colors.append((*interp_color, interp_alpha))
    for i, color in enumerate(gradient_colors):
        ax.add_patch(plt.Rectangle((x1 + width/n * i, y1), width/n, height, color=color, linewidth=0, zorder=0))
    return ax



def simple_beeswarm(y, nbins=None, width = 1.0):
    """
    https://stackoverflow.com/questions/36153410/how-to-create-a-swarm-plot-with-matplotlib
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """
    y = np.asarray(y)
    if nbins is None:
        # nbins = len(y) // 6
        nbins = np.ceil(len(y) / 6).astype(int)

    # Get upper bounds of bins
    x = np.zeros(len(y))

    nn, ybins = np.histogram(y, bins=nbins)
    nmax = nn.max()

    #Divide indices into bins
    ibs = []#np.nonzero((y>=ybins[0])*(y<=ybins[1]))[0]]
    for ymin, ymax in zip(ybins[:-1], ybins[1:]):
        i = np.nonzero((y>ymin)*(y<=ymax))[0]
        ibs.append(i)

    # Assign x indices
    dx = width / (nmax // 2)
    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(yy)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x
