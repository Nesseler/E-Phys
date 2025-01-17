# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 14:31:40 2023

@author: nesseler
"""

import numpy as np
import pandas as pd
import matplotlib as mtl
import matplotlib.pyplot as plt

from functions.functions_useful import calc_time_series, calc_dvdt



def phase_plane_plot(data, axis=None, plot_dict={'c':'k'}, sampling_rate=20e3):
    """
    Calculate and generate phase plane plot for single action potential
    Time is dealt in ms
    
    Parameters:
        data : One-dimensional array with voltage in mV
        axes : Axis to plot phase plane plot on
        sampling_rate : Sampling rate in Hz. Default is 20 kHz.
    Returns:
        v
        dvdt
    """
    v = data
    t = calc_time_series(v, sampling_rate, scale='ms')

    dvdt = calc_dvdt(v,t)
    
    if axis != None:
        axis.plot(v[1:], dvdt, **plot_dict)
        axis.set_ylabel('Rate of membrane potential change\n[mV/ms]')
        axis.set_xlabel('Membrane potential [mV]')

    return v, dvdt



def plot_voltage_v_time(data, axis, plot_dict, v_range = [-100, 20], sampling_rate = 20e3, scale='ms'):
    '''
    Plots the voltage on the given axis with calculated time series.
    Parameters:
        data : One-dimensional array with voltage in mV
        axes : Axis to plot phase plane plot on
        plot_dict : Plotting dictionary that is passed to the matplotlib plot
                    function to specify the plotted line.
        v_range : Range of voltage values covered. Used for the y axis limits.
                  Default is -100 mV to + 20 mV.
        sampling_rate : Sampling rate in Hz. Default is 20 kHz.
        scale : {'ms', 's'} Time scale in which the time series is being 
                calculated. 'ms' Milliseconds is default.
    '''
    #calculate time series
    t = calc_time_series(data, sampling_rate, scale)

    axis.plot(t, data, **plot_dict)

    axis.set_xlabel(f'Time [{scale}]')
    axis.set_xlim([t[0], t[-1]])

    axis.set_ylim(v_range)
    axis.set_ylabel('Membrane potential [mV]')


def get_colorcode(x, y, data_fc, norm=None, cmap='seismic', plot_dict={'c':'k'}):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    return_bool = False

    if norm is None:
        # Create a continuous norm to map from data points to colors
        norm = plt.Normalize(data_fc.min(), data_fc.max())
        return_bool = True
        
        
    lc = mtl.collections.LineCollection(segments, cmap=cmap, norm=norm)
    # Set the values used for colormapping
    lc.set_array(data_fc)
    
    return (norm, lc) if return_bool else lc



def plot_vt_n_dvdtv_colorcoded(v, v_range = [-100, 20], sampling_rate = 20e3, scale='ms', cmap='seismic'):
    t = calc_time_series(v, sampling_rate, scale)
    v, dvdt = phase_plane_plot(v,sampling_rate=sampling_rate)
    data_fc = dvdt
    
    cc_phaseplane, cc_ppaxs = plt.subplots(1,2, layout = 'constrained')
    
    norm = mtl.colors.CenteredNorm(0, data_fc.max())

    lc = get_colorcode(t, v, data_fc, norm = norm, cmap=cmap)
    
    line = cc_ppaxs[0].add_collection(lc)
    cc_phaseplane.colorbar(line, ax=cc_ppaxs[1])
    
    lc = get_colorcode(v[1:], data_fc, data_fc, norm = norm, cmap=cmap)
    
    line = cc_ppaxs[1].add_collection(lc)
    
    cc_ppaxs[0].grid(False)
    cc_ppaxs[0].set_xlim(t.min(), t.max())
    cc_ppaxs[0].set_xlabel(f'Time [{scale}]')
    cc_ppaxs[0].set_ylim(v_range)
    cc_ppaxs[0].set_ylabel('Membrane potential [mV]')
    
    cc_ppaxs[1].set_ylim([-100, 250])
    cc_ppaxs[1].set_xlim(v_range)
    cc_ppaxs[1].set_ylabel('Rate of membrane potential change\n[mV/ms]')
    cc_ppaxs[1].set_xlabel('Membrane potential [mV]')
    
    cc_ppaxs[1].grid(False)
    
    





from functions.functions_extractspike import extract_spike


def calc_vmem_at_spiketrain(t, v, dvdt, spike_idc, min_ISI, SR):
    '''
    Function calculates the average membrane voltage of a voltage trace 
    which includes spikes. The datapoints of these spikes are excluded by
    the dvdt_threshold that is to be crossed.
    Parameters:
        v : Voltage trace.
        dvdt : First derivate of voltage.
        spike_idc : List of spike indices.
        min_ISI : Minimum ISI of spiketrain in ms.
        SR : Sampling Rate. In Hz.
        dvdt_threshold : First derivative threshold that is to be crossed.
    Returns:
        v_mem : Average membran voltage without spikes.
    '''

    ### how to get the v_mem at burst ###
    
    plotting_bool = True
    
    if plotting_bool:
        plt.plot(v, dvdt)
        plt.ylim([-60, 200])
        plt.xlim([-100, 60])
        plt.show()
        plt.pause(0.1)
    
    # exclude datapoints of spike with the dvdt threshold, then take mean
    # define time pre and post of spike to include
    t_pre = 5
    t_post = min_ISI
   
    # loop through each spike and exclude it's datapoints using the dvdt threshold
    for i, idx_spike in enumerate(spike_idc):
        
        # get index of spike relative to spike v array
        idx_spike_rel = int((t_pre * (SR/1e3)))
        
        # get pre and post indices
        idx_pre = idx_spike - idx_spike_rel
        idx_post = idx_spike + int((t_post * (SR/1e3)))
        
        # get v at spike
        spike_v = v[idx_pre:idx_post]
        spike_t = t[idx_pre:idx_post]
        
        # get or calc dvdt
        spike_dvdt = dvdt[idx_pre:idx_post]
        
        # vplot
        if plotting_bool:
            plt.plot(spike_v, spike_dvdt)
            plt.ylim([-20, 20])
            plt.xlim([-100, 60])
            plt.show()
            plt.pause(0.1)
        
        # extract spike
        spike_idc ,_ ,_ ,_ = extract_spike(spike_t, spike_v, spike_dvdt, idx_spike_rel)
           
        # exclude datapoints
        spike_v[spike_idc] = np.nan
        
    # calc v_mem at burst as mean of datapoints without spikes
    v_mem =  np.nanmean(v)
 
    return v_mem


