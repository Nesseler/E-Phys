# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 19:20:26 2023

@author: nesseler
"""

import pandas as pd
import scipy as sc
import matplotlib as mtl
import matplotlib.pyplot as plt
import numpy as np
from useful_functions import calc_time_series, butter_filter
import time

test = sc.io.loadmat('C:/Users/nesseler/Desktop/local E-Phys/AMY-20230905-E004-cc_IF_06.mat')


potential_df = pd.DataFrame()
current_df = pd.DataFrame()
time_in_ms_df = pd.DataFrame()

for i in range(45):
    current_df[i] = pd.DataFrame(np.transpose(test['Trace_2_6_' + str(i+1) + '_2'])[1])
    potential_df[i] = pd.DataFrame(np.transpose(test['Trace_2_6_' + str(i+1) + '_1'])[1])
    time_in_ms_df[i] = pd.DataFrame(np.transpose(test['Trace_2_6_' + str(i+1) + '_1'])[0])


# %%

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

    #calculate the first derivate dv/dt
    dv = np.diff(v)
    dt = np.diff(t)
    dvdt = dv / dt
    
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



test_step = potential_df[21] * 1e3
test_AP = test_step[5000:5401]


phaseplane, ppaxs = plt.subplots(1,2, layout = 'constrained')

plot_voltage_v_time(test_AP, ppaxs[0], {'c':'k', 'label':'unfiltered'})

phase_plane_plot(test_AP, ppaxs[1], {'c':'k', 'label':'unfiltered'})





test_AP_filtered = butter_filter(test_AP, order=1, cutoff=1e3, sampling_rate=20e3)

plot_voltage_v_time(test_AP_filtered, ppaxs[0], {'c':'r', 'label':'filtered'})

phase_plane_plot(test_AP_filtered[20:], ppaxs[1], {'c':'r', 'label':'filtered'})

ppaxs[0].legend(loc = 'upper right')
    

# %%

# filter orders figure


test_step = potential_df[21] * 1e3
test_AP = test_step[5000:5301]


phaseplane, ppaxs = plt.subplots(1,2, layout = 'constrained')

#plot_voltage_v_time(test_AP, ppaxs[0], {'c':'k', 'label':'unfiltered'})

#phase_plane_plot(test_AP, ppaxs[1], {'c':'k', 'label':'unfiltered'})


filter_orders = 5

cmap = plt.get_cmap('plasma', filter_orders+1)

for idx, order in enumerate(np.linspace(0, filter_orders, filter_orders+1)):
      
    test_AP_filtered = butter_filter(test_AP, order=order, cutoff=1e3, sampling_rate=20e3)
    
    plot_voltage_v_time(test_AP_filtered, ppaxs[0], {'c':cmap(idx), 'label':f'filtered ({int(order)} order)'})
    
    phase_plane_plot(test_AP_filtered[20:], ppaxs[1], {'c':cmap(idx), 'label':f'filtered ({int(order)} order)'})


# Normalizer
norm = mtl.colors.Normalize(vmin=0, vmax=filter_orders)
  
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
  

cbar = phaseplane.colorbar(sm, ticks=np.linspace(0.4, filter_orders-0.4, 5+1), label = "Butterworth filter order")
cbar.ax.set_yticklabels(np.linspace(0, filter_orders, 5+1,dtype=np.int16))


# TODO
## new figure with insets at threshold, peak & repolarisation



# %%

test_step = potential_df[21] * 1e3
test_AP = test_step[5000:5401]
t = calc_time_series(test_AP)
v, dvdt = phase_plane_plot(test_AP)



phaseplane, ppaxs = plt.subplots(1,2, layout = 'constrained')



    


#plot_voltage_v_time(test_AP, ppaxs, {'c':'k', 'label':'unfiltered'})


plot_w_colorcode(t, test_AP, t)


ppaxs[0].set_ylim([-100, 20])
ppaxs[0].set_xlim(t.min(), t.max())


#plot_w_colorcode(test_AP[1:], dvdt, t[1:], phaseplane, ppaxs[1])

#phase_plane_plot(test_AP, ppaxs[1], {'c':'k', 'label':'unfiltered'})


#ppaxs[1].plot(v[1:], dvdt)
ppaxs[1].set_ylim([-50, 250])

ppaxs[1].set_xlim([-100, 20])

ppaxs[1].set_ylabel('Rate of membrane potential change\n[mV/ms]')
ppaxs[1].set_xlabel('Membrane potential [mV]')



# %%

plt.style.use('default')

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




test_step = potential_df[14] * 1e3
test_AP = test_step
test_AP_filtered = butter_filter(test_AP, order=1, cutoff=1e3, sampling_rate=20e3)
test_AP_filtered = test_AP_filtered[5150:5401]
t = calc_time_series(test_AP_filtered)
v, dvdt = phase_plane_plot(test_AP_filtered)



phaseplane, ppaxs = plt.subplots(1,2, layout = 'constrained')

t = t
v = test_AP_filtered
dvdt = dvdt
data_fc = t           #data for colorcode
cmap = 'viridis_r'

norm = mtl.colors.CenteredNorm(0, dvdt.max())

lc = get_colorcode(t, v, dvdt, norm = norm)

line = ppaxs[0].add_collection(lc)
phaseplane.colorbar(line, ax=ppaxs[1])
    


lc = get_colorcode(v[1:], dvdt, dvdt, norm = norm)

line = ppaxs[1].add_collection(lc)



ppaxs[0].set_ylim([-100, 20])
ppaxs[0].set_xlim(t.min(), t.max())
ppaxs[0].grid(False)

ppaxs[1].set_ylim([-100, 250])
ppaxs[1].set_xlim([-100, 20])
ppaxs[1].set_ylabel('Rate of membrane potential change\n[mV/ms]')
ppaxs[1].set_xlabel('Membrane potential [mV]')

ppaxs[1].grid(False)


