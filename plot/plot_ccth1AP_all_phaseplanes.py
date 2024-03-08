# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 18:51:47 2024

@author: nesseler
"""

import pandas as pd
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn

# custom directories & parameters
from parameters.directories_win import quant_data_dir, cell_descrip_dir, figure_dir

# custom functions
from functions.functions_plotting import get_figure_size, get_colors, save_figures
from functions.functions_useful import calc_time_series




fstAPs_parameters_df = pd.read_excel(join(cell_descrip_dir, 'ccth1Ap-fst_AP_parameters.xlsx'), index_col='cell_ID')
cell_IDs = fstAPs_parameters_df.index.to_list()



t_post_th = 2
t_pre_th = 0.9


t_df = pd.DataFrame()
v_df = pd.DataFrame()
dvdt_df = pd.DataFrame()

for cell_ID in cell_IDs:
    
    fstAP_df = pd.read_excel(join(quant_data_dir, '1stAP', f'ccth1Ap-1stAP-{cell_ID}.xlsx'), index_col='t')
    
    
    SR_ms = fstAPs_parameters_df.at[cell_ID, 'SR_ms']
    t_pre_idx = int(np.round(fstAPs_parameters_df.at[cell_ID, 't_peaks'] * SR_ms)) - int(t_pre_th * SR_ms)
    t_post_idx = t_pre_idx + int((t_post_th + t_pre_th * 2) * SR_ms) + 1
    
    
    v = fstAP_df['v'].to_list()
    dvdt = fstAP_df['dvdt'].to_list()
    
    v_df[cell_ID] = v[t_pre_idx:t_post_idx]
    dvdt_df[cell_ID] = dvdt[t_pre_idx:t_post_idx]
    t_df[cell_ID] = np.arange(-t_pre_th, (t_post_th + t_pre_th) + 1/SR_ms, 1/SR_ms)
    
    
# %% plot

    






darkmode_bool = True

colors_dict, _ = get_colors(darkmode_bool)


fig_1stall, ax_1stall = plt.subplots(1, 2, 
                                   layout = 'constrained',
                                   figsize = get_figure_size())


ax_1stall[0].plot(t_df, v_df,
                  c = colors_dict['primecolor'])
    
ax_1stall[1].plot(v_df, dvdt_df,
                  c = colors_dict['primecolor'])

    
    
ax_1stall[0].set_xlim([-t_pre_th, t_post_th + t_pre_th])
ax_1stall[0].set_xticks(np.linspace(0, t_post_th, 3))

    
ax_1stall[0].set_ylim([-100, 60])
ax_1stall[0].set_yticks(np.arange(-100, 60 + 1, 20))
ax_1stall[0].set_yticks(np.arange(-100, 60 + 1, 5), minor = True)

ax_1stall[0].set_ylabel('Membrane potential [mV]')

ax_1stall[0].set_xlabel('Time [ms]')


ax_1stall[1].set_ylabel('Rate of membrane potential change\n[mV/ms]')
ax_1stall[1].set_ylim([-150, 250]) 

ax_1stall[1].set_xlabel('Membrane potential [mV]')
ax_1stall[1].set_xlim([-100, 60])   

[ax.grid(False) for ax in ax_1stall]  
    
    
save_figures(fig_1stall, 'cc_th1AP_all_firstAPs', figure_dir, darkmode_bool)
    
plt.show()
    
