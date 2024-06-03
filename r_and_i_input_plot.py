# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 17:03:54 2024

@author: nesseler
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from functions.functions_plotting import get_colors, get_figure_size

darkmode_bool = True
colors_dict, region_colors = get_colors(darkmode_bool)

v_hold = -85 # mV


rinput_low = 100 #MOhm
taumem_low = 7.5 # ms

v_sagx = -160
v_sag = -130 #mV
v_thres = -65 #mV
v_max = -30 #mV
v_maxx = 0 #mV

# ΔI = ΔU / R_input
deltaI_sag = (v_sag - v_hold) / rinput_low * 1e3 # pA
print(deltaI_sag)

r_inputs = np.arange(1, 1600, 1)


deltaIs = pd.DataFrame()


for i, r_input in enumerate(r_inputs):
    deltaIs.at[r_input, 'sagx'] = (v_sagx - v_hold) / r_input * 1e3
    deltaIs.at[r_input, 'sag'] = (v_sag - v_hold) / r_input * 1e3
    deltaIs.at[r_input, 'thres'] = (v_thres - v_hold) / r_input * 1e3
    deltaIs.at[r_input, 'max'] = (v_max - v_hold) / r_input * 1e3
    deltaIs.at[r_input, 'maxx'] = (v_maxx - v_hold) / r_input * 1e3



# low

parameter_low = {'i_start' : -450,
                  'i_delta' : 50,
                  'n_steps' : 22}

iInput_low = np.arange(parameter_low['i_start'], parameter_low['i_start'] + (parameter_low['i_delta'] * parameter_low['n_steps']), parameter_low['i_delta'])

plt.hlines(y = iInput_low,
            xmin = 100, 
            xmax = 250,
            colors = colors_dict['primecolor'],
            alpha = 0.5)


# mid

parameter_mid = {'i_start' : -200,
                  'i_delta' : 20,
                  'n_steps' : 26}

iInput_mid = np.arange(parameter_mid['i_start'], parameter_mid['i_start'] + (parameter_mid['i_delta'] * parameter_mid['n_steps']), parameter_mid['i_delta'])

plt.hlines(y = iInput_mid,
            xmin = 250,     # 200
            xmax = 1000,     # 350
            colors = colors_dict['primecolor'],
            alpha = 0.5)

# high

parameter_high = {'i_start' : -150,
                  'i_delta' : 20,
                  'n_steps' : 21}

iInput_high = np.arange(parameter_high['i_start'], parameter_high['i_start'] + (parameter_high['i_delta'] * parameter_high['n_steps']), parameter_high['i_delta'])

plt.hlines(y = iInput_high,
            xmin = 1000,     # 350
            xmax = 1500,     # 500
            colors = colors_dict['primecolor'],
            alpha = 0.5)

# # 4

# parameter_4 = {'i_start' : -100,
#                   'i_delta' : 10,
#                   'n_steps' : 26}

# iInput_4 = np.arange(parameter_4['i_start'], parameter_4['i_start'] + (parameter_4['i_delta'] * parameter_4['n_steps']), parameter_4['i_delta'])

# plt.hlines(y = iInput_4,
#             xmin = 500, 
#             xmax = 1000,
#             colors = colors_dict['primecolor'],
#             alpha = 0.5)

# # xhigh

# parameter_xhigh = {'i_start' : -50,
#                   'i_delta' : 5,
#                   'n_steps' : 26}

# iInput_xhigh = np.arange(parameter_xhigh['i_start'], parameter_xhigh['i_start'] + (parameter_xhigh['i_delta'] * parameter_xhigh['n_steps']), parameter_xhigh['i_delta'])

# plt.hlines(y = iInput_xhigh,
#             xmin = 1000, 
#             xmax = 1500,
#             colors = colors_dict['primecolor'],
#             alpha = 0.5)






plt.plot(r_inputs, deltaIs['maxx'], '-', label = f'{v_maxx} mV', c = 'r')

plt.plot(r_inputs, deltaIs['max'], '-', label = f'{v_max} mV (~max AP freq.)', c = colors_dict['color1'], linestyle = 'dashed')

plt.plot(r_inputs, deltaIs['thres'], '-', label = f'{v_thres} mV (~threshold)', c = colors_dict['color1'], linestyle = 'solid')

plt.axhline(y = 0, color = colors_dict['color2'], lw = 1, label = '-85 mV (holding)')

plt.plot(r_inputs, deltaIs['sag'], '-', label = f'{v_sag} mV (sag potential)', c = colors_dict['color1'], linestyle = 'dashed')

plt.plot(r_inputs, deltaIs['sagx'], '-', label = f'{v_sagx} mV', c = 'r')


plt.xticks(np.arange(0, 1500+1, 500))
plt.xticks(np.arange(0, 1500+1, 100), minor = True)

plt.xlim([50, 1550])
plt.ylim(-500, 750)

plt.xlabel('R_input [MOhm]')
plt.ylabel('I_input [pA]')

plt.legend()

plt.show()

# %%

from parameters.directories_win import quant_data_dir, cell_descrip_dir, table_file, figure_dir

from os.path import join

import seaborn as sbn

passiv_properties_df = pd.read_excel(join(cell_descrip_dir, 'ccIF-passiv_properties.xlsx'), index_col = 'cell_ID')

# get cell IDs
cell_IDs = passiv_properties_df.index.to_list()

MetaData = pd.read_excel(table_file,
                      sheet_name="MetaData",
                      index_col='cell_ID')

MetaData = MetaData.loc[cell_IDs, :]

plt_df = pd.concat([passiv_properties_df, MetaData], axis = 1)


# initialize figure
fig_comb, axs_comb = plt.subplots(nrows = 2,
                                  ncols = 1,
                                  #figsize = get_figure_size(),
                                  layout = 'constrained',
                                  sharex= True,
                                  height_ratios=[3, 6])

ax = axs_comb[0]

violin = sbn.violinplot(data = plt_df,
                        y = 'Region',
                        x = 'r_input',
                        inner = 'quart',
                        linewidth = 1,
                        ax = ax,
                        size = 0.9,
                        order = ['BAOT/MeA', 'MeA', 'BAOT'],
                        orient = 'h')

for l in violin.lines:
    l.set_color(colors_dict['primecolor'])

for violin in violin.collections:
    violin.set_edgecolor(colors_dict['primecolor'])
    violin.set_facecolor('None')

swarm = sbn.swarmplot(data = plt_df,
                      y = 'Region',
                      x = 'r_input', 
                      ax = ax,
                      hue = 'Region', 
                      palette = region_colors,
                      size = 5,
                      color = colors_dict['primecolor'],
                      order = ['BAOT/MeA', 'MeA', 'BAOT'],
                      orient = 'h')

ax.legend().set_visible(False)



# 1
axs_comb[1].hlines(y = -425, xmin = 100, xmax = 150, colors = colors_dict['primecolor'], alpha = 0.5)

# 2
axs_comb[1].hlines(y = -325, xmin = 150, xmax = 200, colors = colors_dict['primecolor'], alpha = 0.5)

# 3
axs_comb[1].hlines(y = -250, xmin = 200, xmax = 250, colors = colors_dict['primecolor'], alpha = 0.5)

# 4
axs_comb[1].hlines(y = -200, xmin = 250, xmax = 300, colors = colors_dict['primecolor'], alpha = 0.5)

# 5
axs_comb[1].hlines(y = -175, xmin = 300, xmax = 350, colors = colors_dict['primecolor'], alpha = 0.5)

# 6
axs_comb[1].hlines(y = -150, xmin = 350, xmax = 400, colors = colors_dict['primecolor'], alpha = 0.5)

# 7
axs_comb[1].hlines(y = -125, xmin = 400, xmax = 500, colors = colors_dict['primecolor'], alpha = 0.5)

# 8
axs_comb[1].hlines(y = -100, xmin = 500, xmax = 600, colors = colors_dict['primecolor'], alpha = 0.5)

# 9
axs_comb[1].hlines(y = -85, xmin = 600, xmax = 700, colors = colors_dict['primecolor'], alpha = 0.5)

# 10
axs_comb[1].hlines(y = -75, xmin = 700, xmax = 850, colors = colors_dict['primecolor'], alpha = 0.5)


axs_comb[1].plot(r_inputs, deltaIs['maxx'], '-', label = f'{v_maxx} mV', c = 'r')

axs_comb[1].plot(r_inputs, deltaIs['max'], '-', label = f'{v_max} mV (~max AP freq.)', c = colors_dict['color1'], linestyle = 'dashed')

axs_comb[1].plot(r_inputs, deltaIs['thres'], '-', label = f'{v_thres} mV (~threshold)', c = colors_dict['color1'], linestyle = 'solid')

axs_comb[1].axhline(y = 0, color = colors_dict['color2'], lw = 1, label = '-85 mV (holding)')

axs_comb[1].plot(r_inputs, deltaIs['sag'], '-', label = f'{v_sag} mV (sag potential)', c = colors_dict['color1'], linestyle = 'dashed')

axs_comb[1].plot(r_inputs, deltaIs['sagx'], '-', label = f'{v_sagx} mV', c = 'r')


axs_comb[1].set_xticks(np.arange(0, 1500+1, 500))
axs_comb[1].set_xticks(np.arange(0, 1500+1, 100), minor = True)
axs_comb[1].set_xlim([50, 1550])
axs_comb[1].set_xlabel('R_input [MOhm]')

axs_comb[1].set_yticks(np.arange(-500, 750+1, 250))
axs_comb[1].set_yticks(np.arange(-500, 750+1, 50), minor = True)
axs_comb[1].set_ylim(-500, 750)
axs_comb[1].set_ylabel('I_input [pA]')

plt.legend()

plt.show()











