# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 19:46:34 2024

@author: nesseler
"""

# correlation between n_spikes and v_rest

import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sbn

# custom directories & parameters
from directories_win import cell_descrip_file, figure_dir
from functions_plotting import get_colors, set_font_sizes, get_figure_size, save_figures


cells_df = pd.read_excel(cell_descrip_file, index_col='cell_ID')


# %% plotting only spiking

darkmode_bool = True

colors_dict = get_colors(darkmode_bool)

set_font_sizes()

fig_size = get_figure_size()


jointgrid = sbn.jointplot(data = cells_df.query('n_spikes > 0'),
                          x = 'v_rest',
                          y = 'n_spikes',kind = 'reg')

# jointgrid.fig.set_figwidth(fig_size[0])
# jointgrid.fig.set_figheight(fig_size[1])

save_figures(jointgrid, 'correlation-n_spikes&v_rest-regression', figure_dir, darkmode_bool)


# %% plotting only spiking

# set_font_sizes()

jointgrid = sbn.jointplot(data = cells_df,
                          x = 'v_rest',
                          y = 'n_spikes',
                          hue = 'activity')

# jointgrid.fig.set_figwidth(fig_size[0])
# jointgrid.fig.set_figheight(fig_size[1])

save_figures(jointgrid, 'correlation-n_spikes&v_rest', figure_dir, darkmode_bool)


plt.show()





