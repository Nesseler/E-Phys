# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 09:56:57 2025

@author: nesseler
"""

from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import figure_dir, table_file
from parameters.PGFs import vc_Erest_parameters

from functions.functions_import import get_traceIndex_n_file, get_vc_data, get_PSCs_steps
from functions.functions_filter import butter_filter

PGF = 'vc_Erest_3min'

# conditions = ['ctrl', 'AP5_NBQX_washin', 'AP5_NBQX', 'AP5_NBQX_GBZ_washin', 'AP5_NBQX_GBZ']
conditions = ['ctrl', 'GBZ_washin', 'GBZ', 'AP5_NBQX_GBZ_washin', 'AP5_NBQX_GBZ']
blkr1 = 'GBZ'
blkr2 = 'AP5_NBQX'

SR = vc_Erest_parameters['SR']

t = vc_Erest_parameters['t']

n_steps = vc_Erest_parameters['n_steps']

# get cell_IDs
# cell_IDs = get_cell_IDs_one_protocol(PGF = PGF + '_' + 'adaEk', sheet_name= 'PGFs_Syn')

cell_ID = 'E-304' # 305 (AP5-GBZ) # 304 (GBZ-AP5)

condition = 'ctrl'


# %% loading

# initialize dicts and dataframes
i_s = pd.DataFrame(columns=conditions, index = t)


print('loading ...')

for condition in tqdm(conditions):

    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')
    
    # get list of steps to include
    idc_steps = get_PSCs_steps(PGF + '_' + condition, cell_ID, sheet_name = 'PGFs_Syn')
    
    # get data with file path & trace index
    i, _, _, _, nsteps_loaded = get_vc_data(file_path, traceIndex, scale='ms')
    
    # concatenate steps
    i = i.flatten('C')
    
    # filter
    # filter all data with 1kHz cutoff
    i = butter_filter(i, order=3, cutoff=1e3, sampling_rate=SR)
    
    # replace first values with nans to eliminate filter artifact
    i[:100] = np.nan
    
    # add nan values if protocol was stopped before complete recording
    if n_steps != nsteps_loaded:
        # expected datapoints
        n_datapoints = vc_Erest_parameters['dur_steps'] * n_steps * SR
        
        # generate filler array
        filler = np.full(shape = (n_datapoints - len(i)), 
                         fill_value = np.nan)
        
        # add filler
        i = np.concatenate([i, filler])
    
    # split back into steps
    i = np.array_split(i, n_steps)
    
    # initialize i
    i_new = np.full(shape = (n_steps, vc_Erest_parameters['dur_steps'] * SR), 
                    fill_value = np.nan)
    
    # include only listed steps
    for i_step in range(n_steps):
        
        if i_step in idc_steps:
            i_new[i_step] = i[i_step]
    
    # concatenate steps
    i = np.array(i_new).flatten('C')
    
    # write to dict
    i_s[condition] = i




# %%

# init plotting
from functions.initialize_plotting import *

fig, axs = plt.subplots(nrows = 1,
                        ncols = 3,
                        figsize = get_figure_size(height = 60),
                        dpi = 300,
                        layout = 'constrained',
                        sharex=True,
                        sharey=True)


for c_idx, condition in enumerate(conditions[0::2]):
    
    # set axis
    ax = axs[c_idx]

    # plot
    ax.plot(t[::100],
            i_s[condition][::100],
            lw = 0.5,
            c = colors_dict['primecolor'],
            zorder = 2)
    

    

# y
ydict = {'ax_min' : -50,
         'ax_max' : 25,
         'pad' : 0.75,
         'step' : 50,
         'stepminor' : 5,
         'label' : ''}

# edit axis
for ax in axs:
    apply_axis_settings(ax, axis = 'y', **ydict)
    

# blkr indicators   
rect_dict = {'height' : 3,
             'fill' : True,
             'linestyle' : 'solid',
             'lw' : 0,
             'zorder' : 0}
    
# add rectanlge patch
for ax in axs[1:]:
    ax.add_patch(Rectangle(xy = (0, 20), #(x_min, y_min)
                           width = 180,
                           color = blkr_colors[blkr1], ##7b7fae
                           **rect_dict))

axs[2].add_patch(Rectangle(xy = (0, 15), #(x_min, y_min)
                           width = 180,
                           color = blkr_colors[blkr2],
                           **rect_dict))

# x
xdict = {'ax_min' : 0,
         'ax_max' : 180,
         'pad' : 3,
         'step' : 60,
         'stepminor' : 5,
         'label' : ''}

# edit axis
for ax in axs:
    apply_axis_settings(ax, axis = 'x', **xdict)
  
# remove y axis spines   
remove_spines_n_ticks(axs[1:], axis = 'y')

# set sup labels
fig.supylabel('Current [pA]')
fig.supxlabel('Time [s]')

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# align labels
fig.align_labels()

# display figure
plt.show()


# %% washin

from functions.functions_plotting import draw_rectangle_gradient




fig, axs = plt.subplots(nrows = 5,
                        ncols = 1,
                        figsize = get_figure_size(width = 160.5),
                        dpi = 300,
                        layout = 'constrained',
                        sharey=True,
                        sharex = True)


for c_idx, condition in enumerate(conditions):
    
    # set axis
    ax = axs[c_idx]
    
    # set axis title
    ax.set_title(condition,
                 fontsize=9, 
                 loc='left',
                 x = 0.015)   

    # plot
    ax.plot(t,
            i_s[condition],
            lw = 0.5,
            c = colors_dict['primecolor'],
            zorder = 2)
    
# add rectanlge patch
# plt.sca(axs[1])


draw_rectangle_gradient(ax = axs[1],
                        x1 = 0,
                        y1 = 20,
                        width = 180,
                        height = 3,
                        color1 = '#000000',
                        color2 = blkr_colors[blkr1],
                        alpha1 = 1.0,
                        alpha2 = 1.0,
                        n = 100)

for ax in axs[2:]:
    ax.add_patch(Rectangle(xy = (0, 20), #(x_min, y_min)
                           width = 180,
                           color = blkr_colors[blkr1], ##7b7fae
                           **rect_dict))
    
draw_rectangle_gradient(ax = axs[3],
                        x1 = 0,
                        y1 = 15,
                        width = 180,
                        height = 3,
                        color1 = '#000000',
                        color2 = blkr_colors[blkr2],
                        alpha1 = 1.0,
                        alpha2 = 1.0,
                        n = 100)

axs[-1].add_patch(Rectangle(xy = (0, 15), #(x_min, y_min)
                            width = 180,
                            color = blkr_colors[blkr2],
                            **rect_dict))

# for ax in axs[::-1]:
# remove x axis spines   
remove_spines_n_ticks(axs[:-1], axis = 'x')
    
# edit axis
for ax in axs:
    apply_axis_settings(ax, axis = 'x', **xdict)
    apply_axis_settings(ax, axis = 'y', **ydict)

    
# set sup labels
fig.supylabel('Current [pA]', fontsize = 12)
fig.supxlabel('Time [s]', fontsize = 12)

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# align labels
fig.align_labels()

# display figure
plt.show()

