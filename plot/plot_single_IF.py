# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 21:42:50 2025

@author: nesseler
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:27:42 2025

@author: nesseler
"""
# import standard packages
from functions.initialize_packages import *

# custom directories & parameters
from parameters.directories_win import cell_descrip_syn_dir

# spike detection
from parameters.parameters import min_peak_prominence, min_peak_distance, min_max_peak_width

# IF characteristation
from parameters.parameters import t_init_inst_freq


# custom functions
from functions.get_cell_IDs import get_cell_IDs_one_protocol
from functions.functions_import import get_cc_data, get_traceIndex_n_file
from functions.functions_filter import merge_filter_split_steps
from functions.functions_useful import calc_dvdt_padded
from functions.functions_extractspike import get_AP_parameters, extract_spike, get_spiketrain_n_ISI_parameter

# PGF specific
from parameters.PGFs import cc_IF_syn_parameters as PGF_parameters
t = PGF_parameters['t']
SR = PGF_parameters['SR']
    

# define protocol
PGF = 'cc_IF'
sheet_name = 'PGFs_Syn'

# get all cell_IDs for cc_rest
cell_IDs = get_cell_IDs_one_protocol(PGF = PGF, sheet_name = sheet_name)

        
# %% define output

IF           = pd.DataFrame(columns=cell_IDs, index = np.arange(-50, 1000, 1, int))
IF_inst      = pd.DataFrame(columns=cell_IDs, index = np.arange(-50, 1000, 1, int))
IF_inst_init = pd.DataFrame(columns=cell_IDs, index = np.arange(-50, 1000, 1, int))

# set index label
for df in [IF, IF_inst, IF_inst_init]:
    df.index.name = 'i_input'

IF_dict      = pd.DataFrame(columns = ['i_rheobase'       , 'idx_rheobase',
                                       'i_maxfreq'        , 'idx_maxfreq'        , 'maxfreq',      
                                       'i_halfmax'        , 'idx_halfmax',
                                       'i_maxinitinstfreq', 'idx_maxinitinstfreq', 'maxinitinstfreq'],
                            index = cell_IDs)

IF_rheobase = pd.DataFrame(columns = ['rheobase_abs', 'rheobase_rel'], 
                           index = cell_IDs)

adaptation = pd.DataFrame(index = cell_IDs)

# set index label
for df in [IF_dict, IF_rheobase, adaptation]:
    df.index.name = 'cell_ID'


# %% check for new cells to be analyzed

# # load anaylsis worksheet
# from parameters.directories_win import table_file
# analyzed = pd.read_excel(table_file,
#                          sheet_name = 'analyzed',
#                          index_col = 'cell_ID')

# # get list of cell_IDs already analyzed
# analyzed_cell_IDs = analyzed.loc[analyzed[PGF].notna()][PGF].index.to_list()

# # redefine cell_IDs list
# cell_IDs = [cell_ID for cell_ID in cell_IDs if cell_ID not in analyzed_cell_IDs]

# # raise error
# if len(cell_IDs) == 0:
#     raise ValueError('Nothing new to analyze!')


# %% initialize plotting and verificaiton plots

# init plotting
from functions.initialize_plotting import *

# # verification plots
# vplots = True
# if vplots:
#     # load plotting functions
#     from analysis.celldescrip_analysis_Syn.plot_analyze_cc_IF_syn import plot_full_IF, plot_IF_step_spike_detection, plot_rheobase, plot_adaptation

    
# %% load

# 223, 276, 300, 315
cell_IDs = ['E-317']

for cell_ID in tqdm(cell_IDs):

    # load IF protocol
    # get the traceIndex and the file path string for data import functions
    traceIndex, file_path = get_traceIndex_n_file(PGF, cell_ID, sheet_name = sheet_name)
    
    # get data with file path & trace index
    i, v, _, _, n_steps = get_cc_data(file_path, traceIndex, scale='s')
    
    # filter steps
    v = merge_filter_split_steps(v, SR)
    
    # limit time and voltage to stimulation period
    idc_stim = PGF_parameters['idc_stim']
    
    # add 100 ms after stimulation period to include APs at the end of the step
    t_additional = 50 #ms
    idc_additional = np.arange(idc_stim[-1] +1, idc_stim[-1] +1 + (t_additional * (SR/1e3)), dtype = int)
    idc_stim = np.append(idc_stim, idc_additional)
    
    v_full = np.copy(v)
    t_full = np.copy(t)
    v_stim = v_full[:, idc_stim]
    t_stim = t_full[idc_stim]
    
    # construct current array
    from functions.functions_constructors import construct_I_array
    i_calc, i_input = construct_I_array(cell_ID, n_steps,
                                        PGF = PGF,
                                        sheet_name = sheet_name,
                                        parameters = PGF_parameters)
    
    # get i_hold
    i_hold = i_input[0]
    
    # if vplots:
    
        
    
# %%


    
    
# %%
cell_ID
t = t_full
v = v_full[::6]
i = i_calc[::6]
i_input = i_input[::6]
gradient = True

# %%

# y
ydict_v = {'ax_min' : -100,
           'ax_max' : 60,
           'pad' : None,
           'step' : 50,
           'stepminor' : 5,
           'label' : '',
           'limits_n_0' : True}

# x
xdict_tfull = {'ax_min' : 0,
               'ax_max' : 1500,
               'pad' : 10,
               'step' : 1500,
               'stepminor' : 50,
               'label' : ''}

# get number of steps
n_steps = v.shape[0]

# init figure and axes
fig, axs = plt.subplots(nrows = 2,
                        ncols = 1,
                        layout = 'constrained',
                        figsize = get_figure_size(width = 80, height = 65),
                        sharex = True,
                        height_ratios = [1, 5],
                        dpi = 300)

# # set axis title
# axs[0].set_title(f'{cell_ID} cc_IF',
#                  fontsize=9, 
#                  loc='left',
#                  x = 0.02)

# color code
if gradient:
    # specify linecollection settings
    lc_dict = {'lw' : 0.3,
               'linestyle' : 'solid',
               'array' : i_input,
               'cmap' : 'plasma_r',
               'norm' : plt.Normalize(vmin=-100,vmax=700)}
else:
    # specify linecollection settings
    lc_dict = {'lw' : 0.5,
               'linestyle' : 'solid',
               'color' : colors_dict['primecolor']}
    
    
# define line collections
v_collection = LineCollection([np.column_stack([t, v[step]]) for step in range(n_steps)],
                              **lc_dict)



i_collection = LineCollection([np.column_stack([t, i[step]]) for step in range(n_steps)], 
                              **lc_dict) 

# voltage
ax = axs[1]

# add line collection
ax.add_collection(v_collection)

# apply_axis_settings(ax, axis = 'y', **ydict_v)

# current
ax = axs[0]

# add line collection
ax.add_collection(i_collection)

# y
ydict_i = {'ax_min' : -100,
           'ax_max' : 1000,
           'pad' : 10,
           'step' : 1000,
           'stepminor' : 100,
           'label' : '',
           'start_at_0' : True}

# apply_axis_settings(ax, axis = 'y', **ydict_i)

# x 
# apply_axis_settings(axs[-1], axis = 'x', **xdict_tfull)
remove_spines_n_ticks([axs[0]], axis = 'x')

remove_spines_n_ticks([axs[0]], axis = 'y')

remove_spines_n_ticks([axs[1]], axis = 'y')
remove_spines_n_ticks([axs[1]], axis = 'x')


[ax.set_yticklabels([]) for ax in axs]

axs[1].set_xticklabels([])

axs[1].set_xlim([0, 1500])
axs[1].set_ylim([-100, 60])

axs[0].set_ylim([-100, 1000])

# remove spines
[ax.spines[spine].set_visible(False) for ax in axs for spine in ['top', 'right']]

# align labels
fig.align_labels()

# create saving path and save
from parameters.directories_win import figure_dir
path_fig = join(figure_dir, 'temp_figs')
save_figures(fig, f'{cell_ID}-cc_IF', path_fig, darkmode_bool, figure_format='both')

# display figure
plt.show()