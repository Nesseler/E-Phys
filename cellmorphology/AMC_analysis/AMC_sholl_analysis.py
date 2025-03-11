# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:33:42 2025

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# get MetaData
from cellmorphology.AMC_analysis.AMC_analysis_import import get_cells_list
MetaData = get_cells_list()

# get cell_IDs to be analyzed
cell_IDs = MetaData.query('coordinates == "Yes" & paths_checked == "Yes"').index.to_list()

# define sholl radius step size
sholl_step_size = 1 # µm
max_sholl_radius = 590 # µm
radius = np.arange(0, max_sholl_radius + sholl_step_size, sholl_step_size)

# define types to be analyzed
neurite_types = ['neurites', 'dendrites', 'axons']

# create combined sholl profile dataframe
all_sholl_profiles_template = pd.DataFrame(index = radius,
                                           columns = cell_IDs)

# rename index column in template
all_sholl_profiles_template.index.name = 'Radius'

# initialize dictionary for combined sholl profiles
all_sholl_profiles = {'neurites'  : all_sholl_profiles_template.copy(),
                      'dendrites' : all_sholl_profiles_template.copy(),
                      'axons'     : all_sholl_profiles_template.copy()}

# create dataframe for average profiles
avgs = ['mean', 'median', 'std']
avg_columns = [f'{avg}-sholl_profile-{ntype}' for ntype in neurite_types for avg in avgs]
avg_sholl_profiles = pd.DataFrame(columns = avg_columns,
                                  index = radius)
avg_sholl_profiles.index.name = 'Radius'

# create dataframe for sholl metrics
metrics = ['max_intersections', 'critical_radius', 'enclosing_radius']
metrics_columns = [f'{metric}-{ntype}' for ntype in neurite_types for metric in metrics]
sholl_metrics = pd.DataFrame(columns = metrics_columns,
                             index = cell_IDs)
sholl_metrics.index.name = 'cell_ID'


# %% load data

# get directory
from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_sholl_dir

print('loading ...')

for cell_ID in tqdm(cell_IDs):

    # get neurites sholl profile for cell
    sholl_profile_neurites = pd.read_csv(join(AMCs_sholl_dir, 'sholl_profiles_neurites', f'{cell_ID}_profile.csv'), 
                                         index_col = 'Radius',
                                         usecols = ['Radius', 'Inters.'])
    
    # rename column
    sholl_profile_neurites.rename(columns = {'Inters.' : cell_ID},
                                  inplace = True)
    
    # get axon sholl profile for cell
    try:
        sholl_profile_axons = pd.read_csv(join(AMCs_sholl_dir, 'sholl_profiles_axons', f'{cell_ID} [axon]_profile.csv'), 
                                          index_col='Radius',
                                          usecols = ['Radius', 'Inters.'])
        
        # rename column
        sholl_profile_axons.rename(columns = {'Inters.' : cell_ID},
                                   inplace = True)
    
        # calculate dendrites sholl profile
        # as difference between neurites and axons
        sholl_profile_dendrites = sholl_profile_neurites.sub(sholl_profile_axons, fill_value = 0)
        
    except FileNotFoundError:
        # if file does not exsist, no axons was labeled
        # axon profile will be set as nans
        sholl_profile_axons = all_sholl_profiles_template.copy()[[cell_ID]]
        
        # set sholl profile of dendrites as profile of neurites
        sholl_profile_dendrites = sholl_profile_neurites
        
    # write dendrites profile to csv
    sholl_profile_dendrites.rename(columns = {cell_ID : 'Inters.'}).to_csv(join(AMCs_sholl_dir, 'sholl_profiles_dendrites', f'{cell_ID} [dendrites]_profile.csv'),
                                                                           index_label = 'Radius')
    
    # write profiles to dictionary
    all_sholl_profiles['neurites'].loc[:, cell_ID] = sholl_profile_neurites[cell_ID]
    all_sholl_profiles['dendrites'].loc[:, cell_ID] = sholl_profile_dendrites[cell_ID]
    all_sholl_profiles['axons'].loc[:, cell_ID] = sholl_profile_axons[cell_ID]

# save collection of profiles
for ntype in neurite_types:
    all_sholl_profiles[ntype].to_excel(join(AMCs_sholl_dir, 'sholl_combined_profiles', f'sholl_profiles-{ntype}.xlsx'),
                                       index_label = 'Radius')
    
    
# %% calc average profiles

for ntype in neurite_types:
    # mean
    avg_sholl_profiles[f'mean-sholl_profile-{ntype}'] = all_sholl_profiles['neurites'].mean(axis = 1)
    
    # median
    avg_sholl_profiles[f'median-sholl_profile-{ntype}'] = all_sholl_profiles['neurites'].median(axis = 1)
    
    # std
    avg_sholl_profiles[f'std-sholl_profile-{ntype}'] = all_sholl_profiles['neurites'].std(axis = 1)

# save average profiles
avg_sholl_profiles.to_excel(join(AMCs_sholl_dir, 'sholl_combined_profiles', 'sholl_profiles-avg.xlsx'),
                            index_label = 'Radius')


# %% get sholl metrics per type

for ntype in neurite_types:
    # max number of intersections
    sholl_metrics.loc[:, f'max_intersections-{ntype}'] = all_sholl_profiles[ntype].max(axis = 0)
    
    # critical radius
    sholl_metrics.loc[:, f'critical_radius-{ntype}'] = all_sholl_profiles[ntype].dropna(axis = 1, how = 'all').idxmax(axis = 0)
    
    # enclosing radius
    sholl_metrics.loc[:, f'enclosing_radius-{ntype}'] = all_sholl_profiles[ntype].count(axis = 0) -1
    
# replace -1 values in axons with nans
sholl_metrics.astype(float).replace({-1: np.nan}, inplace = True)
    
# save sholl metrics
from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_metrics_dir
sholl_metrics.to_excel(join(AMCs_metrics_dir, 'sholl_metrics.xlsx'),
                       index_label = 'cell_ID')


# %% create plots

# initialize plotting packages
from cellmorphology.AMC_analysis.initialize_AMC_cellmorph_plotting import *

print('plotting ...')

for cell_ID in tqdm(cell_IDs):

    # init figure and axes
    fig, ax = plt.subplots(nrows = 1,
                           ncols = 1,
                           layout = 'constrained',
                           figsize = get_figure_size(width = 160, height = 100),
                           dpi = 300)
    
    # set figure title
    fig.suptitle(f'{cell_ID} sholl profiles',
                 fontsize = 9)
    
    # plot sholl profiles
    for ntype in neurite_types:
        ax.plot(radius,
                all_sholl_profiles[ntype][cell_ID],
                color = neurite_color_dict[ntype],
                lw = 1,
                label = ntype)
    
    # set legend
    ax.legend(title = 'neurite type', 
              frameon = False,
              loc = 'upper right')
    
    # y axis
    ydict = {'ax_min' : 0,
             'ax_max' : 50,
             'pad' : None,
             'step' : 5,
             'stepminor' : 1,
             'label' : 'Number of intersections [#]'}
    
    apply_axis_settings(ax, axis = 'y', **ydict)
    
    # x axis
    xdict = {'ax_min' : 0,
             'ax_max' : 400,
             'pad' : None,
             'step' : 100,
             'stepminor' : 10,
             'label' : 'Radius [µm]'}
    
    apply_axis_settings(ax, axis = 'x', **xdict)
    
    # remove spines
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
   
    # create saving path and save
    from cellmorphology.AMC_analysis.AMC_analysis_directories import AMCs_shollplots
    save_figures(fig, 
                 f'{cell_ID}-sholl_profiles', 
                 AMCs_shollplots, 
                 darkmode_bool, 
                 figure_format='png')
    
    # display figure
    plt.show()