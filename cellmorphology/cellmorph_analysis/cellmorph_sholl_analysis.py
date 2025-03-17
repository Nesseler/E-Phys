# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 10:06:26 2025

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# get directories
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_sholl_dir, cellmorph_analysis_dir, cellmorph_metrics_dir

# get files in traces folder
from functions.get_onlyfiles_list import get_onlyfiles_list
onlyfiles = get_onlyfiles_list(cellmorph_sholl_dir + '/sholl_tables_neurites')

# create filenames dict
filenames = {'E-' + filename[-15:-12] : filename for filename in onlyfiles}

# get cell_IDs from list of all coordinates filenames
cell_IDs = list(filenames.keys())
# cell_IDs = ['E-111', 'E-137', 'E-147', 'E-162']
cell_IDs = cell_IDs[:-20] # for leonie NWG poster

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# set neurite types
neurite_types = ['neurites', 'dendrites', 'axons']

vplots = False

# %% init sholl analysis

# define sholl radius step size
sholl_step_size = 1 # µm
max_sholl_radius = 590 # µm
radius = np.arange(0, max_sholl_radius + sholl_step_size, sholl_step_size)

# create combined sholl profile dataframe
all_sholl_profiles_template = pd.DataFrame(0.,
                                           index = radius,
                                           columns = cell_IDs)
all_sholl_profiles_template_nans = pd.DataFrame(np.nan,
                                                index = radius,
                                                columns = cell_IDs)

# rename index column in template
all_sholl_profiles_template.index.name = 'Radius'

# initialize dictionary for combined sholl profiles
all_sholl_profiles = {'neurites'  : all_sholl_profiles_template.copy(),
                      'dendrites' : all_sholl_profiles_template.copy(),
                      'axons'     : all_sholl_profiles_template.copy()}

# create dataframe for average profiles
avgs = ['mean', 'median', 'std']
avg_columns = [f'{avg}-sholl_profile-{ntype}-{region}' for ntype in neurite_types for region in ['all', 'BAOT', 'MeA'] for avg in avgs]
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

print('loading ...')

for cell_ID in tqdm(cell_IDs):

    # get neurites sholl profile for cell
    sholl_profile_neurites = pd.read_csv(join(cellmorph_sholl_dir, 'sholl_tables_neurites', filenames[cell_ID]), 
                                         index_col = 'Radius',
                                         usecols = ['Radius', 'Inters.'])
    
    # rename column
    sholl_profile_neurites.rename(columns = {'Inters.' : cell_ID},
                                  inplace = True)
    
    # get axon sholl profile for cell
    try:
        sholl_profile_axons = pd.read_csv(join(cellmorph_sholl_dir, 'sholl_tables_axons', filenames[cell_ID][:-12] + ' [axon]' + filenames[cell_ID][-12:]), 
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
        sholl_profile_axons = all_sholl_profiles_template_nans.copy()[[cell_ID]]
        
        # set sholl profile of dendrites as profile of neurites
        sholl_profile_dendrites = sholl_profile_neurites
        
    # write dendrites profile to csv
    sholl_profile_dendrites.rename(columns = {cell_ID : 'Inters.'}).to_csv(join(cellmorph_sholl_dir, 'sholl_tables_dendrites', f'{filenames[cell_ID][:-12]} [dendrites]_profile.csv'),
                                                                           index_label = 'Radius')
    
    # write profiles to dictionary
    all_sholl_profiles['neurites'].loc[sholl_profile_neurites.index, cell_ID] = sholl_profile_neurites[cell_ID]
    all_sholl_profiles['dendrites'].loc[sholl_profile_dendrites.index, cell_ID] = sholl_profile_dendrites[cell_ID]
    all_sholl_profiles['axons'].loc[sholl_profile_axons.index, cell_ID] = sholl_profile_axons[cell_ID]

# save collection of profiles
for ntype in neurite_types:
    all_sholl_profiles[ntype].to_excel(join(cellmorph_analysis_dir, 'sholl-combined_profiles', f'sholl_profiles-{ntype}.xlsx'),
                                        index_label = 'Radius')
    
    
# %% get all cell_IDs

# get cell_IDs
cell_IDs_dict = {'all': cell_IDs,
                 'MeA': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'MeA'].index.to_list() if cell_ID in cell_IDs],
                 'BAOT': [cell_ID for cell_ID in MetaData[MetaData['Region'] == 'BAOT'].index.to_list() if cell_ID in cell_IDs]}
    
# %% calc average profiles

for region in ['all', 'BAOT', 'MeA']:
    
    # get specific cell_IDs:
    region_cellIDs = cell_IDs_dict[region] 
    
    for ntype in neurite_types:
        # mean
        avg_sholl_profiles[f'mean-sholl_profile-{ntype}-{region}'] = all_sholl_profiles[ntype].loc[:, region_cellIDs].mean(axis = 1)
        
        # median
        avg_sholl_profiles[f'median-sholl_profile-{ntype}-{region}'] = all_sholl_profiles[ntype].loc[:, region_cellIDs].median(axis = 1)
        
        # std
        avg_sholl_profiles[f'std-sholl_profile-{ntype}-{region}'] = all_sholl_profiles[ntype].loc[:, region_cellIDs].std(axis = 1)

# save average profiles
avg_sholl_profiles.to_excel(join(cellmorph_analysis_dir, 'sholl-combined_profiles', 'sholl_profiles-avg.xlsx'),
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
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_metrics_dir
sholl_metrics.to_excel(join(cellmorph_metrics_dir, 'sholl_metrics.xlsx'),
                        index_label = 'cell_ID')


# %% create plots

# init plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *

if vplots: 
    
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
                    color = neurite_color_dict['all'][ntype],
                    lw = 1,
                    label = ntype)
        
        # set legend
        ax.legend(title = 'neurite type', 
                  frameon = False,
                  loc = 'upper right')
        
        # y axis
        ydict = {'ax_min' : 0,
                 'ax_max' : 40,
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
        save_figures(fig, 
                     f'{cell_ID}-sholl_profiles', 
                     join(cellmorph_analysis_dir, 'plots-sholl_profile'), 
                     darkmode_bool, 
                     figure_format='png')
        
        # display figure
        plt.show()