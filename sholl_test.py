# -*- coding: utf-8 -*-
"""
Created on Mon May 12 10:06:28 2025

@author: nesseler
"""


# initialize needed packages
from functions.initialize_packages import *

# get directories
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_sholl_dir, cellmorph_analysis_dir, cellmorph_metrics_dir, cellmorph_coordinates_dir

# get files in traces folder
from functions.functions_import import get_onlyfiles_list
onlyfiles = get_onlyfiles_list(cellmorph_sholl_dir + '/sholl_tables_neurites')

# create filenames dict
filenames = {'E-' + filename[-15:-12] : filename for filename in onlyfiles}

# get cell_IDs from list of all coordinates filenames
cell_IDs = list(filenames.keys())
cell_IDs = ['E-162']
# cell_IDs = cell_IDs[:-20] # for leonie NWG poster

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
    
    sholl_profile_neurites.dropna(inplace = True)
    
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


    
    



# %% import cell coordinates

from cellmorphology.cellmorph_functions.cellmorph_parameters import field_of_view

# all coordinates dict for one cell
cell_coordinates = dict.fromkeys(['all_coor', 'end_coor', 'soma_coor', 'path_IDs'])

# import function that cleans up imported coordinates table
from cellmorphology.cellmorph_functions.cellmorph_functions import clean_OnPath_column_to_path_ID_n_label
 
### coordinates
# all coordinates
all_coordinates_path = join(cellmorph_coordinates_dir, f'{cell_ID}-all_coordinates.csv')
cell_allcoordinates = pd.read_csv(all_coordinates_path)
clean_OnPath_column_to_path_ID_n_label(cell_allcoordinates)

# primary last coordinates
primary_coordinates_path = join(cellmorph_coordinates_dir, f'{cell_ID}-primary_last_coordinates.csv')
cell_primarycoordinates = pd.read_csv(primary_coordinates_path)
clean_OnPath_column_to_path_ID_n_label(cell_primarycoordinates) 

# end / last / terminal coordinates
last_coordinates_path = join(cellmorph_coordinates_dir, f'{cell_ID}-terminal_last_coordinates.csv')
cell_endcoordinates = pd.read_csv(last_coordinates_path)
clean_OnPath_column_to_path_ID_n_label(cell_endcoordinates) 
    
# get soma coordinates
soma_coordinates = cell_allcoordinates[cell_allcoordinates['path_ID'] == 1]
 
# check if x coordinates need to be flipped
if MetaData.at[cell_ID, 'to_be_x_flipped']:
    cell_allcoordinates.loc[:, 'X']     = field_of_view - cell_allcoordinates['X']
    cell_primarycoordinates.loc[:, 'X'] = field_of_view - cell_primarycoordinates['X']
    cell_endcoordinates.loc[:, 'X']     = field_of_view - cell_endcoordinates['X']
    soma_coordinates.loc[:, 'X']        = field_of_view - soma_coordinates['X']

# get all path_IDs
path_IDs = cell_allcoordinates['path_ID'].drop_duplicates().astype(int).to_list()

# set dict for all coordinates
cell_coordinates = {'all_coor' : cell_allcoordinates,
                    'pri_coor' : cell_primarycoordinates,
                    'end_coor' : cell_endcoordinates,
                    'soma_coor': soma_coordinates,
                    'path_IDs' : path_IDs}
   

# load plotting function    
from cellmorphology.cellmorph_analysis.plot_cellcoordinates_analysis import plot_cellcoordinates

# plot 
plot_cellcoordinates(cell_ID = cell_ID, 
                     cell_coordinates = cell_coordinates)


# %% pixel size

pixel_size = 0.2885

# import matplotlib.pyplot as plt
import shapely
# import math
# import numpy as np

from shapely import LineString
from shapely.wkt import loads

line = LineString([(0, 0), (2, 2)])

shapely.intersection(line, LineString([(1, 1), (3, 3)]))
# <LINESTRING (1 1, 2 2)>

box1 = shapely.box(0, 0, 2, 2)

box2 = shapely.box(1, 1, 3, 3)

shapely.intersection(box1, box2).normalize()
# <POLYGON ((1 1, 1 2, 2 2, 2 1, 1 1))>

box1 = shapely.box(0.1, 0.2, 2.1, 2.1)

shapely.intersection(box1, box2, grid_size=1)
# <POLYGON ((2 2, 2 1, 1 1, 1 2, 2 2))>

r = 1
n = 64
t = np.linspace(0, 2 * np.pi, n + 1)

x_circle = r * np.cos(t) + 1
y_circle = r * np.sin(t) + 1

# create the linestring of circle's perimeter
wktcode1 = "LINESTRING ("
for i,(x,y) in enumerate(zip(x_circle, y_circle)):
    if i!=len(x_circle)-1:
        wktcode1 += str(x)+" "+str(y)+", "
    else:
        wktcode1 += str(x)+" "+str(y)+")"
    #print(wktcode1)

circle_perim = loads(wktcode1)  #a geometry object

# create another geometry, for the line
wktcode2 = "LINESTRING ("
xs = range(4)
ys = np.array(range(4))*0.42
for i,(x,y) in enumerate(zip(xs, ys)):
    if i!=len(range(4))-1:
        wktcode2 += str(x)+" "+str(y)+", "
    else:
        wktcode2 += str(x)+" "+str(y)+")"
    #print(wktcode2)
    pass

line_str = loads(wktcode2)    #a geometry object

# check if they are intersect
# ixx = circle_perim.intersection(box1)


# check if they are intersect
ixx = circle_perim.intersection(line_str)

intersec = line_str.intersection(circle_perim)
print(intersec.wkt)

# visualization of the intersection
# plot circle
plt.scatter(x_circle, y_circle, s = 1)
    
# x1,y1 = 
plt.plot(*box1.exterior.xy)
plt.plot(*line_str.xy)
# plt.scatter(*intersec.xy)

# # if ixx:
# #     # plot intersection points
for p in list(intersec.geoms):
    plt.scatter(*p.xy, s = 10, marker = 'x', lw = 1, c = 'r')


plt.gca().set_aspect(1)
plt.show()

# %%

# import cv2

plt.scatter(cell_allcoordinates['X'],
            cell_allcoordinates['Y'],
            s = 1,
            marker = 's')

plt.scatter(soma_coordinates['X'],
            soma_coordinates['Y'],
            s = 1,
            c = 'r')

# plt.ylim([field_of_view, 0])


plt.gca().set_aspect(1)

# test site for sholl analysis function


r = 3

stepSize = pixel_size/2

#Generated vertices
positions = []

t = 0
while t < 2 * math.pi:
    positions.append((r * math.cos(t) + soma_coordinates['X'].values[0], r * math.sin(t) + soma_coordinates['Y'].values[0]))
    t += stepSize



    

plt.scatter(*zip(*positions),
            s = 1,
            c = 'r')




plt.xlim([235, 255])
plt.ylim([310, 330])

plt.show()


# %% create plots

# init plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *

    
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