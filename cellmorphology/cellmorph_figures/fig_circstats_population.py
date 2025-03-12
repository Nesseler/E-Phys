# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 15:32:36 2025

@author: nesseler
"""

# initialize needed packages
from functions.initialize_packages import *

# get directory
from cellmorphology.cellmorph_functions.cellmorph_dir import cellmorph_metrics_dir

# re-load orientations dataframe to include all cells
circ_stats = pd.read_excel(join(cellmorph_metrics_dir, 'circ_stats.xlsx'),
                           index_col = 'cell_ID')

# get cell_IDs
cell_IDs = circ_stats.index.to_list()

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)

# set neurite types
neurite_types = ['neurites', 'dendrites', 'axons']

# initialize dataframe for populations
population_circstats = pd.DataFrame(columns = ['circmean', 'circstd', 'circvar'],
                                    index = [f'{ntype}-{region}' for region in ['MeA', 'BAOT'] for ntype in neurite_types])

for region in ['BAOT', 'MeA']:
    
    # get region cell_IDs
    region_cellIDs = MetaData[MetaData['Region'] == region].index.to_list()
    
    # limit to reconstructed cells in region
    region_cellIDs = [cell_ID for cell_ID in cell_IDs if cell_ID in region_cellIDs]
    
    # get region circstats
    region_circstats = circ_stats.loc[region_cellIDs, :]

    for ntype in neurite_types:

        # get all terminal angles
        t_angles_pertype = list()
        
        for cell_ID, row in region_circstats.iterrows(): 
            for angle in row[f'terminal_angles_rad-{ntype}'].strip(' []').split(','):
                if len(angle) > 0:
                    t_angles_pertype.append(float(angle))
        
        # write to dataframe
        population_circstats.at[f'{ntype}-{region}', 'circmean'] = sc.stats.circmean(t_angles_pertype)
        population_circstats.at[f'{ntype}-{region}', 'circstd'] = sc.stats.circstd(t_angles_pertype)
        population_circstats.at[f'{ntype}-{region}', 'circvar'] = sc.stats.circvar(t_angles_pertype)
        
        # print(sc.stats.circmean(t_angles_pertype), region_circstats[f'circmean_rad-{ntype}'].mean())
        
        
        
# %% figure

# init plotting
from cellmorphology.cellmorph_functions.cellmorph_init_plotting import *



fig, axs = plt.subplots(nrows = 3,
                        ncols = 4,
                        subplot_kw = {'projection': 'polar'},
                        figsize = get_figure_size(width = 200, height = 150),
                        dpi = 300, 
                        layout = 'constrained',
                        sharex = True,
                        sharey = 'col')

# flatten axis
axs = axs.flatten()

for n_idx, ntype in enumerate(neurite_types):

    for r_idx, region in enumerate(['BAOT', 'MeA']):
        
        # get region cell_IDs
        region_cellIDs = MetaData[MetaData['Region'] == region].index.to_list()
        
        # limit to reconstructed cells in region
        region_cellIDs = [cell_ID for cell_ID in cell_IDs if cell_ID in region_cellIDs]
        
        # get region circstats
        region_circstats = circ_stats.loc[region_cellIDs, :]
    
        ### points
    
        # set axis
        ax_p = axs[(n_idx*4) + (r_idx*2)]
        ax_l = axs[(n_idx*4) + (r_idx*2) + 1]
        
        # set axis title
        ax_p.set_title(f'{region} {ntype}',
                       fontsize=9, 
                       loc='left',
                       x = -0.1)
        
        # get all terminal angles
        t_angles_pertype = list() 

        for c_idx, (cell_ID, row) in enumerate(region_circstats.iterrows()): 
            
            if not np.isnan(region_circstats.at[cell_ID, f'circmean_rad-{ntype}']):
            
                # get all terminal angles
                t_angles_percell = list() 
                
                for angle in row[f'terminal_angles_rad-{ntype}'].strip(' []').split(','):
                    t_angles_pertype.append(float(angle))
                    t_angles_percell.append(float(angle))
            
                ax_p.scatter(x = t_angles_percell,
                             y = [c_idx+15] * len(t_angles_percell),
                             color = neurite_color_dict[region][ntype],
                             s = 0.5)
    
                # get circ stats
                circmean_rad = region_circstats.at[cell_ID, f'circmean_rad-{ntype}']
                circmean_deg = region_circstats.at[cell_ID, f'circmean_deg-{ntype}']
                circstd_rad = region_circstats.at[cell_ID, f'circstd_rad-{ntype}']
            
                # transform marker to rotate marker
                t = mtl.markers.MarkerStyle(marker='_')
                t._transform = t.get_transform().rotate_deg(circmean_deg)
                
                # plot circular mean and std
                ax_l.errorbar(x = circmean_rad,
                              y = c_idx +15, 
                              xerr = circstd_rad,
                              marker = t,
                              markeredgewidth = 1,
                              markersize = 3,
                              linewidth = 0.4, 
                              color = neurite_color_dict[region][ntype],
                              label = 'circ mean+std')
                
        # y axis
        ymax = np.ceil((region_circstats.shape[0]+15)/5)*5 + 5
        [ax.set_ylim([0, ymax+7]) for ax in [ax_p, ax_l]]
                
        # plot overall means
        points_mean = sc.stats.circmean(t_angles_pertype)
        points_std = sc.stats.circstd(t_angles_pertype)
        
        # transform marker to rotate marker
        t = mtl.markers.MarkerStyle(marker='_')
        t._transform = t.get_transform().rotate_deg(np.rad2deg(points_mean))
        
        # plot circular mean and std
        ax_p.errorbar(x = points_mean,
                      y = [ymax], 
                      xerr = points_std,
                      marker = t,
                      markeredgewidth = 1,
                      markersize = 3,
                      linewidth = 0.5, 
                      color = 'k',
                      label = 'circ mean+std')
        
        # plot overall means
        lines_mean = sc.stats.circmean(region_circstats.loc[:, f'circmean_rad-{ntype}'].dropna())
        lines_std = sc.stats.circmean(region_circstats.loc[:, f'circstd_rad-{ntype}'].dropna())
        # lines_std = np.mean(region_circstats.loc[:, f'circstd_rad-{ntype}'].dropna())
        
        # transform marker to rotate marker
        t = mtl.markers.MarkerStyle(marker='_')
        t._transform = t.get_transform().rotate_deg(np.rad2deg(lines_mean))
        
        # plot circular mean and std
        ax_l.errorbar(x = lines_mean,
                      y = [ymax], 
                      xerr = lines_std,
                      marker = t,
                      markeredgewidth = 1,
                      markersize = 3,
                      linewidth = 0.5, 
                      color = 'k',
                      label = 'circ mean+std')
                


for ax in axs:
    
    # x axis
    ax.set_xticks(np.arange(0, np.pi*2, np.pi / 4))
    ax.set_xticklabels(['p', 'pd', 'd', 'ad', 'a', 'av', 'v', 'pv'])

    # y axis
    ax.set_yticks([])

    # set grid
    ax.grid(True, alpha = 0.1)
    
    # set grid behind plot
    ax.set_axisbelow(True)
    

    
plt.show()