# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 10:21:37 2025

@author: nesseler
"""

### cell morphology ###

# directory for cell morphology measures (SNT traces/Measurements)
cell_morph_parent = 'Z:/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA'
# cell_morph_parent = '/Users/moritznesseler/n2021_MOS_AOS_Integration/cellmorphology_BAOT_MeA'


# OLD
# cell_morph_descrip_dir = cell_morph_parent + '/SNT_traces/'

# cell_morph_measures_dir = cell_morph_parent + '/SNT_traces/Measurements'

# cell_morph_figures_dir = cell_morph_parent + '/cellmorph_figures'

# cell_morph_traces_coordinates_dir = cell_morph_parent + '/SNT_traces/traces_coordinates'

# cell_morph_plots_dir = cell_morph_parent + '/Plots'

# cell_morph_traces_sholl_dir = cell_morph_parent + '/SNT_traces/traces_sholl_tables'


# NEW
cellmorph_analysis_dir = cell_morph_parent + '/cellmorph_analysis'
cellmorph_figures_dir = cell_morph_parent + '/cellmorph_figures'
cellmorph_reconstructions_dir = cell_morph_parent + '/cellmorph_reconstructions'

cellmorph_traces_dir = cellmorph_reconstructions_dir + '/Traces'
cellmorph_coordinates_dir = cellmorph_reconstructions_dir + '/cell_coordinates'

cellmorph_metrics_dir = cellmorph_analysis_dir + '/metrics'