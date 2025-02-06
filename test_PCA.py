# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 17:05:06 2025

@author: nesseler

ressource
https://www.geo.fu-berlin.de/en/v/soga-py/Advanced-statistics/Multivariate-Approaches/Principal-Component-Analysis/PCA-the-basics/Choose-Principal-Components/index.html

"""

from functions.initialize_packages import *

# specific package
from sklearn.decomposition import PCA

# load celldescriptors
from parameters.directories_win import clustering_dir
celldescriptors = pd.read_excel(join(clustering_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')

# get cell_IDs
cell_IDs = celldescriptors.index.to_list()

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)


# %% z-score celldescriptors

# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()


# %% initialize plotting

from functions.initialize_plotting import *


# %% principal component analysis

celldescriptors_PCA = PCA().fit(celldescriptors)

celldescriptors_PCA_components = celldescriptors_PCA.transform(celldescriptors)

celldescriptors_PCA_explained_variance = celldescriptors_PCA.explained_variance_

# # get the eigenvectors and eigenvalues
#     # index: celldescriptors
#     # columns: principal components
celldescriptors_PCA_eigen = pd.DataFrame(celldescriptors_PCA.components_.T,
                                         columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])],
                                         index = celldescriptors.columns)

# celldescriptors_PCA_eigen["eigenvalue"] = celldescriptors_PCA.explained_variance_

prinicpal_components = pd.DataFrame(celldescriptors_PCA_components,
                                    index = cell_IDs,
                                    columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])])




sbn.scatterplot(data = prinicpal_components,
                x = 'PC1',
                y = 'PC2',
                hue = MetaData.loc[cell_IDs, 'Region'],
                palette = region_colors)


# plot the variables as vectors
plt.quiver(
    np.zeros(celldescriptors_PCA_eigen.shape[0]),
    np.zeros(celldescriptors_PCA_eigen.shape[0]),
    celldescriptors_PCA_eigen["PC1"],
    celldescriptors_PCA_eigen["PC2"],
    # celldescriptors_PCA_eigen["eigenvalue"],
    angles="xy",
    # scale_units="xy",
    # scale=0.2,
    color = 'w'
)


# Plot annotations
for i in range(celldescriptors_PCA_eigen.shape[0]):
    plt.text(
        celldescriptors_PCA_eigen["PC1"][i]*200,
        celldescriptors_PCA_eigen["PC2"][i]*200,
        celldescriptors_PCA_eigen.index[i],
        color="red",
    )