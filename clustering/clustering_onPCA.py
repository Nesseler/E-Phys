# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 17:10:43 2025

@author: nesseler
"""

from functions.initialize_packages import *

# specific package
from sklearn.decomposition import PCA

# load celldescriptors
from parameters.directories_win import clustering_dir
celldescriptors = pd.read_excel(join(clustering_dir, 'ePhys_celldescriptors.xlsx'), index_col = 'cell_ID')
celldescriptors_clustered = pd.read_excel(join(clustering_dir, 'ePhys_celldescriptors-clustered.xlsx'), index_col = 'cell_ID')


# get cell_IDs
cell_IDs = celldescriptors.index.to_list()

# load Metadata
from functions.functions_import import get_MetaData
MetaData = get_MetaData(cell_IDs)


# z-score cellmorph matrix
celldescriptors_zscored = (celldescriptors - celldescriptors.mean()) / celldescriptors.std()


# %% principal component analysis

# set up PCA
celldescriptors_PCA = PCA().fit(celldescriptors_zscored)

# perform transfrom
celldescriptors_PCA_components = celldescriptors_PCA.transform(celldescriptors_zscored)

# write to dataframe
prinicpal_components = pd.DataFrame(celldescriptors_PCA_components,
                                    index = cell_IDs,
                                    columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])])

# get eigenvectors
eigenvectors = celldescriptors_PCA.components_

celldescriptors_PCA_eigenvectors = pd.DataFrame(eigenvectors,
                                                columns = celldescriptors.columns,
                                                index = [f'PC{i+1}' for i in range(celldescriptors.shape[1])])

# get the eigenvalues
celldescriptors_PCA_eigenvalues = celldescriptors_PCA.explained_variance_

# get components loading
loadings = eigenvectors.T * np.sqrt(celldescriptors_PCA_eigenvalues)

celldescriptors_PCA_loadings = pd.DataFrame(loadings,
                                            columns = [f'PC{i+1}' for i in range(celldescriptors.shape[1])],
                                            index = celldescriptors.columns)

# combine to one dataframe 
celldescriptors_PCA_metrics = pd.DataFrame(index = [f'PC{i+1}' for i in range(celldescriptors_zscored.shape[1])])

# get eigenvalues and calcuate the explained variance ratio 
celldescriptors_PCA_metrics.loc[:, "eigenvalue"] = celldescriptors_PCA_eigenvalues
celldescriptors_PCA_metrics.loc[:, "explained_variance"] = celldescriptors_PCA_eigenvalues / celldescriptors_zscored.var().sum()

# calc cumulative sum of explained variance ratio
celldescriptors_PCA_metrics.loc[:, "cumsum_explained_variance"] = celldescriptors_PCA_metrics.loc[:, "explained_variance"].cumsum()




colors = [(0.5529411764705883, 0.8274509803921568, 0.7803921568627451),
         (0.996078431372549, 1.0, 0.7019607843137254),
         (0.7490196078431373, 0.7333333333333333, 0.8509803921568627),
         (0.9803921568627451, 0.5058823529411764, 0.4549019607843137),
         (0.5058823529411764, 0.6941176470588235, 0.8235294117647058),
         (0.9921568627450981, 0.7058823529411765, 0.3843137254901961),
         (0.7019607843137254, 0.8705882352941177, 0.4117647058823529),
         (0.7372549019607844, 0.5098039215686274, 0.7411764705882353),
         (0.8, 0.9215686274509803, 0.7686274509803922),
         (1.0, 0.9294117647058824, 0.43529411764705883)]

c_colors = [colors[i] for i in celldescriptors_clustered.loc[cell_IDs, 'hierarchical_cluster'].to_list()]


# get number of clusters
n_clusters = celldescriptors_clustered.loc[cell_IDs, 'hierarchical_cluster'].nunique()




# %% initialize plotting

from functions.initialize_plotting import *


# %%

from sklearn.cluster import KMeans

data = celldescriptors_PCA_components

# 
model = KMeans(n_clusters = 5, init = "k-means++", random_state=7) #5 #7
# centers = np.array(celldescriptors_PCA.cluster_centers_)
label = model.fit_predict(data)
plt.figure(figsize=(10,10))
uniq = np.unique(label)
for i in uniq:
   plt.scatter(data[label == i , 0] , data[label == i , 1] , label = i)
   # plt.scatter(centers[:,0], centers[:,1], marker="x", color='k')
#This is done to find the centroid for each clusters.
plt.legend()
plt.show()

# %%



inertia = []
for k in range(1, 9):
    model = KMeans(n_clusters=k, random_state=1).fit(celldescriptors)
    # label = model.fit(data)
    inertia.append(np.sqrt(model.inertia_))

    # print(kmeans.inertia_)

plt.plot(range(1, 9), inertia, marker='s');
plt.xlabel('$k$')
plt.ylabel('$J(C_k)$');
plt.show()

# %%

from sklearn.cluster import AffinityPropagation


# Configuration options
# num_samples_total = 50
# cluster_centers = [(20,20), (4,4)]
# num_classes = len(cluster_centers)
# 
# Generate data
# X, targets = make_blobs(n_samples = num_samples_total, centers = cluster_centers, n_features = num_classes, center_box=(0, 1), cluster_std = 1)

# np.save('./clusters.npy', X)
X = celldescriptors_PCA_components

# Fit AFfinity Propagation with Scikit
afprop = AffinityPropagation(max_iter=250)
afprop.fit(X)
cluster_centers_indices = afprop.cluster_centers_indices_
n_clusters_ = len(cluster_centers_indices)

# Predict the cluster for all the samples
P = afprop.predict(X)

# Generate scatter plot for training data
colors = [(0.5529411764705883, 0.8274509803921568, 0.7803921568627451),
         (0.996078431372549, 1.0, 0.7019607843137254),
         (0.7490196078431373, 0.7333333333333333, 0.8509803921568627),
         (0.9803921568627451, 0.5058823529411764, 0.4549019607843137),
         (0.5058823529411764, 0.6941176470588235, 0.8235294117647058),
         (0.9921568627450981, 0.7058823529411765, 0.3843137254901961),
         (0.7019607843137254, 0.8705882352941177, 0.4117647058823529),
         (0.7372549019607844, 0.5098039215686274, 0.7411764705882353),
         (0.8, 0.9215686274509803, 0.7686274509803922),
         (1.0, 0.9294117647058824, 0.43529411764705883)]

# plt.scatter(X[:,0], X[:,1], c=colors, marker="o", picker=True)

for ci in range(len(P)):
    plt.scatter(X[ci,0], X[ci,1], 
                c=colors[P[ci]], 
                marker="o")

# for i in uniq:
#    plt.scatter(data[label == i , 0] , data[label == i , 1] , label = i)

plt.title(f'Estimated number of clusters = {n_clusters_}')
# plt.xlabel('Temperature yesterday')
# plt.ylabel('Temperature today')
plt.show()



# %%
