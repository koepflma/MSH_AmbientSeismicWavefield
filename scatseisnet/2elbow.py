"""
This file calculated the distance between the time windows and their cluster center.
The output can be used to create the elbow curve.
The input is provided by 1scatseisnet.py
"""
import sys
import numpy as np
import sklearn
from sklearn.cluster import KMeans

# parameter ========================================================================================================

n_clusters = 20 # define the maximum number of cluster to test

# path to save figures
dirpath = "../output/scatseisnet/"

# Load features
features = np.load(dirpath+"data/independent_components.npz", allow_pickle=True)['features']
times = np.load(dirpath+"data/independent_components.npz", allow_pickle=True)['times']
print('features shape: {}\ntimes shape: {}'.format(features.shape,times.shape))

# clustering for different numbers of clusters =======================================================================

sse = np.zeros(n_clusters) # initialize array

for N_CLUSTERS in range(1,n_clusters+1):

    sys.stdout.write('\r{} of {}'.format(N_CLUSTERS+1, n_clusters+1))
    sys.stdout.flush()
    
    # Perform clustering
    model_cluster = KMeans(n_clusters=N_CLUSTERS, n_init='auto', random_state=4)
    model_cluster.fit(features)
    
    #Sum of squared distances of samples to their closest cluster center, weighted by the sample weights if provided
    sse[N_CLUSTERS-1] = model_cluster.inertia_ 

np.savez(dirpath+"data/sse.npz",sse=sse)
print('\n***DONE***')
