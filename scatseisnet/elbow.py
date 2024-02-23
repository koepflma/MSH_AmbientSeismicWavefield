import os
import io
import sys
import pickle

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

import numpy as np
import pandas as pd
import obspy
from obspy.clients.fdsn.client import Client
client = Client('IRIS')
import sklearn
from sklearn.decomposition import FastICA
from sklearn.cluster import KMeans

from scatseisnet import ScatteringNetwork

# model (octave, resolution,quality)
str_network = 'o4_r4_q1_o5_r2_q2'
# duration
duration = '1month'
month = 'apr'
# seismic station
str_station = 'multistation'

# path to load and save averything exept of figures
dirpath_load = '/data/wsd03/data_manuela/MtStHelens/scatseisnet/{}/{}/{}/{}'.format(str_network,duration,str_station,month)

# path to save figures
dirpath_save = '/home/koepflma/ESS590D_Manuela_Koepfli/plots/scatseisnet/{}/{}/{}/{}'.format(str_network,duration,str_station,month)
os.makedirs(dirpath_save, exist_ok=True) # bild path to save figures if it does not exist yet

# Load features
features = np.load(dirpath_load+"/independent_components.npz", allow_pickle=True)['features']
times = np.load(dirpath_load+"/independent_components.npz", allow_pickle=True)['times']
print('features shape: {}\ntimes shape: {}'.format(features.shape,times.shape))
# Load the dimension reduction model
model_latent = pd.read_pickle(dirpath_load+"/dimension_model.pickle")
model_latent



#-----------------------------------------------------------------------------------------------------------------------------
n_clusters = 50
E = np.zeros(n_clusters)
sse = np.zeros(n_clusters)

def compute_elbow(distances):
    """
    calculates distance between features to cluster mean --> not needed
    distances: list include arrays of distance (cummulative distances of feature sample to cluster mean, per time window) by cluster [np.array(1,2,1,2),np.array(3,4,2,3)] --> 2 clusters, and 4 time windows
    """
    E = 0
    for i in range(0, len(distances)): # loops over of clusters (n_cluster = 5, first for cluster 1 than 2...5)
#         distance = compute_distance(data[clusters == i, :], centers[i, :].reshape(1, -1), 1)
        E = E + np .mean(np.square(distance))
#         E += np.mean(distances[i-1])

    return E

for N_CLUSTERS in range(1,n_clusters+1):

    sys.stdout.write('\r{} of {}'.format(N_CLUSTERS, n_clusters+1))
    sys.stdout.flush()
    
    # Perform clustering
    model_cluster = KMeans(n_clusters=N_CLUSTERS, n_init='auto', random_state=4) #n_init='auto' or 1
    model_cluster.fit(features)
    
    sse[N_CLUSTERS-1] = model_cluster.inertia_ #Sum of squared distances of samples to their closest cluster center, weighted by the sample weights if provided

    # Predict cluster for each sample
    predictions = model_cluster.predict(features)

    # Save the prediction
    np.savez(
        dirpath_load+"/predictions_cl{}_elbow50.npz".format(N_CLUSTERS),
        predictions=predictions,
    )
    
    # Extract distance
    distances = list()
    for cluster in np.unique(predictions):
#         print('Cluster {}'.format(cluster))

        # Calculate the distance of each sample to the cluster mean
        mean = np.mean(features[predictions == cluster], axis=0)
        distance = np.linalg.norm(features[predictions == cluster] - mean, axis=1)
#             closest = times[predictions == cluster][distance.argsort()[:N_WAVEFORMS]]
        distances.append(distance)
    
#     Create a dictionary with named arrays
    distances_dict = {f'array_{i}': array for i, array in enumerate(distances)}
    np.savez(dirpath_load+"/distances50.npz",distances=distances_dict)
    
    E[N_CLUSTERS - 1] = compute_elbow(distances) 
    
#     if N_CLUSTERS != 1:

#         # plot time series-----------------------------------------------------------------------------------------------------
#         SMOOTH_KERNEL = 50

#         # Convert predictions to one-hot encoding
#         one_hot = np.zeros((len(times), N_CLUSTERS + 1))
#         one_hot[np.arange(len(times)), predictions] = 1

#         # Plot the results
#         fig, ax = plt.subplots(N_CLUSTERS,1,sharex=True,figsize=(6.4*2, 4.8*1.5))

#         # Plot each cluster as a separate line
#         for i in range(N_CLUSTERS):

#             # Obtain the detection rate by convolving with a boxcar kernel
#             detection_rate = np.convolve(one_hot[:, i], np.ones(SMOOTH_KERNEL), mode="same") / SMOOTH_KERNEL

#             # Plot the detection rate
#         #     ax[i].plot(times, one_hot[:, i] + i, alpha=0.5)
#             ax[i].plot(times, detection_rate, color="C{}".format(i))
#             ax[i].axvline(datetime.datetime(2004,9,23,9), color='black', linewidth=2, linestyle=':')

#             # Y-Label
#             ax[i].set_yticks([])
#         #     ax[i].set_ylabel("Cluster activity")

#         # X-Label
#         ax[i].set_xlabel("Time")
#         ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))# %H:%M:%S
#         ax[i].xaxis.set_major_locator(mdates.DayLocator(interval=7))

#         plt.savefig(dirpath_save+'/timeseries_subplots_cl{}_elbow50.png'.format(N_CLUSTERS), dpi=300, bbox_inches='tight')


#         centroids = np.abs(model_cluster.cluster_centers_)

#         # Plot the centroids--------------------------------------------------------------------------------------------------
#         fig = plt.figure()
#         ax = plt.axes()

#         # Show the centroids as a heatmap
#         mappable = ax.matshow(centroids.T, cmap="RdPu")

#         # Labels
#         plt.colorbar(mappable).set_label("Amplitude")
#         ax.set_xlabel("Cluster index")
#         ax.set_ylabel("Feature index")

#         # Ticks below
#         ax.xaxis.set_ticks_position("bottom")
#         ax.set_xticks(np.arange(N_CLUSTERS))
#         ax.set_yticks(np.arange(centroids.shape[1]))
#         ax.invert_yaxis()
#         plt.savefig(dirpath_save+'/cluster_feature_cl{}_elbow50.png'.format(N_CLUSTERS), dpi=300, bbox_inches='tight')

#         # plot the number of time windows in each cluster------------------------------------------------------------------------
#         for i, dist in enumerate(distances):
#             plt.bar(i,dist.shape,color='C{}'.format(i))

#         plt.xticks(range(len(distances)), ['Cluster {}'.format(i) for i in range(len(distances))])
#         plt.ylabel('# Time Windows')
#         plt.savefig(dirpath_save+'/barplot_timewindows_cl{}_elbow50.png'.format(N_CLUSTERS), dpi=300, bbox_inches='tight')

#         # plot distances of features for the different clusters-----------------------------------------------------------------
#         fig, ax = plt.subplots(N_CLUSTERS, 1,sharex=True,figsize=(6.4*2,4.8*1.5),constrained_layout=True)
#         for i,dist in enumerate(distances):
#             ax[i].plot(range(len(dist)),dist,color='C{}'.format(i))
#             ax[i].set_xlim(0,len(max(distances, key=len)))
#             ax[i].set_ylim(0,max(max(distances, key=tuple)))
#         fig.supylabel('Feature Distance')
#         fig.supxlabel("# Time Windows")
#         plt.savefig(dirpath_save+'/distance_features_cl{}_elbow50.png'.format(N_CLUSTERS), dpi=300, bbox_inches='tight')

# np.savez(dirpath_load+"/distances50.npz",distances=dist_arr)    
np.savez(dirpath_load+"/elbow50.npz",elbow=E)
np.savez(dirpath_load+"/sse50.npz",sse=sse)
print('\n***DONE***')