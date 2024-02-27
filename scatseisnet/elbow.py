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
duration = '2month'
# month = ''
# seismic station
str_station = 'multistation'

# path to load and save averything exept of figures
dirpath_load = '/data/wsd03/data_manuela/MtStHelens/scatseisnet/{}/{}/{}/'.format(str_network,duration,str_station)


# # path to load and save everything
# dirpath = '../data/scatseisnet/'

# path to save figures
dirpath_save = '../output/data/'
# bild path to save figures if it does not exist yet
os.makedirs(dirpath_save, exist_ok=True)

# Load features
features = np.load(dirpath_load+"/independent_components.npz", allow_pickle=True)['features']
times = np.load(dirpath_load+"/independent_components.npz", allow_pickle=True)['times']
print('features shape: {}\ntimes shape: {}'.format(features.shape,times.shape))

#-----------------------------------------------------------------------------------------------------------------------------
n_clusters = 20
sse = np.zeros(n_clusters)

for N_CLUSTERS in range(1,n_clusters+1):

    sys.stdout.write('\r{} of {}'.format(N_CLUSTERS, n_clusters+1))
    sys.stdout.flush()
    
    # Perform clustering
    model_cluster = KMeans(n_clusters=N_CLUSTERS, n_init='auto', random_state=4)
    model_cluster.fit(features)
    
    sse[N_CLUSTERS-1] = model_cluster.inertia_ #Sum of squared distances of samples to their closest cluster center, weighted by the sample weights if provided

np.savez(dirpath_save+"sse.npz",sse=sse)
print('\n***DONE***')