# Scatering Network for Seismology

import os
import pickle
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import obspy
from obspy.clients.fdsn.client import Client 
import time
import sys
import io
from scatseisnet import ScatteringNetwork

s1  = 'UW-EDM-EHZ'
s2  = 'UW-SHW-EHZ'
s3  = 'UW-HSR-EHZ'

list_stations = [s1,s2,s3] # make a list of all stations
  
# path to save figures
dirpath_save = '../output/data/'
# bild path to save figures if it does not exist yet
os.makedirs(dirpath_save, exist_ok=True)

print('Package imported and folder created.')

stime = time.time() # start time

# load data ========================================================================================================
client = Client("IRIS")

stream = obspy.Stream()
list_netstacha = []

for xx in list_stations:
    # Collect waveforms from the datacenter
    
    netstacha = xx.split('-')
    net = netstacha[0] # network
    sta = netstacha[1] # station
    cha = netstacha[2] # channel
    
    try:
        st = client.get_waveforms(
            network=net,
            station=sta,
            location="*",
            channel=cha,
            starttime=obspy.UTCDateTime("2004-09-01T00:00:00"), #"2004-09-01T00:00"
            endtime=obspy.UTCDateTime("2004-09-03T00:00:00.1"), #"2004-11-01T00:00.1"
        )
#         st.trim(obspy.UTCDateTime("1980-03-31T00:00"),obspy.UTCDateTime("1980-05-19T00:00:00.1"),pad=True, fill_value=0)

        st.detrend("linear")
        st.taper(0.05, type='hann')
        st.resample(50)
        st.merge(fill_value=0)

        st.filter(type="highpass", freq=1.0)
        stream += st
        list_netstacha.append(xx)
        print('{}'.format(netstacha))
    except:
        print('no data for {}'.format(netstacha))
        pass
    
if len(stream) == 0:
    print('empty stream')
    sys.exit()
if len(stream) != 0:
    stream.write(dirpath_save+"scattering_stream.mseed", format="MSEED")
    print(stream)

print('Loading data tooks {} minutes.'.format(round((time.time()-stime)/60,3)))

# create layers ======================================================================================================

segment_duration_seconds = 60.0
sampling_rate_hertz = stream[0].stats.sampling_rate
samples_per_segment = int(segment_duration_seconds * sampling_rate_hertz)
bank_keyword_arguments = (
    {"octaves": 4, "resolution": 4, "quality": 1},
    {"octaves": 5, "resolution": 2, "quality": 2},
)

# create scatnet ======================================================================================================

network = ScatteringNetwork(
    *bank_keyword_arguments,
    bins=samples_per_segment,
    sampling_rate=sampling_rate_hertz,
)

print(network)

# Save the scattering network with Pickle
filepath_save = os.path.join(dirpath_save, "scattering_network.pickle")
with open(filepath_save, "wb") as file_save:
    pickle.dump(network, file_save, protocol=pickle.HIGHEST_PROTOCOL)

# preparation for network transform ===================================================================================

# Extract segment length (from any layer)
segment_duration = network.bins / network.sampling_rate
overlap = 0.5

# Gather list for timestamps and segments
timestamps = list()
segments = list()

# Collect data and timestamps
for traces in stream.slide(segment_duration, segment_duration * overlap):
    timestamps.append(mdates.num2date(traces[0].times(type="matplotlib")[0]))
    segments.append(np.array([trace.data[:-1] for trace in traces]))
    
# scattering transform ===============================================================================================
st = time.time() # start time
scattering_coefficients = network.transform(segments, reduce_type=np.mean)

# Save the features in a pickle file
np.savez(
    dirpath_save+"scattering_coefficients.npz",
    order_1=scattering_coefficients[0],
    order_2=scattering_coefficients[1],
    times=timestamps,
#     list_netstacha=list_stations
    list_netstacha=list_netstacha
)
print('Scattering transform and saving took {} hours.'.format(round((time.time()-st)/3600,3)))
print('Script run for {} hours.'.format(round((time.time()-stime)/3600,3)))