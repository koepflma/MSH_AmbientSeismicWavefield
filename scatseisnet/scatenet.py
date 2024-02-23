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

# Don't forget to make the adjustments for the choosen network!!!

s1  = 'UW-EDM-EHZ'
s2  = 'UW-SHW-EHZ' # skip 
s3  = 'UW-HSR-EHZ'
s4  = 'UW-SOS-EHZ'
s5  = 'UW-JUN-EHZ'
s6  = 'UW-ELK-EHZ'
s7  = 'UW-TDL-EHZ'
s8  = 'UW-SUG-EHZ'
s9  = 'UW-YEL-EHZ'
s10 = 'UW-FL2-EHZ'
s11 = 'UW-CDF-?H?'

s12 = 'UW-SEP-EHZ'
s13 = 'CC-SEP-EHZ'
# s14 = 'UW-STD-EHZ'
s15 = 'CC-STD-BHZ'

s16 = 'CC-VALT-BHZ'
s17 = 'CC-JRO-BHZ'
s18 = 'CC-HOA-BHZ'
s19 = 'CC-LOO-BHZ'
s20 = 'CC-USFR-BHZ'
s21 = 'CC-NED-EHZ'
s22 = 'CC-REM-BHZ'
s23 = 'CC-SWFL-BHZ'
s24 = 'CC-SFW2-BHZ'
s25 = 'CC-MIDE-EHZ'
s26 = 'CC-MIBL-EHZ'
s27 = 'CC-BLIS-EHZ'
s28 = 'CC-RAFT-EHZ'
s29 = 'CC-SPN5-EHZ'
s30 = 'CC-SEND-EHZ'

s100  = 'UW-APE-SHZ'
s110  = 'UW-APE-SHN'
s111  = 'UW-APE-SHE'
s200  = 'UW-CDF-SHZ'
s220  = 'UW-CDF-SHN'
s222  = 'UW-CDF-SHE'
s300  = 'UW-MUD-SHZ'
s330  = 'UW-MUD-SHN'
s333  = 'UW-MUD-SHE'

list_stations = [s100,s110,s111,s200,s220,s222,s300,s330,s333]

# list_stations = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,
#                  s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30]

# list_stations = [s2,s3,s4,s5,s6,s7,s10,s11] # make a list of all stations

if len(list_stations) == 0:
    netstacha = list_stations[0].split('-')
    sta = netstacha[1] # station
    dirpath_save = '/data/wsd03/data_manuela/MtStHelens/scatseisnet/o4_r4_q1_o5_r2_q2/2month/{}'.format(sta)

if len(list_stations) > 0:    
    dirpath_save = '/data/wsd03/data_manuela/MtStHelens/scatseisnet/o4_r4_q1_o5_r2_q2/2month/multistation/apr'

# Create directory to save the results
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
            endtime=obspy.UTCDateTime("2004-11-01T00:00:00.1"), #"2004-11-01T00:00.1"
        )
#         st.trim(obspy.UTCDateTime("1980-03-31T00:00"),obspy.UTCDateTime("1980-05-19T00:00:00.1"),pad=True, fill_value=0)

        st.detrend("linear")
        st.taper(0.05, type='hann')
        st.resample(50)
        st.merge(fill_value=0)

        st.filter(type="highpass", freq=1.0)
#         st.resample(100)
    #     st.plot(rasterized=True);
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
    stream.write(dirpath_save+"/scattering_stream_short_sta_list_apr2.mseed", format="MSEED")
    print(stream)

print('Loading data tooks {} minutes.'.format(round((time.time()-stime)/60,3)))

# create layers ======================================================================================================

segment_duration_seconds = 60.0
sampling_rate_hertz = stream[0].stats.sampling_rate
samples_per_segment = int(segment_duration_seconds * sampling_rate_hertz)
# the network will have 2 layers ------------------------------------------ ajustments hav to be done ---------------
bank_keyword_arguments = (
    {"octaves": 4, "resolution": 4, "quality": 1},
    {"octaves": 5, "resolution": 2, "quality": 2},
#     {"octaves": 5, "resolution": 2, "quality": 2},
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


# # if we already have a stream and a network cerated =================================================================

# reclen = 512
# chunksize = 100000 * reclen # Around 50 MB
# stream = obspy.Stream()

# with io.open("/data/wsd03/data_manuela/MtStHelens/scatseisnet/o4_r4_q1_o5_r2_q2/1month/multistation/apr/scattering_stream_short_sta_list_apr2.mseed", "rb") as fh:
#     while True:
#         with io.BytesIO() as buf:
#             c = fh.read(chunksize)
#             if not c:
#                 break
#             buf.write(c)
#             buf.seek(0, 0)
#             st = obspy.read(buf)
#         # Do something useful!
#         stream += st

# stream = stream.merge()
# print(stream)
        
# network = pd.read_pickle(dirpath_save+"/scattering_network.pickle")

# # preparation for network transform ----------------------------------------------------------------

# # Extract segment length (from any layer)
# segment_duration = network.bins / network.sampling_rate
# overlap = 0.5

# # Gather list for timestamps and segments
# timestamps = list()
# segments = list()

# # Collect data and timestamps
# for traces in stream.slide(segment_duration, segment_duration * overlap):
    
#     samples = [tr.stats.npts for tr in traces]
#     if any(samp != 3001 for samp in samples): # extend traces to lenght of 3001 by filling with 0
#         tr_startt = min([tr.stats.starttime for tr in traces]) # find earliest starttime
#         tr_endt = max([tr.stats.endtime for tr in traces]) # find latest endtime
#         traces.trim(tr_startt, tr_endt,pad=True, fill_value=0) # fill up with 0 if starttime is late or endtime is early
    
#     timestamps.append(mdates.num2date(traces[0].times(type="matplotlib")[0]))
#     segments.append(np.array([trace.data[:-1] for trace in traces]))

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

# Save the features in a pickle file --------------------------------------- ajustments hav to be done ---------------
np.savez(
    dirpath_save+"/scattering_coefficients.npz",
    order_1=scattering_coefficients[0],
    order_2=scattering_coefficients[1],
#     order_3=scattering_coefficients[2],
    times=timestamps,
    list_netstacha=list_stations
#     list_netstacha=list_netstacha,!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
)
print('Scattering transform and saving took {} hours.'.format(round((time.time()-st)/3600,3)))
print('Script run for {} hours.'.format(round((time.time()-stime)/3600,3)))