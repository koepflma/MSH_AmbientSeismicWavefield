import covseisnet as csn
import numpy as np
import obspy
from obspy import UTCDateTime
import time
from obspy.clients.fdsn import Client
import os
from itertools import compress
import warnings
warnings.filterwarnings("ignore")

''' spectral_width is a function to calculate one day of the spectral width
input: jday: integer, interger correspond to julian day
        client: string, waveform client
        list_stations: list with strings, string is 'network-station-channel'
        year: integer, year
        preprocessing_type: list with strings, possibilities: 'no_preprocessing', 'OBT', 'OBS', 'spectral_temporal'
output: save UTCDateTime, frequency, spectral width AND number of stations as .npz (per day) for the selected preprocessing(s)'''

def spectral_width_nsta(jday, client, list_stations, year, preprocessing_type):
    start_time = time.time()
    
    # define time parameter--------------------------------------------------------------------------------------------------

    hour = 0 # start hour
    signal_duration_sec = 24 * 3600 # hour * sec

    window_duration_sec = 30 # Fourier calcualation window in seconds (depends on array opening)
    average = 60 # number of windows to estimate the sample covariance (bigger than number of stations)

    # check which files already exist---------------------------------------------------------------------------------------

    file_path = '/data/whd02/data_manuela/MtStHelens/covariance/{}/{:03d}/'.format(year,jday)
#     file_path = '/data/wsd03/data_manuela/MtStHelens/covariance/tremor/{}/{:03d}/'.format(year,jday)

    file_list = [file_path+'{}_{}_{:03d}_wd{}_av{}.npz'.format(preprocessing,year,jday,int(window_duration_sec), average) 
                 for preprocessing in preprocessing_type]

    if all(list(map(os.path.isfile,file_list))):
        print('all files for {}-{} already exist'.format(year,jday))
        pass

    else:
        boolean_list = ~np.array(list(map(os.path.isfile,file_list))) # inverse boolean list of existing files False if exist
        missing_file_list = list(compress(file_list, boolean_list)) # make list of missing files
        preprocessing_type = [] # overwrite list

        for missing_file in missing_file_list: # loop over missing files
            preprocessing_type.append(missing_file.split('/')[-1].split('_')[0]) # add missing preprocessing to list
    
        if not os.path.exists(file_path): # create path if not exist
            os.makedirs(file_path)

        # read data and synchronize data-------------------------------------------------------------------------------------

        stream = csn.arraystream.ArrayStream()

        t = UTCDateTime(year=year, julday=jday, hour=hour)
        for xx in list_stations:
            netstacha = xx.split('-')
            net = netstacha[0]
            sta = netstacha[1]
            cha = netstacha[2]
            try:
                st = client.get_waveforms(net, sta, '', cha, t, t + signal_duration_sec)

                npts = 0
                for tr in st:
                    npts += tr.stats.npts/tr.stats.sampling_rate
                if npts/(signal_duration_sec) < 0.75: # if stream contains less than 75% of data
                    print('{} {} less than 75% data'.format(netstacha, jday))
                    pass
                else:
                    
                    # correct insrument response
                    inv = obspy.read_inventory('/auto/pnwstore1-wd11/PNWStationXML/{}/{}.{}.xml'.format(net,net,sta))
                    pre_filt = [1e-3, 5e-2, 45, 50]
                    water_level = 60
                    
                    # correct positive dip
                    for tr in st:
                        dip = inv.get_orientation(tr.id, datetime=tr.stats.starttime)['dip']
                        if dip > 0:
                            tr.data *= -1
                            
                    st.merge(fill_value=0)
                    st[0].remove_response(inventory=inv, zero_mean=True,taper=True, taper_fraction=0.05,
                                              pre_filt=pre_filt, output="VEL", water_level=water_level,
                                              plot=False)

                    st.detrend('linear')
                    st.taper(max_percentage=None,max_length=300, type='hann') #max_length in sec
                            
            #         # downsample data to 25 Hz
            #         st[0].resample(25)
                    
                    stream.append(st[0])
                print(netstacha, jday)
                
            except:
                print('no data for {} at {}'.format(netstacha, jday))
                pass

        if len(stream) == 0: # if no station was active or not enough data
            samplingrate = 100
            freq_len = window_duration_sec*2*samplingrate
            time_len = window_duration_sec*0.25*average
            UTC_times = np.array([t+sec for sec in np.arange(0,signal_duration_sec+1,time_len)])
            UTC_times = np.delete(UTC_times,[-3,-2]) # the third and second last time steps are not included (don't know why)
            UTC_times[-1] = UTC_times[-1]-0.01 # last time step is 23, 59, 59, 990000 not 24, 0, 0, 0
            
            # 24 h, wd 60, av 60
#             frequencies = np.linspace(0,100,12000) # start, stop, number of steps
#             spectral_width = np.empty((94, 11999),)
#             spectral_width[:] = np.nan
#             covariances = np.empty((94, 11999, 6, 6),)
#             covariances[:] = np.nan

            frequencies = np.linspace(0,samplingrate,freq_len) # start, stop, number of steps
            spectral_width = np.empty((time_len-1, freq_len-1),)
            spectral_width[:] = np.nan
            covariances = np.empty((time_len-1, freq_len-1, 1, 1),)
            covariances[:] = np.nan

            for preprocessing in preprocessing_type:
                file_name = '{}_{}_{:03d}_wd{}_av{}'.format(preprocessing,year,jday,int(window_duration_sec), average)

                np.savez_compressed(file_path + file_name,UTC_times=UTC_times, frequencies=frequencies, covariances=covariances, spectral_width=spectral_width, n_sta=0)

        else:

            # synchronize data
            stream = stream.synchronize(start=t, duration_sec=signal_duration_sec, method="linear")

            # calculate & save spectral width for different pre-processing---------------------------------------------------------

            if 'NoPreP' in preprocessing_type: 
                    times, frequencies, covariances = csn.covariancematrix.calculate(
                stream, window_duration_sec, average)

                    UTC_times = np.array([t+sec for sec in times]) # convert seconds to UTCDateTime

                    spectral_width = covariances.coherence(kind="spectral_width")

                    file_name = 'NoPreP_{}_{:03d}_wd{}_av{}'.format(year,jday,int(window_duration_sec), average)

                    np.savez_compressed(file_path + file_name,UTC_times=UTC_times, frequencies=frequencies, covariances=covariances, spectral_width=spectral_width, n_sta=len(stream))

            if 'OBS' in preprocessing_type:
                    stream.preprocess()

                    times, frequencies, covariances = csn.covariancematrix.calculate(
                stream, window_duration_sec, average)

                    UTC_times = np.array([t+sec for sec in times]) # convert seconds to UTCDateTime

                    spectral_width = covariances.coherence(kind="spectral_width")

                    file_name = 'OBS_{}_{:03d}_wd{}_av{}'.format(year,jday,int(window_duration_sec), average)

                    np.savez_compressed(file_path + file_name,UTC_times=UTC_times, frequencies=frequencies, covariances=covariances, spectral_width=spectral_width, n_sta=len(stream))    

            if 'OBT' in preprocessing_type:
                    stream.preprocess(domain="temporal", method="onebit")

                    times, frequencies, covariances = csn.covariancematrix.calculate(
                stream, window_duration_sec, average)

                    UTC_times = np.array([t+sec for sec in times]) # convert seconds to UTCDateTime

                    spectral_width = covariances.coherence(kind="spectral_width")

                    file_name = 'OBT_{}_{:03d}_wd{}_av{}'.format(year,jday,int(window_duration_sec), average)

                    np.savez_compressed(file_path + file_name,UTC_times=UTC_times, frequencies=frequencies, covariances=covariances, spectral_width=spectral_width, n_sta=len(stream))

            if 'ST' in preprocessing_type:
                    stream.preprocess(domain="spectral", method="smooth")
                    stream.preprocess(domain="temporal", method="smooth")

                    times, frequencies, covariances = csn.covariancematrix.calculate(
                stream, window_duration_sec, average)

                    UTC_times = np.array([t+sec for sec in times]) # convert seconds to UTCDateTime

                    spectral_width = covariances.coherence(kind="spectral_width")

                    file_name = 'ST_{}_{:03d}_wd{}_av{}'.format(year,jday,int(window_duration_sec), average)

                    np.savez_compressed(file_path + file_name,UTC_times=UTC_times, frequencies=frequencies, covariances=covariances, spectral_width=spectral_width, n_sta=len(stream))

    end_time = time.time()
    print('Day {} calculation-time: {} s'.format(jday,round(end_time-start_time,3)))
    return()