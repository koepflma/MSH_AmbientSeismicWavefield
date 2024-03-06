import numpy as np
import pandas as pd
import obspy
import obspy.signal.filter
import datetime
import scipy
import sys
import os
import time
import warnings
warnings.filterwarnings("ignore")

from obspy.clients.fdsn.client import Client 
client = Client('IRIS')

# Define all functions =======================================================================================

def preprocessing(year,jday, net, sta, cha):
                        
    s_time = obspy.UTCDateTime(year=year, julday=jday)
    e_time = s_time+24*3600

    try:
        # this stream will be used for RSAM and DSAR calculations
        st = client.get_waveforms(
            network=net,
            station=sta,
            location="*",
            channel=cha,
            starttime=s_time,
            endtime=e_time)

        st.detrend('linear')
        st.taper(max_percentage=None, max_length=5, type='hann') #max_length in sec
        
        # correct insrument response
        inv = client.get_stations(network=net, station=sta, location='*', channel=cha,
                                  starttime=s_time, endtime=e_time, level='response')
        pre_filt = [1e-3, 5e-2, 45, 50]
        water_level = 60
        
        for tr in st:
            tr.remove_response(inventory=inv, zero_mean=True,taper=True, taper_fraction=0.05,
                                      pre_filt=pre_filt, output="VEL", water_level=water_level,
                                      plot=False)

            # correct positive dip
            dip = inv.get_orientation(tr.id, datetime=tr.stats.starttime)['dip']
            if dip > 0:
                tr.data *= -1
#         st.merge(fill_value=0)
        print(':) year={}, jday={}, net={}, sta={}, cha={}'.format(year,jday, net, sta, cha))
    except:
        print('pass station {} day {}'.format(sta,jday))
    return(st)
    
def noise_analysis(data, datas, samp_rate, N, Nm):
    rms_list = []
    rmes_list = []
    pgv_list = []
    pga_list = []

    for i in np.arange(0,Nm,N): # start samples (sample, where next 10min starts)
        data_cut = data[i:i+N-1]

        rms = np.sqrt(np.mean(data_cut**2))
        rmes = np.sqrt(np.median(data_cut**2))
        pgv = max(abs(data_cut))

        data_acc = (data_cut.copy()[:-1] - data_cut.copy()[1:]) / (1/samp_rate)
        pga = max(abs(data_acc))

        rms_list.append(rms)
        rmes_list.append(rmes)
        pgv_list.append(pgv)
        pga_list.append(pga)
    datas.append(np.array(rms_list))
    datas.append(np.array(rmes_list))
    datas.append(np.array(pgv_list))
    datas.append(np.array(pga_list))
    return (datas)
    
def RSAM(data, samp_rate, datas, freq, Nm, N):
    filtered_data = obspy.signal.filter.bandpass(data, freq[0], freq[1], samp_rate)
    filtered_data = abs(filtered_data[:Nm])
    datas.append(filtered_data.reshape(-1,N).mean(axis=-1)*1.e9)
    return(datas)

def VSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N):
    # compute ratio between different velocities
    data -= np.mean(data) # detrend('mean')
    j = freqs_names.index('mf')
    mfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    mfd = abs(mfd[:Nm])
    mfd = mfd.reshape(-1,N).mean(axis=-1)
    j = freqs_names.index('hf')
    hfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    hfd = abs(hfd[:Nm])
    hfd = hfd.reshape(-1,N).mean(axis=-1)
    vsar = mfd/hfd
    datas.append(vsar)
    return(datas)

def lhVSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N):
    # compute ratio between different velocities
    data -= np.mean(data) # detrend('mean')
    j = freqs_names.index('rsam')
    mfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    mfd = abs(mfd[:Nm])
    mfd = mfd.reshape(-1,N).mean(axis=-1)
    j = freqs_names.index('hf')
    hfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    hfd = abs(hfd[:Nm])
    hfd = hfd.reshape(-1,N).mean(axis=-1)
    vsar = mfd/hfd
    datas.append(vsar)
    return(datas)

def DSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N):
    # compute dsar
    data = scipy.integrate.cumtrapz(data, dx=1./100, initial=0) # vel to disp
    data -= np.mean(data) # detrend('mean')
    j = freqs_names.index('mf')
    mfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    mfd = abs(mfd[:Nm])
    mfd = mfd.reshape(-1,N).mean(axis=-1)
    j = freqs_names.index('hf')
    hfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    hfd = abs(hfd[:Nm])
    hfd = hfd.reshape(-1,N).mean(axis=-1)
    dsar = mfd/hfd
    datas.append(dsar)
    return(datas)

def lDSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N):
    # compute dsar for low frequencies
    data = scipy.integrate.cumtrapz(data, dx=1./100, initial=0) # vel to disp
    data -= np.mean(data) # detrend('mean')
    j = freqs_names.index('rsam')
    lfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    lfd = abs(lfd[:Nm])
    lfd = lfd.reshape(-1,N).mean(axis=-1)
    j = freqs_names.index('mf')
    mfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    mfd = abs(mfd[:Nm])
    mfd = mfd.reshape(-1,N).mean(axis=-1)
    ldsar = lfd/mfd
    datas.append(ldsar)
    return(datas)

def lhDSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N):
    # compute dsar for low frequencies
    data = scipy.integrate.cumtrapz(data, dx=1./100, initial=0) # vel to disp
    data -= np.mean(data) # detrend('mean')
    j = freqs_names.index('rsam')
    lfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    lfd = abs(lfd[:Nm])
    lfd = lfd.reshape(-1,N).mean(axis=-1)
    j = freqs_names.index('hf')
    mfd = obspy.signal.filter.bandpass(data, freqs[j][0], freqs[j][1], samp_rate)
    mfd = abs(mfd[:Nm])
    mfd = mfd.reshape(-1,N).mean(axis=-1)
    ldsar = lfd/mfd
    datas.append(ldsar)
    return(datas)

def nDSAR(datas):
    dsar = datas[3]
    ndsar = dsar/scipy.stats.zscore(dsar)
    datas.append(ndsar)
    return(datas)    

def create_df(datas, ti, freqs_names, df):
    datas = np.array(datas)
    time = [(ti+j*600).datetime for j in range(datas.shape[1])]
    df_tr = pd.DataFrame(zip(*datas), columns=freqs_names+['rms','rmes','pgv','pga'], index=pd.Series(time))
    df = pd.concat([df, df_tr])
    return(df) 
    
# main function =======================================================================================
def freq_bands(jday, year, netstacha, freqs):   
    ''' 
    calculate and store power in 10 min long time windows for different frequency bands
    sensor measured ground velocity
    freqs: list contains min and max frequency in Hz
    dsar: float represents displacement (integration of)'''
    
    net = netstacha.split('-')[0]
    sta = netstacha.split('-')[1]
    cha = netstacha.split('-')[2]
    
    # path to save data
    file_path = '../output/RSAM_DSAR/data/{}/{}/'.format(year, sta)
    # bild path to save figures if it does not exist yet
    os.makedirs(file_path, exist_ok=True)
    file_name = '{}_{}.csv'.format(sta,jday)
        
    if os.path.isfile(file_path+file_name):
        print('file for {}-{} at {} already exist'.format(year,jday, netstacha))
        pass
    else:    
        start_time = time.time()
        freqs_names = ['rsam','mf','hf','dsar','ldsar', 'lhdsar', 'vsar', 'lhvsar']
        df = pd.DataFrame(columns=freqs_names)
        daysec = 24*3600

        st = preprocessing(year,jday, net, sta, cha)

        if len(st)>0: # if stream not empty
    #         st.resample(50)
            for tr in st:
    #         tr = st[0]
                datas = []
                data = tr.data
                samp_rate = tr.meta['sampling_rate']
                ti = tr.meta['starttime']
                # round start time to nearest 10 min increment
                tiday = obspy.UTCDateTime("{:d}-{:02d}-{:02d} 00:00:00".format(ti.year, ti.month, ti.day)) # date
                ti = tiday+int(np.round((ti-tiday)/600))*600 # nearest 10 min to starttime
                N = int(600*samp_rate)    # 10 minute windows in seconds
                Nm = int(N*np.floor(len(data)/N)) # np.floor rounds always to the smaller number
                # seconds per day (86400) * sampling rate (100) -> datapoints per day

                for freq, frequ_name in zip(freqs, freqs_names[:3]):
                    datas = RSAM(data, samp_rate, datas, freq, Nm, N) # get RSAM for different frequency bands

                datas = DSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N)
                datas = lDSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N)
                datas = lhDSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N)
                datas = VSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N)
                datas = lhVSAR(data, samp_rate, datas, freqs_names, freqs, Nm, N)
    #             datas = nDSAR(datas) # --> add ndsar in freqs_names

                datas = noise_analysis(data, datas, samp_rate, N, Nm)

                df = create_df(datas, ti, freqs_names, df)
                
            df.to_csv(file_path + file_name, index=True, index_label='time')
            print('One day tooks {} seconds.'.format(round(time.time()-start_time),3))
        else:
            print('empty stream station {} day {}'.format(sta,jday))
    return
