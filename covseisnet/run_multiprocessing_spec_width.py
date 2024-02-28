"""
This script calls func_spectral_width_nsta and run a multi-
processing of the function spectral_width_nsta.
The code is adapted from covseisnet.
Output: npz array per day and preprocessing including
        UTC_time, frequencies, covariances, spectral_width, n_sta
"""
import numpy as np
from obspy.clients.fdsn import Client
import time
import argparse
import multiprocessing
from functools import partial
from func_spectral_width_nsta import spectral_width_nsta

# define parameters ======================================================================================= 

# stations as string 'network-station-channel'
s1  = 'UW-EDM-EHZ'
s2  = 'UW-SHW-EHZ'
s3  = 'UW-HSR-EHZ'

list_stations = [s1,s2,s3] # make a list of all stations

preprocessing_type = ['OBS'] # 'NoPreP', 'OBS', 'OBT', 'ST'
# NoPreP: No-preporcessing
# OBS: One-bit spectral withening
# OBT: One-bit temporal normalization
# ST: Spectral and temporal smoothing

#--> python run_multiprocessing_spec_width.py 2004 1 3 =========================================================================
parser = argparse.ArgumentParser(description='Calculate spectral width.')
parser.add_argument('year', type=int, help='Year of interest')
parser.add_argument('start_day', type=int, help='Julian day you want to start')
parser.add_argument('end_day', type=int, help='Julian day you want to end')

args = parser.parse_args()

year = args.year
jdays = range(args.start_day,args.end_day+1)

# multiprocessing =========================================================================
s_time = time.time()

p = multiprocessing.Pool(processes=24)
p.imap_unordered(partial(spectral_width_nsta,list_stations=list_stations,year=year, preprocessing_type=preprocessing_type), jdays)
p.close()
p.join()

print('Total calculation-time: {} s'.format(round(time.time()-s_time,3)))
