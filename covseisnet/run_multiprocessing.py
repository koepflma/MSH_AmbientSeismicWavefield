# import packages
import numpy as np
from obspy.clients.fdsn import Client
import time
import argparse
import multiprocessing
from functools import partial
from func_spectral_width_nsta import spectral_width_nsta
# from func_spectral_width import spectral_width

# define some parameters------------------------------------------------------------------------------------------------

client = Client("IRIS")  

# stations as string 'network-station-channel'
s1  = 'UW-EDM-EHZ'
s2  = 'UW-SHW-EHZ'
s3  = 'UW-HSR-EHZ'
s4  = 'UW-SOS-EHZ'
s5  = 'UW-JUN-EHZ'
s6  = 'UW-ELK-EHZ'
s7  = 'UW-TDL-EHZ'
s8  = 'UW-SUG-EHZ'
s9  = 'UW-YEL-EHZ'
s10 = 'UW-FL2-EHZ'
s11 = 'UW-CDF-EHZ' #-?H?

s12 = 'UW-SEP-EHZ' #-?H?
s13 = 'CC-SEP-EHZ' #-?H?
# s14 = 'UW-STD-EHZ'
s15 = 'CC-STD-BHZ'

s16 = 'CC-VALT-BHZ' #-BH?
s17 = 'CC-JRO-BHZ'
s18 = 'CC-HOA-BHZ' #-BH?
s19 = 'CC-LOO-BHZ' #-BH?
s20 = 'CC-USFR-BHZ' #-BH?
s21 = 'CC-NED-EHZ'
s22 = 'CC-REM-BHZ' #-BH?
s23 = 'CC-SWFL-BHZ' #-BH?
s24 = 'CC-SFW2-BHZ' #-BH?
s25 = 'CC-MIDE-EHZ'
s26 = 'CC-MIBL-EHZ'
s27 = 'CC-BLIS-EHZ'
s28 = 'CC-RAFT-EHZ'
s29 = 'CC-SPN5-EHZ'
s30 = 'CC-SEND-EHZ'

list_stations = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,
                 s12,s13,s15]#,
                 #s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30] # make a list of all stations

# list_stations = [s2,s3,s4,s5,s6,s7,s10,s11,s12] # make a list of all stations
# list_stations = [s1,s2,s3,s4,s5,s6,s7,s10,s11, s12]

# define time parameters
# year = 2003 # year
# jdays = range(1,365+1)

preprocessing_type = ['OBS'] # 'NoPreP', 'OBT', 'OBS', 'ST'


parser = argparse.ArgumentParser(description='Calculate spectral width.')
parser.add_argument('year', type=int, help='Year of interest')
parser.add_argument('start_day', type=int, help='Julian day you want to start')
parser.add_argument('end_day', type=int, help='Julian day you want to end')

args = parser.parse_args()

year = args.year
jdays = range(args.start_day,args.end_day+1)

#--> python run_multiprocessing.py 2004 1 3 ['NoPreP', 'OBT', 'OBS', 'ST']

# import packages for multi processing and the main function-------------------------------------------------------------


s_time = time.time()

p = multiprocessing.Pool(processes=24)
p.imap_unordered(partial(spectral_width_nsta,client=client,list_stations=list_stations,year=year, preprocessing_type=preprocessing_type), jdays)
p.close()
p.join()

print('Total calculation-time: {} s'.format(round(time.time()-s_time,3)))