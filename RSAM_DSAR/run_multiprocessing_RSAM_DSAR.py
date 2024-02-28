import time
import argparse
import warnings
warnings.filterwarnings("ignore")

import multiprocessing
from functools import partial

from func_RSAM_DSAR import freq_bands

# define parameters =======================================================================================

# stations as string 'network-station-channel'
s1  = 'UW-EDM-EHZ'
s2  = 'UW-SHW-EHZ'
s3  = 'UW-HSR-EHZ'

list_stations = [s1,s2,s3] # make a list of all stations

freqs = [[2,5], [4.5,8], [8,16]] # frequenciy bands Hz

#--> python RSAM_DSAR.py 2004 2 3 =========================================================================

parser = argparse.ArgumentParser(description='Calculate different frequency bands of RSMA, DSAR, RMS and more.')
parser.add_argument('year', type=int, help='Year of interest')
parser.add_argument('start_day', type=int, help='Julian day you want to start')
parser.add_argument('end_day', type=int, help='Julian day you want to end')

args = parser.parse_args()

year = args.year
jdays = ['{:03d}'.format(jday) for jday in range(args.start_day,args.end_day+1)]

# multiprocessing =========================================================================

for netstacha in list_stations:
    print('Station {}'.format(netstacha))
    stime = time.time()
    p = multiprocessing.Pool(processes=24)
    p.imap_unordered(partial(freq_bands,year=year, netstacha=netstacha, freqs=freqs), jdays)
    p.close()
    p.join()
    print('Calculation tooks {} seconds.'.format(round(time.time()-stime),3))

# for netstacha in list_stations:
#     for jday in jdays:
#         freq_bands(jday, year, netstacha, freqs)
