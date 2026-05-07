import os
import numpy as np
import netCDF4 as nc
import datetime as dt
import sys
from glob import glob
import matplotlib.pyplot as plt
import calendar
import numpy.ma as ma
import string
import matplotlib.cm as cm
from scipy import integrate
from scipy import interpolate
from datetime import datetime
from beam_model import eval_qs,eval_q,eval_theta_e
# limit of number of COSMIC-2 profiles to process
lim = 250000


# IMERG and COSMIC-2 directories
imerg_dir = '/neelin2020/COSMIC2/COSMIC2-IMERGfiles20221112095135/COSMIC2_IMERG/'
cosmic_dir = '/neelin2020/COSMIC2/data/'

# Sort through COSMIC-2 files
cosmic_glob = sorted(glob(cosmic_dir+'*'))
cosmic_match = np.empty(lim,dtype=object)
count = 0
for a,i in enumerate(cosmic_glob):
    z = len(sorted(glob(i+'/*nc')))
    if z==0:
        continue
    if count + z > lim:
        break
    cosmic_match[count:count+z] = sorted(glob(i+'/*nc'))
    count = count + z
cosmic_match = cosmic_match[:count]

# Detect bad days
date_format = '%Y-%m-%d_%H:%M:%S.%f'
dates = np.empty(count,dtype=object)


badinds = np.zeros(count,dtype=bool) 
for a,i in enumerate(cosmic_match):
    f = nc.Dataset(i,'r')
    if f.errstr == 'No wetPf2 netCDF file created':
        badinds[a] = True
        continue    
    dates[a] = datetime.strptime(f.date, date_format)
    f.close()

cosmic_match_good = np.copy(cosmic_match)[~badinds]
