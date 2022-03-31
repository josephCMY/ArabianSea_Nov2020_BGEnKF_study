''' 
Sanity checking the QC functions.
'''

import pickle
import numpy as np
from netCDF4 import Dataset as ncopen
import datetime
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import wrf
import time
from scipy.interpolate import interp1d, griddata
import sys
import math
from funclib_qc_gts import *


# User controls
ctrl_params = {}
ctrl_params['datestr'] = "201705310000"
ctrl_params['date'] = datetime.datetime.strptime( ctrl_params['datestr'] , "%Y%m%d%H%M")
ctrl_params['n_ens'] =50 
ctrl_params['obs_vlist'] = ['pres','wspd','wdir','hgt','t','dwpt', 'rh']
ctrl_params['era5_vlist'] = ['u','v','q','t','z']
ctrl_params['nSamp'] = 1000                 # number of monte carlo samples to interpolate observations

# sanity check flags
ctrl_params['check_interp'] = False
#ctrl_params['check_location'] = False

# Prepare observation fname
obs_gts = {}
obs_gts['obs_raw_fname'] = "gts_raw/obs_gts_2017-05-31_00:00:00.3DVAR"
obs_gts['era5_fname'] = "era5_data/era5_reanalysis_2017-05-31_00UTC.nc"
obs_gts['obs_out_fname'] = "obs_gts_2017-05-31_00:00:00.3DVAR"

'''
TEST 1: READ WRITE TEST
read_gts_obs( ctrl_params, obs_gts )
write_gts_obs( ctrl_params, obs_gts )
'''

'''
TEST 2: CHECK VERTICAL INTERPOLATION
read_gts_obs( ctrl_params, obs_gts )
ctrl_params['check_interp'] = True
era_soundings( ctrl_params, obs_gts )
'''






