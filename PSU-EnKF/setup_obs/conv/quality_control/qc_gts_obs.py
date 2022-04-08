''' 
Script to QC GTS observations
'''

import numpy as np
import datetime
import sys
from funclib_qc_gts import *

# Date controls
date_st = datetime.datetime.strptime( sys.argv[1], "%Y%m%d%H%M" ) 
date_ed = datetime.datetime.strptime( sys.argv[2], "%Y%m%d%H%M" )
t_int = 60
date = date_st + datetime.timedelta( minutes=0 )

while date <= date_ed:
  # QC controls
  ctrl_params = {}
  ctrl_params['datestr'] = date.strftime( "%Y%m%d%H%M" )
  ctrl_params['date'] = date
  ctrl_params['obs_vlist'] = ['pres','wspd','wdir','z','t','dwpt', 'rh']
  ctrl_params['era5_vlist'] = ['u','v','q','t','z']
  ctrl_params['nSamp'] = 2000                 # number of monte carlo samples to interpolate observations
  
  # sanity check flags
  ctrl_params['check_interp'] = False
  
  # Prepare observation fname
  obs_gts = {}
  obs_gts['obs_raw_fname'] = date.strftime( "raw_obs/obs_gts_%Y-%m-%d_%H:%M:00.3DVAR" )
  obs_gts['era5_fname'] =    date.strftime( "era5_data/era5_reanalysis_%Y%m%d%H%M.nc" )
  obs_gts['obs_out_fname'] = date.strftime( "qc_obs/obs_gts_%Y-%m-%d_%H:%M:00.3DVAR" )
  
  # Running QC system
  read_obs_gts(  ctrl_params, obs_gts )
  era5_soundings( ctrl_params, obs_gts )
  compare_era5_obs( ctrl_params, obs_gts )
  write_obs_gts( ctrl_params, obs_gts )

  date += datetime.timedelta( minutes = t_int )
