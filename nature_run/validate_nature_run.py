'''
Python script to validate nature run against Meteosat-8 IR observations 

Stuff to check:
    1) Hovmoller diagram of Window-BT -- just to check for MJO initiation

Two-phase approach:
    1) Compute the Hovmoller averages first, output a pickle file
    2) Load pickle file and plot the Hovmoller diagram.
'''

import numpy as np
from netCDF4 import Dataset as ncopen
import pickle
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta


# =====================================================================
# USER INPUTS
# ---------------------------------------------------------------------

# Date controls
date_st = datetime.strptime('202011070000', '%Y%m%d%H%M')
date_ed = datetime.strptime('202011110300', '%Y%m%d%H%M')
date_interval = 60 # In minutes

# Switch to determine if Hovmoller needs calculating
flag_compute_hovmoller = True

# Hovmoller latitude range to average over
lat_range=[-10,10]

# Path format for nature run BT ncfiles (will do date -> string conversion later)
nature_run_format = 'nature_run/met8_bt_%Y%m%d%H%M.nc'
seviri_obs_format = 'met8_ncfiles/seviri_%Y-%m-%d_%H%MUTC.nc'

# Dictionary of channels to look at, including names of the channels from the 
# various datasets
ch_dict = { 'window' : { 'seviri': 'IR_108', 'nature': 'seviri_m08_ch010'},
            'wv'     : { 'seviri': 'WV_062', 'nature': 'seviri_m08_ch005'}
          }





# ======================================================================
# Compute Hovmoller averages 
# ----------------------------------------------------------------------

# Predefine dictionary to hold Hovmoller average
hov_data = {}

# Estimating how many dates there are in the date range
n_seconds = (date_ed - date_st).total_seconds()
n_dates = int( n_seconds / (date_interval * 60 ) ) + 1


# Dealing with SEVIRI data
# ------------------------
date_nw = date_st
dd = 0   # Index for handling date dimension (see later)
while ( date_nw <= date_ed ):
 
  print( date_nw.strftime('Processing SEVIRI file on %Y-%m-%d %H:%M') )
  # Load SEVIRI file
  fname = date_nw.strftime( seviri_obs_format )
  f = ncopen(fname, 'r')

  # Construct latitude mask
  latmask = f.variables['latitude']
  latmask = ( latmask <= lat_range[1]) * ( latmask >= lat_range[0] )
  latmask = np.invert( latmask )

  # On-the-fly memory allcoation for first time iteration
  if ( date_nw == date_st ):
    shp2d = latmask.shape
    hov_data['seviri'] = {}
    hov_data['seviri']['lon'] = np.zeros( [n_dates, shp2d[1]] ) + np.nan
    hov_data['seviri']['date'] = np.zeros( [n_dates, shp2d[1] ) + date_nw
    for ch in ch_list.keys():
      hov_data['seviri'][ch] = np.zeros( [n_dates, shp2d[1]] ) + np.nan
  # --- End of on-the-fly memory allocation

  # Compute number of valid pixels
  n_pix = np.sum( np.invert(latmask), axis=0)

  # For each channel, compute average over latitude range
  for ch in ch_dict.keys():
    # Construct masked array of bts
    bt = np.ma.array( f.variables[ ch_dict[ch]['seviri'] ], mask = latmask )
    # Take average over non-masked values
    bt = bt.sum( axis=0 ) / n_pix
    # Save hovmoller value of bt
    hov_data['seviri'][ch][dd] = bt*1.

  # Compute average longitude
  lon = np.ma.array( f.variables[ 'longitude' ], mask = latmask )
  lon = lon.sum( axis=0) / n_pix
  hov_data['seviri']['lon'][dd] = lon*1.

  # Store date
  hov_data['seviri']['date'][dd] = date_nw


  # Increment date and date index
  date_nw += timedelta( minutes=date_interval)
  dd += 1

# --- End of loop over SEVIRI files

  


# Dealing with nature run data
# -----------------------------
date_nw = date_st
dd = 0   # Index for handling date dimension (see later)
while ( date_nw <= date_ed ):
 
  print( date_nw.strftime('Processing nature run on %Y-%m-%d %H:%M') )

  # Load nature run file
  fname = date_nw.strftime( nature_run_format )
  f = ncopen(fname, 'r')

  # Construct latitude mask
  latmask = f.variables['latitude']
  latmask = ( latmask <= lat_range[1]) * ( latmask >= lat_range[0] )
  latmask = np.invert( latmask )

  # On-the-fly memory allcoation for first time iteration
  if ( date_nw == date_st ):
    shp2d = latmask.shape
    hov_data['nature'] = {}
    hov_data['nature']['lon'] = np.zeros( [n_dates, shp2d[1]] ) + np.nan
    hov_data['nature']['date'] = np.zeros( [n_dates, shp2d[1] ) + date_nw
    for ch in ch_list.keys():
      hov_data['nature'][ch] = np.zeros( [n_dates, shp2d[1]] ) + np.nan
  # --- End of on-the-fly memory allocation

  # Compute number of valid pixels
  n_pix = np.sum( np.invert(latmask), axis=0)

  # For each channel, compute average over latitude range
  for ch in ch_dict.keys():
    # Construct masked array of bts
    bt = np.ma.array( f.variables[ ch_dict[ch]['nature'] ], mask = latmask )
    # Take average over non-masked values
    bt = bt.sum( axis=0 ) / n_pix
    # Save hovmoller value of bt
    hov_data['nature'][ch][dd] = bt*1.

  # Compute average longitude
  lon = np.ma.array( f.variables[ 'longitude' ], mask = latmask )
  lon = lon.sum( axis=0) / n_pix
  hov_data['nature']['lon'][dd] = lon*1.

  # Store date
  hov_data['nature']['date'][dd] = date_nw


  # Increment date and date index
  date_nw += timedelta( minutes=date_interval)
  dd += 1

# --- End of loop over nature run files


# Storing Hovmoller data
f = open( 'hov_data.pkl', 'wb')
pickle.dump( hov_data, f )
f.close()



print("Finished computing Hovmoller data")
