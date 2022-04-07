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
flag_compute_hovmoller = False

# Hovmoller latitude range to average over
lat_range=[ -10, 10]

# Path format for nature run BT ncfiles (will do date -> string conversion later)
path_fmt = { 'nature': 'nature_run/met8_bt_%Y%m%d%H%M.nc',
             'seviri': 'met8_ncfiles/seviri_%Y-%m-%d_%H%MUTC.nc' }

# Dictionary of channels to look at, including names of the channels from the 
# various datasets
ch_dict = { 'window' : { 'seviri': 'IR_108', 'nature': 'seviri_m08_ch010'},
            'wv'     : { 'seviri': 'WV_062', 'nature': 'seviri_m08_ch005'}
          }





# ======================================================================
# Compute Hovmoller averages 
# ----------------------------------------------------------------------

# Switch controlling whether to run Hovmoller computation
if flag_compute_hovmoller:
  # Predefine dictionary to hold Hovmoller average
  hov_data = {}
  
  # Estimating how many dates there are in the date range
  n_seconds = (date_ed - date_st).total_seconds()
  n_dates = int( n_seconds / (date_interval * 60 ) ) + 1
  
  
  # Dealing with SEVIRI and nature run data
  # ---------------------------------------
  for dset in ['seviri','nature']:
    date_nw = date_st
    dd = 0   # Index for handling date dimension (see later)
    while ( date_nw <= date_ed ):
     
      print( date_nw.strftime('Processing '+dset+' file on %Y-%m-%d %H:%M') )
  
      # Load file
      fname = date_nw.strftime( path_fmt[dset] )
      f = ncopen(fname, 'r')
    
      # Construct latitude mask
      latmask = np.array(f.variables['latitude'])
      latmask = ( latmask <= lat_range[1]) * ( latmask >= lat_range[0] )
      latmask = np.invert( latmask )
    
      # On-the-fly memory allcoation for first time iteration
      if ( date_nw == date_st ):
        shp2d = latmask.shape
        hov_data[dset] = {}
        hov_data[dset]['lon'] = np.zeros( [n_dates, shp2d[1]] ) + np.nan
        hov_data[dset]['date'] = []
        for ch in ch_dict.keys():
          hov_data[dset][ch] = np.zeros( [n_dates, shp2d[1]] ) + np.nan
      # --- End of on-the-fly memory allocation
    
      # Compute number of valid pixels
      n_pix = np.sum( np.invert(latmask), axis=0)
    
      # For each channel, compute average over latitude range
      for ch in ch_dict.keys():
        # Construct masked array of bts
        bt = np.ma.array( np.array( f.variables[ ch_dict[ch][dset] ] ), 
                          mask = latmask )
        # Take average over non-masked values
        bt = bt.sum( axis=0 ) / n_pix
        # Save hovmoller value of bt
        hov_data[dset][ch][dd] = bt*1.
    
      # Compute average longitude
      lon = np.ma.array( np.array( f.variables[ 'longitude' ] ), mask = latmask )
      lon = lon.sum( axis=0) / n_pix
      hov_data[dset]['lon'][dd] = lon*1.
    
      # Store date
      hov_data[dset]['date'].append( date_nw )

      # Close file
      f.close()
    
    
      # Increment date and date index
      date_nw += timedelta( minutes=date_interval)
      dd += 1
    
    # --- End of loop over BT netcdf files
  # --- End of loop over different datasets

  print("Finished computing Hovmoller data")

  # Storing data
  f = open( 'hov_data.pkl','wb')
  pickle.dump( hov_data, f )
  f.close()

# -- End of Hovmoller computation if statement



# ====================================================================
# Making Hovmoller plots 
# --------------------------------------------------------------------

# Load pickle data
f = open( 'hov_data.pkl', 'rb' )
hov_data = pickle.load( f )


# Generate figure object
fig, axs = plt.subplots( ncols = 2, nrows = 1, figsize=(6,6) )

crange = np.linspace(220,240,11)
cmap = 'jet_r'

# Plot seviri stuff
print( len(hov_data['seviri']['date']) )
print( np.mean( hov_data['seviri']['lon'], axis=0 ).shape)
print( hov_data['seviri']['window'].shape)
axs[0].contourf( np.mean( hov_data['seviri']['lon'], axis=0 ),
                 hov_data['seviri']['date'] ,
                 hov_data['seviri']['window'], crange, cmap = cmap )
axs[0].set_title('SEVIRI Window-BT')

# Plot nature run stuff
axs[1].contourf( np.mean( hov_data['nature']['lon'], axis=0 ),
                 hov_data['seviri']['date'] ,
                 hov_data['nature']['window'], crange, cmap = cmap )
axs[1].set_title('Nature Window-BT')

plt.savefig('hovmoller_comparison.png' )
plt.close()

