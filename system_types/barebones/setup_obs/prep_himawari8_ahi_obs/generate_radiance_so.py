''' 
Script to construct text files containing Him8 AHI radiance data.
'''

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import datetime
import random
from netCDF4 import Dataset as ncopen



# USER SETTINGS
# --------------
hroi_d  = 100		        # HROI for EnKF's update to dynamical fields
hroi_hy = 100		        # HROI for EnKF's update to hydrometeor fields (deprecated)
thin = 6		        # Thinning distance in terms of AHI pixels 
sat = 'ahi_himawari8'           # Satellite name according to CRTM
channel=8                       # Sensor channel to use
lat0, lat1 = -20.0, 15.0        # Latitude boundaries of domain 
lon0, lon1 = 85, 130            # Longitude boundaries of domain
out_dir ='TARGET_DIRECTORY'     # Location to store text files
raw_dir ='RAW_FILES_DIRECTORY'  # Location containing raw Him8 netCDF files.

# First date of experiment 
itime = datetime.datetime.strptime('201705300100', '%Y%m%d%H%M')

# Final data of experiment
ftime = datetime.datetime.strptime('201705300200', '%Y%m%d%H%M')

# Time interval between successive scans (in minutes)
tint = 60



# ACTUAL DATA PROCESSING
# ------------------------

# time range
time = itime

# Iterate over time
while time <= ftime:

  # Reading in data
  fname = raw_dir+'/NC_H08_' + time.strftime('%Y%m%d_%H%M_') + '%s_FLDK.%05d_%05d.nc' % ('R21', 2401, 2401)

  f = ncopen( fname, 'r')

  Tbch8 = f['tbb_08'][:]
  lons  = f['longitude'][:]
  lats  = f['latitude'][:]
  Tbch14 = f['tbb_14'][:]

  f.close()

  ymax = lats.shape[0]
  xmax = lons.shape[0]

  # Set up vertical localization pressures depending on clear or empty
  pres = np.zeros([ymax, xmax])
  # Only dealing differentiating btwn tall clouds and lack thereof.
  Tbch14_flag = Tbch14 > 260
  pres[ Tbch14_flag ] = 400*100.0
  # If less than 12 deg C, ch2 should be seeing cloud
  pres[ np.invert( Tbch14_flag ) ] = 250*100.0
  Tb_done = Tbch8[:,:]

  # Now outputting in an appropriate format
  file_output = out_dir+'/Tb_d01_'+time.strftime('%Y%m%d%H%M')+'_so'
  f = open(file_output,'w')

  # Writting data
  for i in range(thin, xmax, thin):
    for j in range(thin, ymax, thin):
       error = 3.0
       Tbout = Tb_done[j,i]
       flag = True
       lon_flag = ( lons[i] < lon1 ) * ( lons[i] > lon0 )
       lat_flag = ( lats[j] < lat1 ) * ( lats[j] > lat0 )
       if lon_flag and lat_flag and Tbout>160 and Tbout<330:
          dataset = time.strftime('%Y%m%d%H%M')+"{0:>12}".format(sat)+"{0:>12}".format(channel)+"{0:12.3f}".format(float(lats[j]))
          dataset = dataset+"{0:12.3f}".format(float(lons[i]))+"{0:12.3f}".format(float(Tbout))+"{0:>12}".format(hroi_hy)+"{0:>12}".format(hroi_d)+"{0:12.3f}".format(error)+"{0:>12.3f}".format(pres[j,i])+'\n'
          f.write(dataset)


  print( "  Done processing :"+time.strftime('%Y%m%d%H%M'))
  time = time + datetime.timedelta(minutes = tint) 


